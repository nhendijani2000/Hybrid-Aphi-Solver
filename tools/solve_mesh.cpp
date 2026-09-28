// Assemble, factorize and solve one input file, and report what happened.
//
// `docs/SOLVER_PLAN.md` step 5 of §9. The test suite cannot do this: factorizing
// a real mesh takes tens of seconds at Stage 1's scalar speed, which is fine for
// a tool and not fine for 1700 checks. Every step from here -- equilibration,
// static pivoting, refinement, the DC milestone -- needs somewhere to watch those
// numbers move, and this is it.
//
//     solve_mesh examples/loop_sweep.aphi
//
// It reports the residual on the ORIGINAL system, which is the only number that
// says whether the answer is any good, plus the pivot range, which is what says
// whether it was close to not being.

#include <algorithm>
#include <chrono>
#include <cstring>
#include <filesystem>
#include <system_error>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "aphi_solver/assembly.hpp"
#include "aphi_solver/factorization.hpp"
#include "aphi_solver/gmsh_reader.hpp"
#include "aphi_solver/input_file.hpp"
#include "aphi_solver/postprocess.hpp"
#include "aphi_solver/version.hpp"

using namespace aphi_solver;
using Clock = std::chrono::steady_clock;
using Complex = std::complex<double>;

namespace {

double ms(Clock::time_point a, Clock::time_point b) {
    return std::chrono::duration<double, std::milli>(b - a).count();
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2 || argc > 3) {
        std::cerr << "A-Phi solver " << kVersion << " -- assemble, factorize, solve one file\n\n"
                  << "usage: solve_mesh <input-file> [natural|rcm|amd]\n\n"
                  << "The ordering defaults to amd. The conditioning comes from the file's\n"
                  << "[solver] section; only the symmetric ones can be solved yet, since the\n"
                  << "unsymmetric LU is step 8 of docs/SOLVER_PLAN.md.\n\n"
                  << "try: solve_mesh examples/loop_sweep.aphi\n";
        return 2;
    }

    Ordering ordering = Ordering::ApproximateMinimumDegree;
    if (argc == 3 && !ordering_from_keyword(argv[2], ordering)) {
        std::cerr << "unknown ordering '" << argv[2] << "'; expected natural, rcm or amd\n";
        return 2;
    }

    try {
        const ParseResult parsed = parse_input_file(argv[1]);
        const Problem& p = parsed.problem;

        if (!conditioning_needs_ac(p.conditioning)) {
            std::cerr << argv[1] << ": conditioning = " << conditioning_keyword(p.conditioning)
                      << " gives an UNSYMMETRIC matrix, and the LU path is step 8 of\n"
                      << "docs/SOLVER_PLAN.md. Use row_scaled or scaled_phi for now -- both need\n"
                      << "a non-zero frequency.\n";
            return 1;
        }

        // Make the output directory before anything long-running, so a missing
        // or unwritable path fails in a second rather than after the solve.
        if (!p.output_dir.empty()) {
            std::error_code ec;
            std::filesystem::create_directories(p.output_dir, ec);
            if (ec && !std::filesystem::is_directory(p.output_dir)) {
                std::cerr << "output: could not create \"" << p.output_dir
                          << "\": " << ec.message() << "\n";
                return 1;
            }
            std::cout << "output dir    " << p.output_dir << "\n";
        }

        Mesh mesh = read_gmsh_msh(p.mesh_file);
        scale_mesh_to_metres(mesh, p.length_unit);
        const BoundProblem bound = bind_to_mesh(p, mesh);
        const DofMap dofs = build_dof_map(bound, mesh);
        const SparsityPattern pattern =
            build_sparsity(dofs, bound, mesh, SparsityStorage::UpperTriangle);

        std::cout << "A-Phi solver " << kVersion << " -- solve_mesh\n\n";
        std::cout << "input         " << argv[1] << "\n";
        std::cout << "mesh          " << p.mesh_file << "   " << mesh.num_tets() << " tets\n";
        std::cout << "conditioning  " << conditioning_keyword(p.conditioning)
                  << "   (symmetric, upper triangle stored)\n";
        std::cout << "unknowns      " << dofs.num_total << "   " << pattern.nnz()
                  << " stored nonzeros\n";
        for (const std::string& w : bound.warnings) std::cout << "warning       " << w << "\n";

        // Symbolic, once: it does not depend on frequency.
        auto t0 = Clock::now();
        const SolverAnalysis analysis = analyze(pattern, ordering);
        auto t1 = Clock::now();
        std::cout << "\nordering      " << ordering_keyword(ordering) << "   nnz(L) = "
                  << analysis.predicted_nnz << "   "
                  << analysis.predicted_nnz * 20.0 / 1048576.0 << " MB   [" << ms(t0, t1)
                  << " ms]\n\n";

        SymmetricSystem system = make_symmetric_system(pattern, dofs.num_total);
        const std::vector<double> frequencies =
            p.type == AnalysisType::DC ? std::vector<double>{0.0} : p.frequencies;

        // Each frequency is solved TWICE -- with equilibration and without -- so
        // the two residuals sit side by side. A number only ever seen scaled
        // cannot say whether the scaling helped.
        std::cout << std::right << "  " << std::setw(12) << "frequency" << std::setw(14)
                  << "resid plain" << std::setw(14) << "resid equil" << std::setw(14) << "backward err" << std::setw(12)
                  << "d max/min" << std::setw(13) << "min |D|" << std::setw(13) << "max |D|"
                  << std::setw(12) << "max |L|" << std::setw(11) << "factor ms" << "\n";

        int failures = 0;
        for (double f : frequencies) {
            const double omega = 2.0 * 3.14159265358979323846 * f;
            refill(system, bound, mesh, dofs, pattern, omega, p.conditioning);

            SolveOptions plain;
            plain.ordering = ordering;
            plain.equilibration_iterations = 0;
            SolveOptions equil;
            equil.ordering = ordering;
            equil.equilibration_iterations = 10;

            std::vector<Complex> x_plain;
            std::vector<Complex> x_equil;
            SolveReport rp;
            SolveReport re;
            const bool ok_plain = solve_symmetric(system.matrix, system.rhs, x_plain, rp, plain);
            const bool ok_equil = solve_symmetric(system.matrix, system.rhs, x_equil, re, equil);

            std::ostringstream label;
            label << std::setprecision(4) << f << " Hz";
            std::cout << "  " << std::setw(12) << label.str();
            if (!ok_plain || !ok_equil) {
                std::cout << "   FAILED at column "
                          << (ok_plain ? re.columns_done : rp.columns_done) << " of "
                          << dofs.num_total << "\n";
                ++failures;
                continue;
            }
            std::cout << std::setw(14) << rp.residual << std::setw(14) << re.residual
                      << std::setw(14) << re.backward_error
                      << std::setw(12) << re.equilibration_max / re.equilibration_min
                      << std::setw(13) << re.smallest_pivot << std::setw(13) << re.largest_pivot
                      << std::setw(12) << re.largest_multiplier << std::setw(11) << re.factorize_ms
                      << "\n";
            if (!(re.backward_error < 1e-12)) ++failures;

            // The potential, written for every frequency. One file when there is
            // only one solve, numbered otherwise, so a sweep does not silently
            // overwrite itself.
            Solution sol;
            sol.x = x_equil;
            sol.conditioning = p.conditioning;
            sol.omega = omega;
            const NodalPotential nodal = potential_at_nodes(mesh, bound, dofs, sol);
            std::string out_path = "potential.out";
            if (frequencies.size() > 1) {
                out_path = "potential_" + std::to_string(&f - frequencies.data()) + ".out";
            }
            if (!p.output_dir.empty()) out_path = p.output_dir + "/" + out_path;
            const WriteStats ws = write_potential(out_path, nodal, sol);
            const auto fields_started = std::chrono::steady_clock::now();
            const FieldOutput fields = compute_fields(mesh, bound, dofs, sol);
            const double fields_ms =
                std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() -
                                                          fields_started)
                    .count();
            const std::string stem = out_path.substr(0, out_path.size() - 4);
            std::string vtk_path = stem + ".vtk";
            const WriteStats vs =
                write_vtk(vtk_path, mesh, nodal, sol, &fields, FieldSet::Potential);
            std::cout << "                wrote " << out_path << "   " << ws.nodes << " nodes, "
                      << ws.bytes / 1024 << " KB, " << ws.milliseconds << " ms\n";
            std::cout << "                fields A,B,H,E at " << fields.num_nodes()
                      << " nodes in " << fields_ms << " ms";
            if (fields.num_interface_nodes > 0) {
                std::cout << "   " << fields.num_interface_nodes << " on a material interface";
            }
            if (fields.num_orphan_vertices > 0) {
                std::cout << "   (" << fields.num_orphan_vertices << " orphan vertices)";
            }
            std::cout << "\n";
            std::cout << "                wrote " << vtk_path << "   " << vs.bytes / 1024
                      << " KB, " << vs.milliseconds << " ms   (Phi)\n";

            // One pair of files per field. Opening B_field.vtk gives a source
            // with B on it and nothing else to pick through. There is no
            // combined file: to see two fields together, open two of these in
            // the same ParaView session -- they share a mesh, so the views line
            // up, and nothing has to be written twice to allow it.
            const std::string dir = stem.substr(0, stem.find_last_of("/\\") + 1);
            const std::string suffix = stem.substr(dir.size());
            const FieldSet sets[] = {FieldSet::A, FieldSet::B, FieldSet::H, FieldSet::E,
                                     FieldSet::J};
            for (const FieldSet fset : sets) {
                const std::string base =
                    dir + field_set_name(fset) + "_field" +
                    (suffix == "potential" ? "" : suffix.substr(std::strlen("potential")));
                const WriteStats ts = write_solution(base + ".out", nodal, fields, sol, fset);
                const WriteStats vv =
                    write_vtk(base + ".vtk", mesh, nodal, sol, &fields, fset);
                std::cout << "                wrote " << base << ".out / .vtk   "
                          << ts.bytes / 1024 << " + " << vv.bytes / 1024 << " KB, "
                          << (ts.milliseconds + vv.milliseconds) << " ms\n";
            }
        }

        std::cout << "\nNothing is EXTRACTED from these solutions yet: no currents, voltages, R or\n"
                  << "L come out of this. The residual says the linear system was solved; whether\n"
                  << "the physics is right is the DC milestone, step 9 of docs/SOLVER_PLAN.md.\n";
        return failures == 0 ? 0 : 1;
    } catch (const GmshReadError& e) {
        std::cerr << "mesh: " << e.what() << "\n";
        return 1;
    } catch (const InputError& e) {
        std::cerr << e.what() << "\n";
        return 1;
    } catch (const std::exception& e) {
        std::cerr << argv[1] << ": " << e.what() << "\n";
        return 1;
    }
}
