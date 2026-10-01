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
#include "aphi_solver/mumps_backend.hpp"
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
    if (argc < 2) {
        std::cerr << "A-Phi solver " << kVersion << " -- assemble, factorize, solve one file\n\n"
                  << "usage: solve_mesh <input-file> [natural|rcm|amd] [--compare-plain]\n\n"
                  << "The ordering defaults to amd. The conditioning comes from the file's\n"
                  << "[solver] section; only the symmetric ones can be solved yet, since the\n"
                  << "unsymmetric LU is step 8 of docs/SOLVER_PLAN.md.\n\n"
                  << "--compare-plain solves each frequency a SECOND time without\n"
                  << "equilibration, so the two residuals can be compared. That is a whole\n"
                  << "extra factorization and it doubles the runtime; the answer written to\n"
                  << "disk is the equilibrated one either way.\n\n"
                  << "try: solve_mesh examples/loop_sweep.aphi\n";
        return 2;
    }

    // Equilibration is on for the solve that produces the answer. The
    // unequilibrated one exists only to put `resid plain` beside `resid equil`,
    // and it costs a second factorization -- 97 % of a large run is
    // factorization, so that diagnostic doubles the wall time. Opt-in.
    Ordering ordering = Ordering::ApproximateMinimumDegree;
    bool compare_plain = false;
    for (int i = 2; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--compare-plain") {
            compare_plain = true;
        } else if (!ordering_from_keyword(argv[i], ordering)) {
            std::cerr << "unknown argument '" << arg
                      << "'; expected natural, rcm, amd or --compare-plain\n";
            return 2;
        }
    }

    const Clock::time_point run_started = Clock::now();
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

        // Gauge completeness, printed because it is the one property whose
        // failure is invisible downstream: a tree missing edges leaves part of
        // the gradient nullspace unpinned, so A can grow in those directions
        // while B = curl A stays correct and every residual looks fine.
        // A genuine spanning tree has exactly (nodes - components) edges.
        int tree = 0, dirichlet = 0, free_edges = 0;   // kept for RunInfo
        {
            std::vector<bool> on_dirichlet(static_cast<std::size_t>(mesh.num_nodes()), false);
            for (int e = 0; e < mesh.num_edges(); ++e) {
                switch (dofs.edge_state[static_cast<std::size_t>(e)]) {
                    case EdgeDof::Tree: ++tree; break;
                    case EdgeDof::Free: ++free_edges; break;
                    case EdgeDof::Dirichlet: {
                        ++dirichlet;
                        const std::pair<int, int>& ends = mesh.edges[static_cast<std::size_t>(e)];
                        on_dirichlet[static_cast<std::size_t>(ends.first)] = true;
                        on_dirichlet[static_cast<std::size_t>(ends.second)] = true;
                        break;
                    }
                }
            }
            int boundary_nodes = 0;
            for (bool b : on_dirichlet) {
                if (b) ++boundary_nodes;
            }

            // Count in the GROUP view, which is the one the interior tree
            // spans: each Dirichlet surface is one group (its own surface tree
            // edges hold it together, and those are classified Dirichlet here,
            // not Tree), every other node is a singleton. One connected mesh
            // with one Dirichlet surface therefore wants groups - 1 interior
            // tree edges. Comparing against nodes - 1 instead counts the
            // surface tree twice over and reports a deficit that is not there.
            const int groups = 1 + (mesh.num_nodes() - boundary_nodes);
            std::cout << "edges         " << mesh.num_edges() << " total = " << free_edges
                      << " free + " << tree << " tree + " << dirichlet << " Dirichlet\n";
            std::cout << "gauge         " << groups << " groups (" << boundary_nodes
                      << " boundary nodes as one) wants " << (groups - 1)
                      << " interior tree edges, has " << tree;
            std::cout << (tree == groups - 1 ? "   complete\n" : "   *** INCOMPLETE ***\n");
        }
        const SparsityPattern pattern =
            build_sparsity(dofs, bound, mesh, SparsityStorage::UpperTriangle);

        std::cout << "A-Phi solver " << kVersion << " -- solve_mesh\n\n";
        std::cout << "input         " << argv[1] << "\n";
        std::cout << "mesh          " << p.mesh_file << "   " << mesh.num_tets() << " tets\n";
        std::cout << "backend       " << solver_backend_keyword(p.backend)
                  << (p.backend == SolverBackend::Internal ? "   (this project's LDL^T)\n"
                                                           : "   (MUMPS)\n");
        std::cout << "conditioning  " << conditioning_keyword(p.conditioning)
                  << "   (symmetric, upper triangle stored)\n";
        std::cout << "unknowns      " << dofs.num_total << "   " << pattern.nnz()
                  << " stored nonzeros\n";
        for (const std::string& w : bound.warnings) std::cout << "warning       " << w << "\n";

        // Symbolic, once: it does not depend on frequency.
        //
        // SKIPPED ENTIRELY for MUMPS, which does its own analysis and never
        // sees ours. On case 02 this ordering costs 9 s -- more than MUMPS's
        // whole factorization -- so computing it only to print nnz(L) would be
        // 40 % of the run. MUMPS's real fill-in comes back in INFOG(29).
        //
        // For the internal solver it is computed once here and PASSED IN, so
        // solve_symmetric does not repeat it. It used to be done twice on
        // every run: once for this progress line, once inside the solver.
        const bool need_analysis = p.backend != SolverBackend::Mumps;
        SolverAnalysis analysis;
        if (need_analysis) {
            auto t0 = Clock::now();
            analysis = analyze(pattern, ordering);
            auto t1 = Clock::now();
            std::cout << "\nordering      " << ordering_keyword(ordering) << "   nnz(L) = "
                      << analysis.predicted_nnz << "   "
                      << analysis.predicted_nnz * 20.0 / 1048576.0 << " MB   [" << ms(t0, t1)
                      << " ms]\n\n";
        } else {
            std::cout << "\nordering      by MUMPS; ours skipped (it would be discarded)\n\n";
        }

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
        double write_ms = 0.0, fields_total_ms = 0.0, assemble_total_ms = 0.0;
        std::size_t last_factor_nnz = 0;  // MUMPS INFOG(29); ours comes from analysis
        double factor_total_ms = 0.0, solve_total_ms = 0.0;
        for (double f : frequencies) {
            const double omega = 2.0 * 3.14159265358979323846 * f;
            const Clock::time_point asm_t0 = Clock::now();
            refill(system, bound, mesh, dofs, pattern, omega, p.conditioning);
            const double assemble_ms = ms(asm_t0, Clock::now());

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
            // The second factorization, and the reason a run can take twice as
            // long as it needs to. Skipped unless asked for.
            const bool ok_plain =
                !compare_plain ||
                solve_symmetric(system.matrix, system.rhs, x_plain, rp, plain, &analysis);
            const bool ok_equil =
                p.backend == SolverBackend::Mumps
#ifdef APHI_WITH_MUMPS
                    ? solve_symmetric_mumps(system.matrix, system.rhs, x_equil, re)
#else
                    ? false  // unreachable: rejected at parse time
#endif
                    : solve_symmetric(system.matrix, system.rhs, x_equil, re, equil, &analysis);

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
            if (compare_plain) {
                std::cout << std::setw(14) << rp.residual;
            } else {
                std::cout << std::setw(14) << "--";
            }
            std::cout << std::setw(14) << re.residual
                      << std::setw(14) << re.backward_error
                      << std::setw(12) << re.equilibration_max / re.equilibration_min
                      << std::setw(13) << re.smallest_pivot << std::setw(13) << re.largest_pivot
                      << std::setw(12) << re.largest_multiplier << std::setw(11) << re.factorize_ms
                      << "\n";
            if (!(re.backward_error < 1e-12)) ++failures;

            // The potential, written for every frequency. One file when there is
            // only one solve, numbered otherwise, so a sweep does not silently
            // overwrite itself.
            // What produced the files about to be written. Built once the solve
            // has reported, so it can carry the convergence numbers too.
            RunInfo run;
            run.input_file = argv[1];
            run.mesh_file = p.mesh_file;
            run.ordering = ordering_keyword(ordering);
            run.num_tets = mesh.num_tets();
            run.unknowns = dofs.num_total;
            run.free_edges = free_edges;
            run.tree_edges = tree;
            run.dirichlet_edges = dirichlet;
            run.stored_nonzeros = pattern.nnz();
            run.factor_nnz = re.factor_nnz;
            run.backward_error = re.backward_error;
            run.residual = re.residual;
            run.assemble_ms = assemble_ms;
            run.analyze_ms = re.analyze_ms;
            run.factorize_ms = re.factorize_ms;
            run.solve_ms = re.solve_ms;

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
            const WriteStats ws = write_potential(out_path, nodal, sol, &run);
            write_ms += ws.milliseconds;
            const auto fields_started = std::chrono::steady_clock::now();
            const FieldOutput fields = compute_fields(mesh, bound, dofs, sol);
            const double fields_ms =
                std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() -
                                                          fields_started)
                    .count();
            last_factor_nnz = re.factor_nnz;
            fields_total_ms += fields_ms;
            assemble_total_ms += assemble_ms;
            // Each frequency is factorized TWICE -- plain and equilibrated -- so the
            // two residuals can sit side by side. Both cost real time, so both count;
            // charging only the equilibrated one left half the runtime in "other".
            factor_total_ms += (compare_plain ? rp.factorize_ms : 0.0) + re.factorize_ms;
            solve_total_ms += (compare_plain ? rp.solve_ms : 0.0) + re.solve_ms;
            const std::string stem = out_path.substr(0, out_path.size() - 4);
            std::string vtk_path = stem + ".vtk";
            const WriteStats vs =
                write_vtk(vtk_path, mesh, nodal, sol, &fields, FieldSet::Potential, &run);
            write_ms += vs.milliseconds;
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
                const WriteStats ts = write_solution(base + ".out", nodal, fields, sol, fset, &run);
                const WriteStats vv =
                    write_vtk(base + ".vtk", mesh, nodal, sol, &fields, fset, &run);
                std::cout << "                wrote " << base << ".out / .vtk   "
                          << ts.bytes / 1024 << " + " << vv.bytes / 1024 << " KB, "
                          << (ts.milliseconds + vv.milliseconds) << " ms\n";
                write_ms += ts.milliseconds + vv.milliseconds;
            }

            // The [postprocess] requests, for tools/postprocess.py. Written
            // beside the fields it refers to and with the same suffix, so a
            // sweep's solves do not read each other's pictures. Nothing is
            // written when the file asked for none.
            const std::string manifest =
                dir + "postprocess" +
                (suffix == "potential" ? "" : suffix.substr(std::strlen("potential"))) + ".json";
            const std::string msuffix =
                suffix == "potential" ? std::string() : suffix.substr(std::strlen("potential"));
            const WriteStats ms =
                write_postprocess_manifest(manifest, p, bound, f, msuffix, &run);
            if (ms.bytes > 0) {
                std::cout << "                wrote " << manifest << "   " << ms.nodes
                          << " request(s), " << ms.bytes << " bytes\n";
                write_ms += ms.milliseconds;
            }
        }

        // Where the time went. Printed unconditionally because a run that is
        // slow for the wrong reason -- paging, or an ordering that blew up --
        // looks identical to a healthy one in every other line of output.
        //
        // Built once into a string, then sent BOTH to the console and to
        // output/run_summary.txt. Console-only was the first version, and it
        // meant the one number people actually go looking for afterwards --
        // where did the time go -- survived only as long as the scrollback.
        {
            const double total = ms(run_started, Clock::now());
            const double accounted = assemble_total_ms + factor_total_ms + solve_total_ms +
                                     fields_total_ms + write_ms;
            std::ostringstream s;
            s << "\n" << std::left << std::setw(22) << "timing" << std::right
              << std::setw(12) << "seconds" << std::setw(10) << "percent" << "\n";
            const auto row = [&](const char* name, double v) {
                s << std::left << std::setw(22) << name << std::right << std::fixed
                  << std::setprecision(2) << std::setw(12) << v / 1000.0
                  << std::setw(9) << (total > 0.0 ? 100.0 * v / total : 0.0) << " %\n";
            };
            row("  matrix assembly", assemble_total_ms);
            row(compare_plain ? "  factorization x2" : "  factorization", factor_total_ms);
            row("  triangular solve", solve_total_ms);
            row("  field computation", fields_total_ms);
            row("  writing fields", write_ms);
            row("  mesh, dofs, other", total - accounted);
            s << std::left << std::setw(22) << "  TOTAL" << std::right << std::fixed
              << std::setprecision(2) << std::setw(12) << total / 1000.0
              << std::setw(9) << 100.0 << " %\n" << std::defaultfloat;
            if (compare_plain) {
                s << "\n  factorization is counted TWICE per frequency: --compare-plain ran\n"
                     "  the solve plain as well as equilibrated so the two residuals could\n"
                     "  be compared. Half of that line is diagnostic.\n";
            } else {
                s << "\n  the solve ran ONCE, equilibrated. Pass --compare-plain to also solve\n"
                     "  it unequilibrated and compare residuals -- that is a second\n"
                     "  factorization and roughly doubles the runtime.\n";
            }
            std::cout << s.str();

            if (!p.output_dir.empty()) {
                const std::string sum_path = p.output_dir + "/run_summary.txt";
                if (std::FILE* sf = std::fopen(sum_path.c_str(), "wb")) {
                    std::ostringstream h;
                    h << "A-Phi solver " << kVersion << " -- run summary\n\n"
                      << "input         " << argv[1] << "\n"
                      << "mesh          " << p.mesh_file << "   " << mesh.num_tets() << " tets\n"
                      << "conditioning  " << conditioning_keyword(p.conditioning) << "\n"
                      << "unknowns      " << dofs.num_total << "   (A: " << free_edges
                      << " free edges, " << tree << " gauged to zero, " << dirichlet
                      << " Dirichlet)\n"
                      << "matrix        " << pattern.nnz() << " stored nonzeros, ordering "
                      << ordering_keyword(ordering) << "\n"
                      << "fill-in       "
                      << (need_analysis ? "nnz(L) = " : "by MUMPS, INFOG(29) = ")
                      << (need_analysis ? analysis.predicted_nnz : last_factor_nnz) << "\n";
                    const std::string text = h.str() + s.str();
                    std::fwrite(text.data(), 1, text.size(), sf);
                    std::fclose(sf);
                    std::cout << "\n                wrote " << sum_path << "\n";
                } else {
                    std::cout << "\n                could not write " << sum_path << "\n";
                }
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
