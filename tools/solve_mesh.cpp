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

        std::cout << std::right << "  " << std::setw(12) << "frequency" << std::setw(12)
                  << "factor ms" << std::setw(12) << "solve ms" << std::setw(13) << "residual"
                  << std::setw(13) << "min |D|" << std::setw(13) << "max |D|" << "\n";

        int failures = 0;
        for (double f : frequencies) {
            const double omega = 2.0 * 3.14159265358979323846 * f;
            refill(system, bound, mesh, dofs, pattern, omega, p.conditioning);

            SymmetricFactor factor;
            FactorStats stats;
            const bool ok = factorize_ldlt(system.matrix, analysis, factor, stats);

            std::ostringstream label;
            label << std::setprecision(4) << f << " Hz";
            std::cout << "  " << std::setw(12) << label.str() << std::setw(12)
                      << stats.milliseconds;
            if (!ok) {
                std::cout << "   FAILED at column " << stats.columns_done << " of "
                          << dofs.num_total << ", pivot " << stats.smallest_pivot << "\n";
                ++failures;
                continue;
            }

            auto s0 = Clock::now();
            const std::vector<Complex> x = solve(factor, analysis.permutation, system.rhs);
            auto s1 = Clock::now();
            const double residual = relative_residual(system.matrix, x, system.rhs);

            std::cout << std::setw(12) << ms(s0, s1) << std::setw(13) << residual << std::setw(13)
                      << stats.smallest_pivot << std::setw(13) << stats.largest_pivot << "\n";
            if (!(residual < 1e-6)) ++failures;
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
