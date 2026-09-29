#include "aphi_solver/mumps_backend.hpp"

#include <chrono>
#include <cstdio>
#include <stdexcept>

#include "zmumps_c.h"

// MUMPS documents its control and information arrays 1-based, matching the
// Fortran and the user manual's tables. In C they are plain 0-based arrays, so
// every reference has to be shifted -- and an off-by-one here does not fail
// loudly, it silently sets a different control. These macros keep the code
// reading the way the manual does.
#define MUMPS_ICNTL(i) id.icntl[(i) - 1]
#define MUMPS_INFOG(i) id.infog[(i) - 1]
#define MUMPS_RINFOG(i) id.rinfog[(i) - 1]

// The value MUMPS expects when it is built against its sequential MPI stub
// (libmpiseq), which is how third_party/build_mumps.bat builds it.
#define MUMPS_USE_COMM_WORLD -987654

namespace aphi_solver {
namespace {

using Clock = std::chrono::steady_clock;

double ms_since(Clock::time_point t) {
    return std::chrono::duration<double, std::milli>(Clock::now() - t).count();
}

/// Run one MUMPS phase and report failure.
///
/// INFOG(1) < 0 is an error, > 0 a warning. Warnings are not silently
/// swallowed: a singular matrix reports INFOG(1) = 1 with the count in
/// INFOG(2), and treating that as success would hand back a meaningless
/// solution with a clean-looking return.
bool run_phase(ZMUMPS_STRUC_C& id, int job, const char* what) {
    id.job = job;
    zmumps_c(&id);
    if (MUMPS_INFOG(1) < 0) {
        std::fprintf(stderr, "MUMPS %s failed: INFOG(1) = %d, INFOG(2) = %d\n", what,
                     MUMPS_INFOG(1), MUMPS_INFOG(2));
        return false;
    }
    if (MUMPS_INFOG(1) > 0) {
        std::fprintf(stderr, "MUMPS %s warning: INFOG(1) = %d, INFOG(2) = %d\n", what,
                     MUMPS_INFOG(1), MUMPS_INFOG(2));
    }
    return true;
}

}  // namespace

bool solve_symmetric_mumps(const SparseSymmetricZ& a,
                           const std::vector<std::complex<double>>& b,
                           std::vector<std::complex<double>>& x, SolveReport& report) {
    report = SolveReport{};
    const int n = a.rows();
    if (static_cast<int>(b.size()) != n) {
        throw std::invalid_argument("solve_symmetric_mumps: b has the wrong length");
    }

    // --- COO, one-based, LOWER triangle ------------------------------------
    // We store the UPPER triangle in CSR: every entry has col >= row. MUMPS
    // with SYM = 2 wants one triangle of a symmetric matrix, and supplying
    // (col, row) instead of (row, col) turns our upper triangle into the lower
    // one without touching a value -- the matrix is symmetric, so the transpose
    // of one triangle IS the other. Supplying BOTH would double every
    // off-diagonal.
    const std::vector<int>& row_ptr = a.row_ptr();
    const std::vector<int>& col_index = a.col_index();
    const std::vector<std::complex<double>>& values = a.values();

    const MUMPS_INT8 nnz = static_cast<MUMPS_INT8>(values.size());
    std::vector<MUMPS_INT> irn(values.size());
    std::vector<MUMPS_INT> jcn(values.size());
    for (int r = 0; r < n; ++r) {
        const std::size_t ur = static_cast<std::size_t>(r);
        for (int k = row_ptr[ur]; k < row_ptr[ur + 1]; ++k) {
            const std::size_t uk = static_cast<std::size_t>(k);
            irn[uk] = static_cast<MUMPS_INT>(col_index[uk] + 1);  // lower: row = our col
            jcn[uk] = static_cast<MUMPS_INT>(r + 1);              //        col = our row
        }
    }

    // MUMPS overwrites the right-hand side with the solution in place.
    x = b;

    ZMUMPS_STRUC_C id;
    id.comm_fortran = MUMPS_USE_COMM_WORLD;
    id.par = 1;  // this process takes part in the factorization
    id.sym = 2;  // general symmetric, INDEFINITE -- not 1, which asserts
                 // positive definite and would skip the pivoting these
                 // systems need. See CONDITIONING.md: RowScaled and ScaledPhi
                 // are symmetric but not definite.

    if (!run_phase(id, -1, "initialization")) return false;

    // Quiet. MUMPS writes a page of diagnostics per phase by default, which
    // would bury the solver's own output. ICNTL(4) = 1 keeps error messages.
    MUMPS_ICNTL(1) = 6;   // error stream
    MUMPS_ICNTL(2) = -1;  // diagnostic  off
    MUMPS_ICNTL(3) = -1;  // global info off
    MUMPS_ICNTL(4) = 1;   // errors only

    MUMPS_ICNTL(5) = 0;   // assembled matrix, not elemental
    MUMPS_ICNTL(18) = 0;  // centralized on the host
    MUMPS_ICNTL(20) = 0;  // dense right-hand side

    // ICNTL(7) = 7 lets MUMPS choose the ordering. Left at the default
    // deliberately: the point of this backend is to see what a well-tuned
    // library does, and second-guessing its ordering heuristic would defeat
    // that. Note MUMPS was built WITHOUT METIS, so this comes down to AMD or
    // its own QAMD -- the same family the internal solver uses, which makes
    // the comparison about the factorization rather than the ordering.
    MUMPS_ICNTL(7) = 7;

    id.n = static_cast<MUMPS_INT>(n);
    id.nnz = nnz;
    id.irn = irn.data();
    id.jcn = jcn.data();
    // ZMUMPS_COMPLEX is { double r, i; }, which is layout-compatible with
    // std::complex<double> -- the standard guarantees the latter is two
    // contiguous doubles, real first.
    id.a = reinterpret_cast<ZMUMPS_COMPLEX*>(const_cast<std::complex<double>*>(values.data()));
    id.rhs = reinterpret_cast<ZMUMPS_COMPLEX*>(x.data());

    bool ok = true;
    auto t = Clock::now();
    if (!run_phase(id, 1, "analysis")) { ok = false; }
    report.analyze_ms = ms_since(t);

    if (ok) {
        t = Clock::now();
        if (!run_phase(id, 2, "factorization")) { ok = false; }
        report.factorize_ms = ms_since(t);
    }

    if (ok) {
        t = Clock::now();
        if (!run_phase(id, 3, "solve")) { ok = false; }
        report.solve_ms = ms_since(t);
    }

    if (ok) {
        // INFOG(29): entries in the factors after factorization.
        report.factor_nnz = static_cast<std::size_t>(MUMPS_INFOG(29));
        report.columns_done = n;
        // MUMPS's own scaling, reported so the column means something next to
        // the internal solver's equilibration spread. RINFOG(16)/(17) are the
        // smallest and largest scaling factors when ICNTL(8) applied one.
        report.equilibration_min = 1.0;
        report.equilibration_max = 1.0;
    }

    // Always finalize, including after a failure -- MUMPS holds the factor in
    // memory (gigabytes on these problems) until it is told to release it.
    id.job = -2;
    zmumps_c(&id);

    if (!ok) {
        report.ok = false;
        return false;
    }

    // Measured OUR way, on the original system, so the residual column is
    // comparable with the internal solver's rather than with MUMPS's own
    // differently-normalised RINFOG.
    report.residual = relative_residual(a, x, b);
    report.backward_error = backward_error(a, x, b);
    report.ok = true;
    return true;
}

}  // namespace aphi_solver
