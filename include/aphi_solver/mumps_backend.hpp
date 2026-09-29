#pragma once

/// MUMPS as an optional factorization backend.
///
/// Declared unconditionally so `solve_mesh` can dispatch without `#ifdef` in
/// the caller; DEFINED only when the build was configured with
/// `-DAPHI_WITH_MUMPS=ON`. Without it the translation unit is not compiled and
/// nothing can reach here anyway, because `[solver] backend = mumps` is
/// rejected at parse time (`solver_backend_available`).
///
/// See `docs/MUMPS_SETUP.md` for what has to be installed, and
/// `docs/SOLVER_PLAN.md` Sec. 12 for why this is an option and not a
/// dependency.

#include <complex>
#include <vector>

#include "aphi_solver/factorization.hpp"
#include "aphi_solver/sparse_symmetric.hpp"

namespace aphi_solver {

/// Solve `a x = b` with MUMPS, filling `report` the same way
/// `solve_symmetric` does.
///
/// The point of this backend is a like-for-like timing comparison, so the
/// numbers it reports are measured the SAME WAY as the internal solver's:
/// `residual` and `backward_error` come from `relative_residual` and
/// `backward_error` applied to the original `a`, `x`, `b` -- not from MUMPS's
/// own `RINFOG`, which uses a different normalisation and would make the two
/// columns incomparable.
///
/// `analyze_ms`, `factorize_ms` and `solve_ms` are taken from MUMPS's own
/// phases (JOB 1, 2, 3 run separately rather than JOB 6) so they line up with
/// the internal solver's breakdown. `factor_nnz` is MUMPS's `INFOG(29)`.
///
/// Equilibration is left to MUMPS (`ICNTL(8)`, its default automatic choice)
/// rather than applying ours first: scaling twice would be measuring neither.
///
/// Returns false with `report.ok = false` if MUMPS reports an error; the
/// message is printed to stderr with its `INFOG(1)` and `INFOG(2)` codes,
/// which is what the MUMPS manual's error table is indexed by.
bool solve_symmetric_mumps(const SparseSymmetricZ& a,
                           const std::vector<std::complex<double>>& b,
                           std::vector<std::complex<double>>& x, SolveReport& report);

}  // namespace aphi_solver
