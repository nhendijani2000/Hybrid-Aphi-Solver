#pragma once

#include <vector>

#include "aphi_solver/complex_matrix.hpp"

namespace aphi_solver {

/// Numerical (not physics-motivated) preconditioning step: diagonal equilibration
/// of a dense complex matrix before factorization. This is a distinct job from the
/// transforms in conditioning.hpp -- those fix the *structural*/frequency-domain
/// behavior of the A-Phi system (row scaling by j*omega, the Phi = j*omega*Phi'
/// substitution); this fixes the *numerical scale* of whatever system you hand it,
/// physics-transformed or not. Apply both, in either order, as needed -- see
/// docs/LINEAR_SOLVER.md.
///
/// The scaling here is deliberately symmetric: a single diagonal `d` applied as
/// diag(d) * A * diag(d), rather than independent row and column scalings. A is
/// complex symmetric (A = A^T, not Hermitian -- see conditioning.hpp), and
/// symmetric scaling is what preserves that structure (diag(d) * A * diag(d) is
/// still symmetric because diag(d) equals its own transpose; independent row/
/// column scaling would generally destroy the symmetry). Keeping A symmetric
/// matters downstream: it's what makes symmetric-indefinite factorization and
/// COCG/COCR-type iterative solvers applicable in the first place (Phase 05).

/// These operate on the sparse system type as of Phase 03.5
/// (`docs/ROADMAP.md`). Equilibration runs on the assembled matrix
/// immediately before factorization, so keeping it on the dense
/// `ComplexMatrix` would have forced a dense round-trip in the middle of the
/// production solve path -- the one place that cannot afford it. Both
/// operations are naturally sparse anyway: `diag(d) A diag(d)` scales each
/// stored value by `d[row] * d[col]` and never changes the sparsity pattern.

/// Computes a real, positive diagonal scaling vector `d` (length A.rows()) via
/// symmetric Ruiz-style iteration: repeatedly rescale so each row/column's
/// largest-magnitude entry approaches 1. A zero row is left unscaled (d = 1 there)
/// rather than dividing by zero. `iterations` of ~5-10 is normally enough for the
/// scaling to converge; more rarely helps.
std::vector<double> compute_symmetric_equilibration(const SparseMatrixZ& A, int iterations = 10);

/// Returns diag(d) * A * diag(d).
SparseMatrixZ apply_symmetric_equilibration(const SparseMatrixZ& A, const std::vector<double>& d);

/// Returns diag(d) * b (elementwise scaling of a right-hand-side vector).
std::vector<Complex> scale_rhs(const std::vector<Complex>& b, const std::vector<double>& d);

/// Recovers the physical solution x = diag(d) * y after solving the equilibrated
/// system (diag(d) A diag(d)) y = diag(d) b.
std::vector<Complex> recover_equilibrated_solution(const std::vector<Complex>& y, const std::vector<double>& d);

}  // namespace aphi_solver
