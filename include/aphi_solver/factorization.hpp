#pragma once

#include <complex>
#include <vector>

#include "aphi_solver/sparse_matrix.hpp"
#include "aphi_solver/sparse_symmetric.hpp"
#include "aphi_solver/symbolic.hpp"

namespace aphi_solver {

/// Numeric factorization. Step 4 of `docs/SOLVER_PLAN.md` §9.
///
/// `P A Pᵀ = L D Lᵀ` with `L` **unit** lower triangular and `D` diagonal. Both
/// complex, neither positive: the matrix is complex symmetric (`A = Aᵀ`, not
/// `A = Aᴴ`) and indefinite, so Cholesky does not apply and nothing here may
/// conjugate. A Hermitian factorization applied to one of these matrices is
/// silently wrong rather than an error, which is why one test builds a matrix
/// whose reconstruction only comes out right if no conjugation happened.
///
/// Up-looking: row `k` of `L` is computed from rows already done, by walking the
/// elimination tree (Davis, *Direct Methods for Sparse Linear Systems*, SIAM
/// 2006, §4.8). Implemented from the published description. Scalar, no blocking:
/// Stage 1 is allowed to be slow.
///
/// **No pivoting at this step.** `LDLᵀ` without pivoting can meet a zero pivot
/// on a perfectly nonsingular complex symmetric matrix -- `[[0,1],[1,0]]` is the
/// smallest example -- so `factorize` reports that rather than dividing by it.
/// Step 7 adds static pivoting with perturbation and iterative refinement.

/// `L` and `D`. `L`'s unit diagonal is not stored, so `col_ptr` describes the
/// **strictly** lower triangle and holds `predicted_nnz - rows` entries.
struct SymmetricFactor {
    int rows = 0;
    std::vector<int> col_ptr;    ///< rows + 1
    std::vector<int> row_index;  ///< strictly-lower row indices, ascending per column
    std::vector<std::complex<double>> value;
    std::vector<std::complex<double>> diagonal;  ///< D

    std::size_t nnz() const { return value.size() + diagonal.size(); }
};

/// What the factorization did, reported rather than inferred.
struct FactorStats {
    bool ok = false;               ///< false when a pivot was too small to use
    int columns_done = 0;          ///< where it stopped, if it did
    double smallest_pivot = 0.0;   ///< min |D[j]| over the columns completed
    double largest_pivot = 0.0;
    std::size_t nnz = 0;           ///< L's stored entries plus D
    double milliseconds = 0.0;
};

/// The permuted matrix's **lower** triangle in CSR -- which is the same data as
/// its upper triangle by columns, and is exactly what the up-looking algorithm
/// consumes row by row.
///
/// `a` holds A's upper triangle in the original ordering. Every stored entry is
/// mapped to its permuted position and normalised below the diagonal, so this is
/// where a `perm`/`iperm` confusion would show up -- and one test checks the
/// result against `a.entry()` queried directly.
SparseMatrixZ permuted_lower(const SparseSymmetricZ& a, const Permutation& p);

/// Factorizes `P A Pᵀ` into `L D Lᵀ` using the structure in `analysis`.
///
/// `analysis` must have been produced from `a`'s pattern: the column counts lay
/// out `L`'s storage and the elimination tree drives the row walks, so an
/// analysis of a different matrix would write into the wrong places.
///
/// Returns `stats.ok == false`, having stopped, if `|D[j]|` falls at or below
/// `pivot_floor` -- with no pivoting there is nothing useful to do with such a
/// column, and continuing would produce a factor whose entries are meaningless.
/// The caller sees where it stopped and how small the pivot was.
bool factorize_ldlt(const SparseSymmetricZ& a, const SolverAnalysis& analysis,
                    SymmetricFactor& out, FactorStats& stats, double pivot_floor = 0.0);

/// Solves `L D Lᵀ y = y` in place, in the FACTORIZED ordering.
///
/// Three sweeps: forward through `L` (unit lower, so no division), then by
/// `D`, then backward through `Lᵀ`. The last uses the transpose and **not** the
/// conjugate transpose -- these matrices are complex symmetric, and a stray
/// conjugation here gives a wrong answer with a small-looking residual, since
/// the residual would be computed against the same wrong operator.
///
/// Throws std::invalid_argument on a length mismatch, and std::logic_error if a
/// diagonal entry is zero -- which cannot happen for a factor `factorize_ldlt`
/// returned true for, so it means the factor came from somewhere else.
void solve_in_place(const SymmetricFactor& f, std::vector<std::complex<double>>& y);

/// Solves `A x = b` for the original, unpermuted `A`.
///
/// The factorization is of `P A Pᵀ`, so with `b~[i] = b[perm[i]]` and
/// `x[perm[i]] = x~[i]`:
///
///     P A Pᵀ (P x) = P b
///
/// Both directions use **`perm`**, one on the way in and one on the way out.
/// Using `iperm` for either, or permuting the right-hand side and forgetting the
/// solution, gives a smooth plausible field that solves nothing -- which is why
/// one test solves the same system under two different orderings and requires
/// the same answer.
std::vector<std::complex<double>> solve(const SymmetricFactor& f, const Permutation& p,
                                       const std::vector<std::complex<double>>& b);

/// `||A x - b|| / ||b||` in the 2-norm, on the ORIGINAL system: unpermuted, and
/// unscaled by anything the solver did internally. That is what makes it a true
/// backward error rather than a statement about the factorization's own
/// arithmetic. Returns `||A x||` when `b` is zero.
double relative_residual(const SparseSymmetricZ& a, const std::vector<std::complex<double>>& x,
                         const std::vector<std::complex<double>>& b);

/// `L D Lᵀ` expanded back to a full matrix, for checking. Only for tests and
/// diagnostics -- it is quadratic in a column's length and defeats the whole
/// point of a sparse factor.
SparseMatrixZ reconstruct(const SymmetricFactor& f);

}  // namespace aphi_solver
