#pragma once

#include <complex>
#include <vector>

#include "aphi_solver/equilibration.hpp"
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

    /// The largest |L(i,j)| produced. **This, not the pivot range, is what says
    /// whether an unpivoted factorization was stable.** A well-scaled pivot can
    /// still divide into a large numerator, and the resulting multiplier amplifies
    /// every later update -- which is exactly what pivoting exists to bound.
    /// Measured because the pivot range alone proved misleading on a real system
    /// (docs/SOLVER_PLAN.md Sec. 18).
    double largest_multiplier = 0.0;
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

/// The standard relative BACKWARD error,
///
///     ||A x - b|| / ( max|A_ij| * ||x|| + ||b|| )
///
/// and the number to judge a solve by. `relative_residual` above divides by
/// `||b||` alone, which is fine when the right-hand side is a comparable size to
/// the rest of the problem and badly misleading when it is not.
///
/// It is not here for tidiness. On the real loop system at 100 MHz under
/// `row_scaled`, `||b||` is 1.6e-9 with exactly ONE nonzero entry -- the port
/// row, carrying `r * I = I / (j*omega)` -- while the matrix reaches 4.4e12. The
/// ratio is 2.8e21, so dividing by `||b||` turns pure rounding noise into a
/// reported 6.7 %, while the backward error is 4.2e-25. One of those is a
/// statement about the solve and the other is a statement about the
/// normalisation. See `docs/SOLVER_PLAN.md` Sec. 18.
///
/// Note the dependence this exposes: `||b||` scales with the conditioning choice
/// (`Natural` would give `I = 1` in that row, `RowScaled` gives `I/(j*omega)`),
/// so a `||b||`-relative residual is not even comparable between formulations of
/// the same problem. The backward error is.
double backward_error(const SparseSymmetricZ& a, const std::vector<std::complex<double>>& x,
                     const std::vector<std::complex<double>>& b);


/// Options for the whole solve. `docs/SOLVER_PLAN.md` §7.
struct SolveOptions {
    Ordering ordering = Ordering::ApproximateMinimumDegree;

    /// Ruiz iterations for `diag(d) A diag(d)`. **0 disables scaling**, which is
    /// what the negative control for it uses.
    int equilibration_iterations = 10;

    /// A pivot at or below this, relative to the largest |D| seen so far, stops
    /// the factorization. Steps 7 onwards perturb instead.
    double pivot_floor = 0.0;
};

/// What one solve did, reported rather than inferred. `docs/SOLVER_PLAN.md` §5.
struct SolveReport {
    bool ok = false;

    /// `||A x - b|| / ( max|A| ||x|| + ||b|| )` on the ORIGINAL system: the
    /// number to judge the solve by. See `backward_error`.
    double backward_error = 0.0;

    /// `||A x - b|| / ||b||`, kept because it is what a user expects to see, and
    /// misleading whenever `||b||` is small compared with `||A|| ||x||`.
    double residual = 0.0;

    /// The residual the same solve would have had without equilibration is not
    /// computed here -- run it twice with `equilibration_iterations = 0` to get
    /// that, which is what `tools/solve_mesh.cpp` does.
    double equilibration_min = 1.0;  ///< smallest d
    double equilibration_max = 1.0;  ///< largest d

    double smallest_pivot = 0.0;
    double largest_pivot = 0.0;
    double largest_multiplier = 0.0;  ///< max |L(i,j)|; see FactorStats, and Sec. 18
    std::size_t factor_nnz = 0;
    int columns_done = 0;

    double analyze_ms = 0.0;
    double factorize_ms = 0.0;
    double solve_ms = 0.0;
};

/// Equilibrate, order, factorize, solve, un-scale, and measure.
///
/// The order is **scale, then order, then factorize**: ordering is structural so
/// scaling cannot affect it, and scaling improves pivot quality so it has to come
/// first.
///
/// The substitution and its undoing, which is the step that is easy to get wrong:
///
///     A~ = D A D,   solve  A~ y = D b,   then  x = D y
///
/// Both times a **multiplication** by `d` -- not `x = y`, and not `x = D⁻¹ y`.
/// Getting it wrong yields a smooth, plausible, completely wrong field, so one
/// test compares against an unscaled solve of the same system and one control
/// drops the recovery.
///
/// `report.residual` is always measured against the **original** `a` and `b`,
/// never the equilibrated pair, which is what makes it a true backward error.
///
/// Returns false, with `report` filled in as far as it got, if the factorization
/// hit a pivot at or below `options.pivot_floor`.
bool solve_symmetric(const SparseSymmetricZ& a, const std::vector<std::complex<double>>& b,
                     std::vector<std::complex<double>>& x, SolveReport& report,
                     const SolveOptions& options = {},
                     /// Already-analysed structure for this exact pattern, to avoid
                     /// repeating the ordering. `analyze` reads only the pattern, and
                     /// scaling does not change it, so reuse is exact -- not an
                     /// approximation. nullptr analyses internally as before.
                     const SolverAnalysis* precomputed = nullptr);

/// `L D Lᵀ` expanded back to a full matrix, for checking. Only for tests and
/// diagnostics -- it is quadratic in a column's length and defeats the whole
/// point of a sparse factor.
SparseMatrixZ reconstruct(const SymmetricFactor& f);

}  // namespace aphi_solver
