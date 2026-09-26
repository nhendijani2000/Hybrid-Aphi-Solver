#pragma once

#include <cstddef>
#include <vector>

#include "aphi_solver/ordering.hpp"
#include "aphi_solver/sparse_symmetric.hpp"
#include "aphi_solver/sparsity.hpp"

namespace aphi_solver {

/// The symbolic pass: everything about a factorization that does not depend on
/// a single value. `docs/SOLVER_PLAN.md` §2 and step 2 of §9.
///
/// This is the same split that paid for itself in assembly, for the same
/// reason. The matrix's structure does not depend on frequency and its values
/// do -- measured on the cylinder, 1557013 of 1557522 entries differ between
/// 1 kHz and 1 MHz while `row_ptr` and `col_index` are identical. So ordering,
/// elimination tree and fill pattern are computed **once per mesh** and reused
/// at every frequency.
///
/// Algorithms are the standard ones: the elimination tree by path compression,
/// and row patterns of `L` by the reachable-set search on that tree (Davis,
/// *Direct Methods for Sparse Linear Systems*, SIAM 2006; column counts after
/// Gilbert, Ng & Peyton). Implemented from the published descriptions.

/// The elimination tree of a symmetric pattern, already permuted.
///
/// `parent[j]` is the smallest `i > j` with `L(i,j) != 0`, or -1 when column
/// `j` has no off-diagonal entry -- which makes `j` a root. There is one root
/// per connected component, and an isolated unknown (a voltage port's
/// constraint row is one) is a root of its own.
///
/// `parent[j] > j` always, so the tree is a forest oriented towards higher
/// indices; that is what makes the reachable-set search below terminate.
std::vector<int> elimination_tree(const SparsityPattern& permuted);

/// Nonzeros in each row of `L`, including the diagonal, from the tree.
///
/// Row `k` of `L` is the set of nodes reachable from `A(k, 0:k-1)` by walking
/// up the elimination tree, plus `k` itself. Computing the counts without
/// storing the pattern is the point: on the cylinder with the natural ordering
/// the pattern is hundreds of millions of entries, and the counts are what say
/// so **before** anything tries to allocate it.
std::vector<int> factor_row_counts(const SparsityPattern& permuted,
                                   const std::vector<int>& parent);

/// What a factorization of `pattern` under `p` would cost, in structure only.
/// O(n) memory whatever the fill is, so it is safe to call on an ordering that
/// would be ruinous -- which is exactly when the answer matters.
struct FactorSize {
    std::size_t nnz = 0;   ///< nonzeros in L, including the diagonal
    int num_roots = 0;     ///< connected components of the pattern
    int max_row = 0;       ///< longest row of L
    double mean_row = 0.0;

    /// Bytes a complex factor would need: one complex value and one int column
    /// index per nonzero, plus the row pointers.
    std::size_t bytes(int rows) const {
        return nnz * (sizeof(double) * 2 + sizeof(int)) +
               (static_cast<std::size_t>(rows) + 1) * sizeof(int);
    }
};

FactorSize predict_factor_size(const SparsityPattern& pattern, const Permutation& p);

/// Everything the numeric factorization needs, computed once per mesh and
/// ordering.
struct SolverAnalysis {
    Ordering ordering = Ordering::Natural;
    Permutation permutation;

    /// The pattern actually factorized: `P A Pᵀ`.
    SparsityPattern permuted;

    std::vector<int> parent;     ///< elimination tree of `permuted`
    std::vector<int> row_count;  ///< nonzeros per row of L, diagonal included

    /// Nonzeros per COLUMN of L, diagonal included. Different from
    /// `row_count` -- row k of L holds columns j <= k, column j holds rows
    /// i >= j -- and it is the one an up-looking factorization needs, because
    /// that is what lays out its column storage before any value is known.
    /// Both sum to `predicted_nnz`, which is a check on each.
    std::vector<int> col_count;

    /// `L`'s pattern, **lower** triangular, in permuted indices. Each row
    /// ascending, so a numeric factorization can binary-search it.
    SparsityPattern factor;

    std::size_t predicted_nnz = 0;  ///< sum of `row_count`; checked against the built pattern
    int num_roots = 0;
};

/// Orders, permutes, builds the elimination tree, counts, and builds `L`'s
/// pattern.
///
/// `max_factor_nnz` is a budget checked **after** the counts and **before**
/// the pattern is allocated: exceeding it throws, naming the number, rather
/// than attempting an allocation of hundreds of megabytes. 0 means no limit.
/// This is not hypothetical -- on `cylinder_box.msh` the natural ordering
/// predicts a factor far larger than the machine should be asked for, which is
/// the whole argument for reordering and is reported as a number in
/// `SOLVER_PLAN.md` §14.
SolverAnalysis analyze(const SparsityPattern& pattern, Ordering ordering,
                       std::size_t max_factor_nnz = 0);

/// A matrix's own pattern, so a caller does not have to reconstruct one.
SparsityPattern pattern_of(const SparseSymmetricZ& a);
SparsityPattern pattern_of(const SparseMatrixZ& a);

/// Convenience: analyse a matrix directly. Equivalent to analysing its
/// pattern, and the form every caller actually wants.
inline SolverAnalysis analyze(const SparseSymmetricZ& a, Ordering ordering,
                             std::size_t max_factor_nnz = 0) {
    return analyze(pattern_of(a), ordering, max_factor_nnz);
}
inline SolverAnalysis analyze(const SparseMatrixZ& a, Ordering ordering,
                             std::size_t max_factor_nnz = 0) {
    return analyze(pattern_of(a), ordering, max_factor_nnz);
}
inline Permutation compute_ordering(const SparseSymmetricZ& a, Ordering ordering) {
    return compute_ordering(pattern_of(a), ordering);
}

}  // namespace aphi_solver
