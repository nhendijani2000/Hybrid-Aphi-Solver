#include "aphi_solver/symbolic.hpp"

#include <algorithm>
#include <stdexcept>
#include <string>

namespace aphi_solver {

namespace {

/// Neighbours of `k` that come before it, i.e. the entries of column `k` above
/// the diagonal, which is what both algorithms below consume.
///
/// Taken from the symmetric adjacency rather than from the pattern's rows, so
/// an upper-triangle pattern works as well as a full one: row `k` of an
/// upper-only pattern does not hold its own earlier columns, but the adjacency
/// does.
inline void earlier_neighbours(const SymmetricAdjacency& adj, int k, std::vector<int>& out) {
    out.clear();
    for (int t = adj.offset[static_cast<std::size_t>(k)];
         t < adj.offset[static_cast<std::size_t>(k) + 1]; ++t) {
        const int i = adj.neighbor[static_cast<std::size_t>(t)];
        if (i < k) out.push_back(i);
    }
}

/// Row `k` of `L`, as the set of tree nodes reachable from `A(k, 0:k-1)`.
///
/// For each earlier neighbour, walk up the elimination tree marking nodes until
/// one already marked for this row is met. Every node visited is a nonzero of
/// row `k`; the marking is what keeps the whole thing linear in the size of the
/// row rather than quadratic.
///
/// Returns the count and, when `out` is non-null, the column indices --
/// **unsorted**, since the walk produces tree order. The caller sorts if it
/// needs ascending order.
int row_reach(const SymmetricAdjacency& adj, const std::vector<int>& parent, int k,
              std::vector<int>& earlier, std::vector<int>& marked, std::vector<int>* out) {
    earlier_neighbours(adj, k, earlier);
    marked[static_cast<std::size_t>(k)] = k;
    int count = 1;  // the diagonal
    if (out != nullptr) out->push_back(k);

    for (int start : earlier) {
        int i = start;
        while (i >= 0 && marked[static_cast<std::size_t>(i)] != k) {
            marked[static_cast<std::size_t>(i)] = k;
            ++count;
            if (out != nullptr) out->push_back(i);
            i = parent[static_cast<std::size_t>(i)];
        }
    }
    return count;
}

}  // namespace

std::vector<int> elimination_tree(const SparsityPattern& permuted) {
    if (permuted.rows != permuted.cols) {
        throw std::invalid_argument("elimination_tree: the pattern must be square");
    }
    const int n = permuted.rows;
    const SymmetricAdjacency adj = SymmetricAdjacency::build(permuted);

    std::vector<int> parent(static_cast<std::size_t>(n), -1);
    std::vector<int> ancestor(static_cast<std::size_t>(n), -1);
    std::vector<int> earlier;
    earlier.reserve(64);

    // The standard path-compressing construction: for each k, walk from each
    // earlier neighbour towards the root, and the first node whose walk has not
    // yet reached anything gets k as its parent. `ancestor` is the compression,
    // and it is what keeps this near-linear rather than quadratic.
    for (int k = 0; k < n; ++k) {
        earlier_neighbours(adj, k, earlier);
        for (int start : earlier) {
            int i = start;
            while (i != -1 && i < k) {
                const int next = ancestor[static_cast<std::size_t>(i)];
                ancestor[static_cast<std::size_t>(i)] = k;
                if (next == -1) parent[static_cast<std::size_t>(i)] = k;
                i = next;
            }
        }
    }
    return parent;
}

std::vector<int> factor_row_counts(const SparsityPattern& permuted,
                                   const std::vector<int>& parent) {
    const int n = permuted.rows;
    if (static_cast<int>(parent.size()) != n) {
        throw std::invalid_argument(
            "factor_row_counts: the elimination tree has " + std::to_string(parent.size()) +
            " entries for a pattern of " + std::to_string(n) + " rows");
    }
    const SymmetricAdjacency adj = SymmetricAdjacency::build(permuted);

    std::vector<int> count(static_cast<std::size_t>(n), 0);
    std::vector<int> marked(static_cast<std::size_t>(n), -1);
    std::vector<int> earlier;
    earlier.reserve(64);
    for (int k = 0; k < n; ++k) {
        count[static_cast<std::size_t>(k)] = row_reach(adj, parent, k, earlier, marked, nullptr);
    }
    return count;
}

FactorSize predict_factor_size(const SparsityPattern& pattern, const Permutation& p) {
    const SparsityPattern permuted = permute_pattern(pattern, p);
    const std::vector<int> parent = elimination_tree(permuted);
    const std::vector<int> count = factor_row_counts(permuted, parent);

    FactorSize s;
    for (int c : count) {
        s.nnz += static_cast<std::size_t>(c);
        s.max_row = std::max(s.max_row, c);
    }
    for (int j = 0; j < static_cast<int>(parent.size()); ++j) {
        if (parent[static_cast<std::size_t>(j)] < 0) ++s.num_roots;
    }
    s.mean_row = permuted.rows > 0 ? static_cast<double>(s.nnz) / permuted.rows : 0.0;
    return s;
}

SolverAnalysis analyze(const SparsityPattern& pattern, Ordering ordering,
                       std::size_t max_factor_nnz) {
    SolverAnalysis a;
    a.ordering = ordering;
    a.permutation = compute_ordering(pattern, ordering);
    a.permuted = permute_pattern(pattern, a.permutation);
    a.parent = elimination_tree(a.permuted);
    a.row_count = factor_row_counts(a.permuted, a.parent);

    for (int c : a.row_count) a.predicted_nnz += static_cast<std::size_t>(c);
    for (int j = 0; j < static_cast<int>(a.parent.size()); ++j) {
        if (a.parent[static_cast<std::size_t>(j)] < 0) ++a.num_roots;
    }

    // Checked here, before the pattern is allocated: the counts are O(n) memory
    // whatever the fill, so an ordering that would be ruinous says so instead of
    // trying.
    if (max_factor_nnz != 0 && a.predicted_nnz > max_factor_nnz) {
        throw std::invalid_argument(
            "analyze: the " + std::string(ordering_keyword(ordering)) + " ordering predicts " +
            std::to_string(a.predicted_nnz) + " nonzeros in L (" +
            std::to_string(FactorSize{a.predicted_nnz, a.num_roots, 0, 0.0}.bytes(a.permuted.rows) /
                           1048576) +
            " MB as a complex factor), over the budget of " + std::to_string(max_factor_nnz) +
            ". Refusing before allocating it.");
    }

    // Build L's pattern. Two passes, as `build_sparsity` does: the row counts
    // above already give `row_ptr` exactly, so `col_index` is allocated once at
    // its final size rather than grown.
    const SymmetricAdjacency adj = SymmetricAdjacency::build(a.permuted);
    const int n = a.permuted.rows;
    a.factor.rows = n;
    a.factor.cols = n;
    a.factor.upper_only = false;  // L is LOWER triangular; this flag means "upper"
    a.factor.row_ptr.assign(static_cast<std::size_t>(n) + 1, 0);
    for (int k = 0; k < n; ++k) {
        a.factor.row_ptr[static_cast<std::size_t>(k) + 1] =
            a.factor.row_ptr[static_cast<std::size_t>(k)] + a.row_count[static_cast<std::size_t>(k)];
    }
    a.factor.col_index.assign(static_cast<std::size_t>(a.predicted_nnz), 0);

    std::vector<int> marked(static_cast<std::size_t>(n), -1);
    std::vector<int> earlier;
    earlier.reserve(64);
    std::vector<int> row;
    row.reserve(64);
    for (int k = 0; k < n; ++k) {
        row.clear();
        const int count = row_reach(adj, a.parent, k, earlier, marked, &row);
        if (count != a.row_count[static_cast<std::size_t>(k)]) {
            throw std::logic_error(
                "analyze: row " + std::to_string(k) + " of L counted " +
                std::to_string(a.row_count[static_cast<std::size_t>(k)]) +
                " entries but produced " + std::to_string(count) +
                ". The counting and building passes disagree.");
        }
        // The reach comes out in tree order; a numeric factorization
        // binary-searches these rows, so they are sorted here.
        std::sort(row.begin(), row.end());
        std::copy(row.begin(), row.end(),
                  a.factor.col_index.begin() + a.factor.row_ptr[static_cast<std::size_t>(k)]);
    }

    // Column counts, by transposing what was just built. An up-looking
    // factorization appends to columns as it walks rows, so it needs these to
    // lay out its storage before the first value exists.
    a.col_count.assign(static_cast<std::size_t>(n), 0);
    for (int j : a.factor.col_index) ++a.col_count[static_cast<std::size_t>(j)];
    return a;
}

namespace {

template <typename Matrix>
SparsityPattern pattern_from(const Matrix& a, bool upper_only) {
    SparsityPattern p;
    p.rows = a.rows();
    p.cols = a.cols();
    p.upper_only = upper_only;
    p.row_ptr = a.row_ptr();
    p.col_index = a.col_index();
    return p;
}

}  // namespace

SparsityPattern pattern_of(const SparseSymmetricZ& a) { return pattern_from(a, true); }
SparsityPattern pattern_of(const SparseMatrixZ& a) { return pattern_from(a, false); }

}  // namespace aphi_solver
