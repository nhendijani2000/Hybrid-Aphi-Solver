#include "aphi_solver/ordering.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

namespace aphi_solver {

const char* ordering_keyword(Ordering ordering) {
    switch (ordering) {
        case Ordering::Natural: return "natural";
        case Ordering::ReverseCuthillMcKee: return "rcm";
        case Ordering::ApproximateMinimumDegree: return "amd";
    }
    return "natural";
}

bool ordering_from_keyword(const std::string& word, Ordering& out) {
    if (word == "natural") {
        out = Ordering::Natural;
    } else if (word == "rcm") {
        out = Ordering::ReverseCuthillMcKee;
    } else if (word == "amd") {
        out = Ordering::ApproximateMinimumDegree;
    } else {
        return false;
    }
    return true;
}

bool Permutation::is_valid() const {
    if (perm.size() != iperm.size()) return false;
    const int n = static_cast<int>(perm.size());
    std::vector<char> hit(static_cast<std::size_t>(n), 0);
    for (int i = 0; i < n; ++i) {
        const int old = perm[static_cast<std::size_t>(i)];
        if (old < 0 || old >= n) return false;
        if (hit[static_cast<std::size_t>(old)]) return false;  // not injective
        hit[static_cast<std::size_t>(old)] = 1;
        if (iperm[static_cast<std::size_t>(old)] != i) return false;
    }
    return true;
}

Permutation Permutation::identity(int n) {
    Permutation p;
    p.perm.resize(static_cast<std::size_t>(std::max(n, 0)));
    std::iota(p.perm.begin(), p.perm.end(), 0);
    p.iperm = p.perm;
    return p;
}

namespace {

/// Undirected adjacency, no self-loops, each neighbour list ascending.
///
/// Built from whatever the pattern stores by treating `(i, j)` as an edge both
/// ways, so a full pattern, an upper triangle, and an unsymmetric pattern all
/// give the adjacency of the symmetric structure. The diagonal is dropped: a
/// self-loop is not a connection to anywhere and would distort every degree by
/// one.
struct Adjacency {
    std::vector<int> offset;    ///< rows + 1
    std::vector<int> neighbor;  ///< ascending within each row

    int degree(int i) const {
        return offset[static_cast<std::size_t>(i) + 1] - offset[static_cast<std::size_t>(i)];
    }
};

Adjacency build_adjacency(const SparsityPattern& pattern) {
    const int n = pattern.rows;
    Adjacency adj;
    adj.offset.assign(static_cast<std::size_t>(n) + 1, 0);

    // Count both directions.
    for (int r = 0; r < n; ++r) {
        for (int k = pattern.row_ptr[static_cast<std::size_t>(r)];
             k < pattern.row_ptr[static_cast<std::size_t>(r) + 1]; ++k) {
            const int c = pattern.col_index[static_cast<std::size_t>(k)];
            if (c == r) continue;
            ++adj.offset[static_cast<std::size_t>(r) + 1];
            ++adj.offset[static_cast<std::size_t>(c) + 1];
        }
    }
    for (int r = 0; r < n; ++r) {
        adj.offset[static_cast<std::size_t>(r) + 1] += adj.offset[static_cast<std::size_t>(r)];
    }

    std::vector<int> cursor(adj.offset.begin(), adj.offset.end() - 1);
    adj.neighbor.assign(static_cast<std::size_t>(adj.offset.back()), 0);
    for (int r = 0; r < n; ++r) {
        for (int k = pattern.row_ptr[static_cast<std::size_t>(r)];
             k < pattern.row_ptr[static_cast<std::size_t>(r) + 1]; ++k) {
            const int c = pattern.col_index[static_cast<std::size_t>(k)];
            if (c == r) continue;
            adj.neighbor[static_cast<std::size_t>(cursor[static_cast<std::size_t>(r)]++)] = c;
            adj.neighbor[static_cast<std::size_t>(cursor[static_cast<std::size_t>(c)]++)] = r;
        }
    }

    // Sort and unique each list. A full pattern supplies every edge twice, and
    // an upper triangle once, so duplicates are normal rather than an error.
    std::vector<int> compacted;
    compacted.reserve(adj.neighbor.size());
    std::vector<int> new_offset(static_cast<std::size_t>(n) + 1, 0);
    for (int r = 0; r < n; ++r) {
        const auto begin = adj.neighbor.begin() + adj.offset[static_cast<std::size_t>(r)];
        const auto end = adj.neighbor.begin() + adj.offset[static_cast<std::size_t>(r) + 1];
        std::sort(begin, end);
        const auto last = std::unique(begin, end);
        new_offset[static_cast<std::size_t>(r)] = static_cast<int>(compacted.size());
        compacted.insert(compacted.end(), begin, last);
    }
    new_offset[static_cast<std::size_t>(n)] = static_cast<int>(compacted.size());
    adj.offset = std::move(new_offset);
    adj.neighbor = std::move(compacted);
    return adj;
}

/// One BFS over the component containing `root`, returning the visiting order
/// and the number of levels (the rooted level structure's depth).
///
/// Neighbours are taken in increasing degree, ties by index, so the result is
/// deterministic -- which matters because this project compares results
/// bitwise elsewhere and a nondeterministic ordering would make a factorization
/// nondeterministic too.
int bfs_component(const Adjacency& adj, int root, std::vector<char>& seen,
                  std::vector<int>& order, std::vector<int>& level) {
    const std::size_t first = order.size();
    seen[static_cast<std::size_t>(root)] = 1;
    order.push_back(root);
    level[static_cast<std::size_t>(root)] = 0;
    int depth = 0;

    for (std::size_t head = first; head < order.size(); ++head) {
        const int v = order[head];
        const int lv = level[static_cast<std::size_t>(v)];
        const std::size_t begin = order.size();
        for (int k = adj.offset[static_cast<std::size_t>(v)];
             k < adj.offset[static_cast<std::size_t>(v) + 1]; ++k) {
            const int w = adj.neighbor[static_cast<std::size_t>(k)];
            if (seen[static_cast<std::size_t>(w)]) continue;
            seen[static_cast<std::size_t>(w)] = 1;
            level[static_cast<std::size_t>(w)] = lv + 1;
            depth = std::max(depth, lv + 1);
            order.push_back(w);
        }
        // Within one vertex's newly discovered neighbours, visit lower degree
        // first. This is Cuthill-McKee's rule and it is what makes the
        // bandwidth reduction work.
        std::sort(order.begin() + static_cast<std::ptrdiff_t>(begin), order.end(),
                  [&](int x, int y) {
                      const int dx = adj.degree(x), dy = adj.degree(y);
                      return dx != dy ? dx < dy : x < y;
                  });
    }
    return depth;
}

/// A pseudo-peripheral vertex of the component containing `start`, by the
/// rooted-level-structure heuristic of George & Liu: BFS, take a lowest-degree
/// vertex of the last level, repeat while the depth keeps growing.
///
/// RCM's quality depends materially on where it starts, and a minimum-degree
/// start can be much worse than this. Capped at a few passes because the
/// heuristic can otherwise oscillate.
int pseudo_peripheral(const Adjacency& adj, int start, std::vector<char>& scratch_seen,
                      std::vector<int>& scratch_order, std::vector<int>& scratch_level) {
    int best = start;
    int best_depth = -1;
    for (int pass = 0; pass < 5; ++pass) {
        std::fill(scratch_seen.begin(), scratch_seen.end(), 0);
        scratch_order.clear();
        const int depth = bfs_component(adj, best, scratch_seen, scratch_order, scratch_level);
        if (depth <= best_depth) break;
        best_depth = depth;

        // A lowest-degree vertex in the deepest level, ties by index.
        int candidate = best;
        int candidate_degree = -1;
        for (int v : scratch_order) {
            if (scratch_level[static_cast<std::size_t>(v)] != depth) continue;
            const int d = adj.degree(v);
            if (candidate_degree < 0 || d < candidate_degree) {
                candidate = v;
                candidate_degree = d;
            }
        }
        if (candidate == best) break;
        best = candidate;
    }
    return best;
}

std::vector<int> reverse_cuthill_mckee(const Adjacency& adj, int n) {
    std::vector<char> seen(static_cast<std::size_t>(n), 0);
    std::vector<int> level(static_cast<std::size_t>(n), 0);
    std::vector<int> order;
    order.reserve(static_cast<std::size_t>(n));

    std::vector<char> scratch_seen(static_cast<std::size_t>(n), 0);
    std::vector<int> scratch_order;
    scratch_order.reserve(static_cast<std::size_t>(n));
    std::vector<int> scratch_level(static_cast<std::size_t>(n), 0);

    // Components are taken in order of their lowest-numbered vertex, so the
    // whole thing is deterministic. Isolated vertices -- a voltage port's
    // constraint row is one -- are just components of size 1 and need no
    // special case.
    for (int s = 0; s < n; ++s) {
        if (seen[static_cast<std::size_t>(s)]) continue;
        const int root = pseudo_peripheral(adj, s, scratch_seen, scratch_order, scratch_level);
        bfs_component(adj, root, seen, order, level);
    }

    std::reverse(order.begin(), order.end());
    return order;
}

}  // namespace

Permutation compute_ordering(const SparsityPattern& pattern, Ordering ordering) {
    if (pattern.rows != pattern.cols) {
        throw std::invalid_argument("compute_ordering: the pattern must be square");
    }
    if (ordering == Ordering::Natural) {
        return Permutation::identity(pattern.rows);
    }
    if (ordering == Ordering::ApproximateMinimumDegree) {
        throw std::invalid_argument(
            "compute_ordering: the approximate-minimum-degree ordering is not implemented yet "
            "(step 3 of docs/SOLVER_PLAN.md Sec. 9). Refusing rather than substituting another "
            "ordering, which would report a fill figure belonging to something else.");
    }

    const Adjacency adj = build_adjacency(pattern);
    Permutation p;
    p.perm = reverse_cuthill_mckee(adj, pattern.rows);
    p.iperm.assign(p.perm.size(), 0);
    for (int i = 0; i < static_cast<int>(p.perm.size()); ++i) {
        p.iperm[static_cast<std::size_t>(p.perm[static_cast<std::size_t>(i)])] = i;
    }
    return p;
}

SparsityPattern permute_pattern(const SparsityPattern& pattern, const Permutation& p) {
    if (p.size() != pattern.rows) {
        throw std::invalid_argument(
            "permute_pattern: the permutation is sized " + std::to_string(p.size()) +
            " but the pattern has " + std::to_string(pattern.rows) + " rows");
    }

    SparsityPattern out;
    out.rows = pattern.rows;
    out.cols = pattern.cols;
    out.upper_only = pattern.upper_only;
    out.row_ptr.assign(static_cast<std::size_t>(out.rows) + 1, 0);

    // Entry (i, j) of the result is entry (perm[i], perm[j]) of the input, i.e.
    // an input entry (r, c) lands at (iperm[r], iperm[c]).
    const auto place = [&](int r, int c, int& row, int& col) {
        row = p.iperm[static_cast<std::size_t>(r)];
        col = p.iperm[static_cast<std::size_t>(c)];
        // A permutation moves entries across the diagonal, so an upper-only
        // pattern has to be normalised back or it would stop being upper-only.
        if (out.upper_only && col < row) std::swap(row, col);
    };

    for (int r = 0; r < pattern.rows; ++r) {
        for (int k = pattern.row_ptr[static_cast<std::size_t>(r)];
             k < pattern.row_ptr[static_cast<std::size_t>(r) + 1]; ++k) {
            int row = 0, col = 0;
            place(r, pattern.col_index[static_cast<std::size_t>(k)], row, col);
            ++out.row_ptr[static_cast<std::size_t>(row) + 1];
        }
    }
    for (int r = 0; r < out.rows; ++r) {
        out.row_ptr[static_cast<std::size_t>(r) + 1] += out.row_ptr[static_cast<std::size_t>(r)];
    }

    std::vector<int> cursor(out.row_ptr.begin(), out.row_ptr.end() - 1);
    out.col_index.assign(static_cast<std::size_t>(out.row_ptr.back()), 0);
    for (int r = 0; r < pattern.rows; ++r) {
        for (int k = pattern.row_ptr[static_cast<std::size_t>(r)];
             k < pattern.row_ptr[static_cast<std::size_t>(r) + 1]; ++k) {
            int row = 0, col = 0;
            place(r, pattern.col_index[static_cast<std::size_t>(k)], row, col);
            out.col_index[static_cast<std::size_t>(cursor[static_cast<std::size_t>(row)]++)] = col;
        }
    }
    for (int r = 0; r < out.rows; ++r) {
        std::sort(out.col_index.begin() + out.row_ptr[static_cast<std::size_t>(r)],
                  out.col_index.begin() + out.row_ptr[static_cast<std::size_t>(r) + 1]);
    }
    return out;
}

BandwidthStats bandwidth_stats(const SparsityPattern& pattern) {
    BandwidthStats s;
    for (int r = 0; r < pattern.rows; ++r) {
        const int begin = pattern.row_ptr[static_cast<std::size_t>(r)];
        const int end = pattern.row_ptr[static_cast<std::size_t>(r) + 1];
        if (begin == end) continue;
        int furthest = 0;
        for (int k = begin; k < end; ++k) {
            furthest = std::max(furthest, std::abs(pattern.col_index[static_cast<std::size_t>(k)] - r));
        }
        s.bandwidth = std::max(s.bandwidth, furthest);
        s.profile += furthest;
    }
    s.mean_bandwidth = pattern.rows > 0 ? static_cast<double>(s.profile) / pattern.rows : 0.0;
    return s;
}

}  // namespace aphi_solver
