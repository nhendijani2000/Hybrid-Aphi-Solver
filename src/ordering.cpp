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

SymmetricAdjacency SymmetricAdjacency::build(const SparsityPattern& pattern) {
    const int n = pattern.rows;
    SymmetricAdjacency adj;
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

namespace {

/// One BFS over the component containing `root`, returning the visiting order
/// and the number of levels (the rooted level structure's depth).
///
/// Neighbours are taken in increasing degree, ties by index, so the result is
/// deterministic -- which matters because this project compares results
/// bitwise elsewhere and a nondeterministic ordering would make a factorization
/// nondeterministic too.
int bfs_component(const SymmetricAdjacency& adj, int root, std::vector<char>& seen,
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
int pseudo_peripheral(const SymmetricAdjacency& adj, int start, std::vector<char>& scratch_seen,
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

/// Approximate minimum degree, after Amestoy, Davis & Duff (SIAM J. Matrix
/// Anal. Appl. 17(4):886-905, 1996; Algorithm 837, ACM TOMS 2004).
///
/// Minimum degree orders by repeatedly eliminating a vertex of least degree in
/// the *elimination* graph -- the graph that gains a clique among a vertex's
/// neighbours each time one is removed. Doing that literally is ruinous,
/// because those cliques are dense. The quotient graph avoids materialising
/// them: an eliminated vertex becomes an **element** holding the list of
/// variables it connects, and a variable's neighbourhood is then its remaining
/// variables plus the union of its elements' lists. Nothing is ever expanded.
///
/// The "approximate" is the degree. The true external degree of a variable
/// costs a set union over its elements to compute; AMD instead uses a bound
/// that is exact when the variable has at most two elements and an upper bound
/// otherwise, which is cheap and in practice orders about as well.
///
/// **What this implementation leaves out of the paper**, deliberately, per
/// Stage 1's correctness-over-speed brief:
///
///  - ~~Supervariables~~ -- IMPLEMENTED 2026-09-29. Was listed here as "mostly
///    a speed optimisation", which was wrong: it cut case 02 fill 211.5 M to
///    173.0 M and the factorization 550 s to 335 s. Eliminating
///    indistinguishable variables separately updates degrees from stale data,
///    so the ordering itself degrades. (Old note follows.)
///  - Supervariables (detecting indistinguishable variables by hashing and
///    eliminating them together). Mostly a speed optimisation, and FEM meshes
///    do have many indistinguishable vertices, so this is where the time goes
///    if it ever becomes a problem.
///  - **Aggressive absorption** (absorbing an element as soon as its list is
///    contained in the new one, not only when it is adjacent to the pivot).
///    Also mostly speed.
///  - The paper's packed integer workspace with garbage collection, in favour
///    of per-vertex vectors. Costs memory and some locality, and is far easier
///    to read and to be sure of.
///
/// Leaving these out can only make the ordering worse, never wrong: **any**
/// permutation yields a correct factorization, so a defect here shows up as
/// fill, which `predict_factor_size` measures exactly. That is why this is
/// implementable with confidence at this stage.
std::vector<int> approximate_minimum_degree(const SymmetricAdjacency& adj, int n) {
    enum class State : char { Variable, Element, Absorbed };

    std::vector<State> state(static_cast<std::size_t>(n), State::Variable);
    std::vector<std::vector<int>> vars(static_cast<std::size_t>(n));  // A_i: adjacent variables
    std::vector<std::vector<int>> elems(static_cast<std::size_t>(n)); // E_i: adjacent elements
    std::vector<std::vector<int>> members(static_cast<std::size_t>(n));  // L_e for an element
    std::vector<int> degree(static_cast<std::size_t>(n), 0);

    // --- supervariables ----------------------------------------------------
    // `nv[i]` is how many ORIGINAL variables the principal variable i stands
    // for, and `merged[i]` names them. Two variables are indistinguishable when
    // their adjacency is identical; eliminating them together is not merely
    // faster, it orders BETTER, because eliminating them one at a time updates
    // every degree using stale information about the other. Finite-element
    // meshes are full of indistinguishable vertices, which is why leaving this
    // out cost a factor of two in fill here -- see SOLVER_PLAN.md Sec. 15.
    std::vector<int> nv(static_cast<std::size_t>(n), 1);
    std::vector<std::vector<int>> merged(static_cast<std::size_t>(n));
    // Weighted |L_e|, in original variables. Set when an element is created and
    // valid for its lifetime: a merge moves weight between two members that are
    // BOTH in the element (indistinguishable variables share their elements), so
    // the total does not move; and a variable's elimination absorbs every
    // element it belonged to.
    std::vector<int> elem_weight(static_cast<std::size_t>(n), 0);

    for (int i = 0; i < n; ++i) {
        vars[static_cast<std::size_t>(i)].assign(
            adj.neighbor.begin() + adj.offset[static_cast<std::size_t>(i)],
            adj.neighbor.begin() + adj.offset[static_cast<std::size_t>(i) + 1]);
        degree[static_cast<std::size_t>(i)] = static_cast<int>(vars[static_cast<std::size_t>(i)].size());
    }

    // Degree buckets: a doubly linked list per degree, so finding the minimum is
    // a scan from the last known minimum rather than over every variable. A
    // linear scan each step would be O(n^2), which at 27000 unknowns is 7e8
    // operations before any real work.
    std::vector<int> head(static_cast<std::size_t>(n) + 1, -1);
    std::vector<int> next(static_cast<std::size_t>(n), -1);
    std::vector<int> prev(static_cast<std::size_t>(n), -1);
    const auto bucket_insert = [&](int i) {
        const int d = std::min(degree[static_cast<std::size_t>(i)], n);
        next[static_cast<std::size_t>(i)] = head[static_cast<std::size_t>(d)];
        prev[static_cast<std::size_t>(i)] = -1;
        if (head[static_cast<std::size_t>(d)] >= 0) {
            prev[static_cast<std::size_t>(head[static_cast<std::size_t>(d)])] = i;
        }
        head[static_cast<std::size_t>(d)] = i;
    };
    const auto bucket_remove = [&](int i) {
        const int d = std::min(degree[static_cast<std::size_t>(i)], n);
        const int p = prev[static_cast<std::size_t>(i)];
        const int q = next[static_cast<std::size_t>(i)];
        if (p >= 0) {
            next[static_cast<std::size_t>(p)] = q;
        } else if (head[static_cast<std::size_t>(d)] == i) {
            head[static_cast<std::size_t>(d)] = q;
        }
        if (q >= 0) prev[static_cast<std::size_t>(q)] = p;
        prev[static_cast<std::size_t>(i)] = -1;
        next[static_cast<std::size_t>(i)] = -1;
    };
    for (int i = 0; i < n; ++i) bucket_insert(i);

    std::vector<int> order;
    order.reserve(static_cast<std::size_t>(n));
    std::vector<int> stamp(static_cast<std::size_t>(n), -1);  // membership of the current L_p
    std::vector<int> external(static_cast<std::size_t>(n), 0);  // |L_e \ L_p| per element
    std::vector<int> touched_elements;
    std::vector<int> pivot_list;
    // Scratch for supervariable detection, allocated once.
    std::vector<long long> hash_of(static_cast<std::size_t>(n), 0);
    std::vector<int> hash_head(static_cast<std::size_t>(n), -1);
    std::vector<int> hash_next(static_cast<std::size_t>(n), -1);
    std::vector<int> hash_used;
    std::vector<int> mark(static_cast<std::size_t>(n), -1);
    int min_degree = 0;
    int eliminated = 0;
    int step = 0;

    while (eliminated < n) {
        // --- pick a variable of least (approximate) degree --------------
        while (min_degree <= n && head[static_cast<std::size_t>(min_degree)] < 0) ++min_degree;
        if (min_degree > n) break;  // nothing left, which happens only if n == 0
        const int p = head[static_cast<std::size_t>(min_degree)];
        bucket_remove(p);

        // A supervariable is eliminated as a block: it and everything merged
        // into it leave together, in that order.
        order.push_back(p);
        for (int q : merged[static_cast<std::size_t>(p)]) order.push_back(q);
        eliminated += nv[static_cast<std::size_t>(p)];
        const int k = step++;

        // --- L_p: the pivot's neighbourhood, principal variables only ---
        pivot_list.clear();
        for (int j : vars[static_cast<std::size_t>(p)]) {
            if (state[static_cast<std::size_t>(j)] != State::Variable) continue;
            if (stamp[static_cast<std::size_t>(j)] == k) continue;
            stamp[static_cast<std::size_t>(j)] = k;
            pivot_list.push_back(j);
        }
        for (int e : elems[static_cast<std::size_t>(p)]) {
            if (state[static_cast<std::size_t>(e)] != State::Element) continue;
            for (int j : members[static_cast<std::size_t>(e)]) {
                if (state[static_cast<std::size_t>(j)] != State::Variable) continue;
                if (j == p || stamp[static_cast<std::size_t>(j)] == k) continue;
                stamp[static_cast<std::size_t>(j)] = k;
                pivot_list.push_back(j);
            }
        }

        // --- the pivot's elements are absorbed into it ------------------
        // Their lists are subsets of L_p now, so keeping them would only make
        // every later neighbourhood scan longer.
        for (int e : elems[static_cast<std::size_t>(p)]) {
            if (state[static_cast<std::size_t>(e)] != State::Element) continue;
            state[static_cast<std::size_t>(e)] = State::Absorbed;
            members[static_cast<std::size_t>(e)].clear();
            members[static_cast<std::size_t>(e)].shrink_to_fit();
        }

        // --- prune each neighbour's lists ------------------------------
        // A_i loses p and anything now reachable through element p; E_i loses
        // absorbed elements and gains p.
        for (int i : pivot_list) {
            std::vector<int>& ai = vars[static_cast<std::size_t>(i)];
            std::size_t w = 0;
            for (std::size_t t = 0; t < ai.size(); ++t) {
                const int j = ai[t];
                if (j == p) continue;
                if (state[static_cast<std::size_t>(j)] != State::Variable) continue;
                if (stamp[static_cast<std::size_t>(j)] == k) continue;  // in L_p already
                ai[w++] = j;
            }
            ai.resize(w);

            std::vector<int>& ei = elems[static_cast<std::size_t>(i)];
            std::size_t v = 0;
            for (std::size_t t = 0; t < ei.size(); ++t) {
                const int e = ei[t];
                if (state[static_cast<std::size_t>(e)] != State::Element) continue;
                ei[v++] = e;
            }
            ei.resize(v);
            ei.push_back(p);
        }

        // --- |L_e \ L_p| for every element still adjacent to L_p --------
        // Weighted, in original variables. Start each at the element's weight,
        // then subtract nv[i] once per member that is also in L_p. Every i in
        // L_p belongs to every element in its own E_i, so one pass suffices.
        touched_elements.clear();
        for (int i : pivot_list) {
            for (int e : elems[static_cast<std::size_t>(i)]) {
                if (e == p) continue;
                if (stamp[static_cast<std::size_t>(e)] != -2 - k) {
                    stamp[static_cast<std::size_t>(e)] = -2 - k;
                    external[static_cast<std::size_t>(e)] = elem_weight[static_cast<std::size_t>(e)];
                    touched_elements.push_back(e);
                }
                external[static_cast<std::size_t>(e)] -= nv[static_cast<std::size_t>(i)];
            }
        }

        // --- the approximate degree ------------------------------------
        // Every count is weighted: a degree is a number of ORIGINAL variables,
        // not of supervariables, or the bucket order would compare unlike
        // things.
        long long lp_weight = 0;
        for (int i : pivot_list) lp_weight += nv[static_cast<std::size_t>(i)];

        for (int i : pivot_list) {
            long long sum = 0;
            for (int e : elems[static_cast<std::size_t>(i)]) {
                if (e == p) continue;
                sum += std::max(0, external[static_cast<std::size_t>(e)]);
            }
            long long ai_weight = 0;
            for (int j : vars[static_cast<std::size_t>(i)]) {
                if (state[static_cast<std::size_t>(j)] == State::Variable) {
                    ai_weight += nv[static_cast<std::size_t>(j)];
                }
            }
            const long long approx = ai_weight + (lp_weight - nv[static_cast<std::size_t>(i)]) + sum;

            // Three bounds, whichever is tightest. The first two are exact
            // facts about any elimination graph; the third is AMD's estimate.
            const long long remaining = static_cast<long long>(n) - eliminated;
            const long long grew = static_cast<long long>(degree[static_cast<std::size_t>(i)]) +
                                   lp_weight - nv[static_cast<std::size_t>(i)];
            const long long d = std::min(remaining, std::min(grew, approx));

            bucket_remove(i);
            degree[static_cast<std::size_t>(i)] = static_cast<int>(std::max(0LL, d));
            bucket_insert(i);
            min_degree = std::min(min_degree, degree[static_cast<std::size_t>(i)]);
        }
        for (int e : touched_elements) stamp[static_cast<std::size_t>(e)] = -1;

        // --- the pivot becomes an element ------------------------------
        state[static_cast<std::size_t>(p)] = State::Element;
        members[static_cast<std::size_t>(p)] = pivot_list;
        elem_weight[static_cast<std::size_t>(p)] = static_cast<int>(lp_weight);
        vars[static_cast<std::size_t>(p)].clear();
        vars[static_cast<std::size_t>(p)].shrink_to_fit();
        elems[static_cast<std::size_t>(p)].clear();
        elems[static_cast<std::size_t>(p)].shrink_to_fit();

        // --- supervariable detection within L_p -------------------------
        // Only variables in L_p can have become indistinguishable this step:
        // nothing outside it had its adjacency changed. Hash (A_i union E_i),
        // then compare exactly within a hash bucket -- the hash only narrows
        // the candidates, it never decides a merge on its own.
        hash_used.clear();
        for (int i : pivot_list) {
            if (state[static_cast<std::size_t>(i)] != State::Variable) continue;
            long long h = 0;
            for (int j : vars[static_cast<std::size_t>(i)]) h += j;
            for (int e : elems[static_cast<std::size_t>(i)]) h += e;
            h %= n;
            hash_of[static_cast<std::size_t>(i)] = h;
            const std::size_t hb = static_cast<std::size_t>(h);
            if (hash_head[hb] < 0) hash_used.push_back(static_cast<int>(h));
            hash_next[static_cast<std::size_t>(i)] = hash_head[hb];
            hash_head[hb] = i;
        }

        for (int h : hash_used) {
            for (int i = hash_head[static_cast<std::size_t>(h)]; i >= 0;
                 i = hash_next[static_cast<std::size_t>(i)]) {
                if (state[static_cast<std::size_t>(i)] != State::Variable) continue;

                // Mark i's neighbourhood once, then test each later candidate
                // against it. Comparing counts first rejects most pairs before
                // any set work.
                for (int j : vars[static_cast<std::size_t>(i)]) mark[static_cast<std::size_t>(j)] = i;
                for (int e : elems[static_cast<std::size_t>(i)]) mark[static_cast<std::size_t>(e)] = i;

                for (int j = hash_next[static_cast<std::size_t>(i)]; j >= 0;
                     j = hash_next[static_cast<std::size_t>(j)]) {
                    if (state[static_cast<std::size_t>(j)] != State::Variable) continue;
                    if (vars[static_cast<std::size_t>(j)].size() !=
                            vars[static_cast<std::size_t>(i)].size() ||
                        elems[static_cast<std::size_t>(j)].size() !=
                            elems[static_cast<std::size_t>(i)].size()) {
                        continue;
                    }
                    bool same = true;
                    for (int q : vars[static_cast<std::size_t>(j)]) {
                        if (mark[static_cast<std::size_t>(q)] != i) { same = false; break; }
                    }
                    if (same) {
                        for (int e : elems[static_cast<std::size_t>(j)]) {
                            if (mark[static_cast<std::size_t>(e)] != i) { same = false; break; }
                        }
                    }
                    if (!same) continue;

                    // Absorb j into i. j leaves the graph; every list that
                    // names it skips it on the State check, exactly as it
                    // already does for eliminated variables.
                    bucket_remove(j);
                    nv[static_cast<std::size_t>(i)] += nv[static_cast<std::size_t>(j)];
                    merged[static_cast<std::size_t>(i)].push_back(j);
                    for (int q : merged[static_cast<std::size_t>(j)]) {
                        merged[static_cast<std::size_t>(i)].push_back(q);
                    }
                    merged[static_cast<std::size_t>(j)].clear();
                    merged[static_cast<std::size_t>(j)].shrink_to_fit();
                    nv[static_cast<std::size_t>(j)] = 0;
                    state[static_cast<std::size_t>(j)] = State::Absorbed;
                    vars[static_cast<std::size_t>(j)].clear();
                    vars[static_cast<std::size_t>(j)].shrink_to_fit();
                    elems[static_cast<std::size_t>(j)].clear();
                    elems[static_cast<std::size_t>(j)].shrink_to_fit();
                }
            }
            hash_head[static_cast<std::size_t>(h)] = -1;
        }
    }

    // Any variable never selected -- which cannot happen, but a silent
    // truncation here would shrink the system without saying so.
    if (static_cast<int>(order.size()) != n) {
        throw std::logic_error("approximate_minimum_degree: ordered " +
                               std::to_string(order.size()) + " of " + std::to_string(n) +
                               " unknowns. Every variable must be eliminated exactly once.");
    }
    return order;
}

std::vector<int> reverse_cuthill_mckee(const SymmetricAdjacency& adj, int n) {
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
    const SymmetricAdjacency adj = SymmetricAdjacency::build(pattern);
    Permutation p;
    p.perm = ordering == Ordering::ApproximateMinimumDegree
                 ? approximate_minimum_degree(adj, pattern.rows)
                 : reverse_cuthill_mckee(adj, pattern.rows);
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
