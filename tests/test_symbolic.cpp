// Tests for the symbolic pass -- step 2 of `docs/SOLVER_PLAN.md` §9.
//
// The strongest test here is an INDEPENDENT reference: a naive dense symbolic
// elimination, which simulates fill by explicit set union in O(n^3) and shares
// no code with the elimination-tree path. On small patterns the two must agree
// on nnz(L) exactly. Every clever thing about the real implementation -- path
// compression, the reachable-set walk, the marking -- is a chance to be subtly
// wrong in a way that still produces a plausible number, and this is what
// catches that.
//
// Alongside it, three patterns whose factor is known by hand: diagonal (no
// fill), tridiagonal (no fill), and dense (complete fill).

#include <algorithm>
#include <iostream>
#include <numeric>
#include <random>
#include <set>
#include <string>
#include <vector>

#include "aphi_solver/gmsh_reader.hpp"
#include "aphi_solver/problem_binding.hpp"
#include "aphi_solver/symbolic.hpp"

using namespace aphi_solver;

namespace {

int g_checks = 0;
int g_failures = 0;

void check(bool condition, const std::string& what) {
    ++g_checks;
    if (!condition) {
        ++g_failures;
        std::cerr << "FAIL: " << what << "\n";
    }
}

SparsityPattern pattern_from_edges(int n, const std::vector<std::pair<int, int>>& edges) {
    std::vector<std::set<int>> rows(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) rows[static_cast<std::size_t>(i)].insert(i);
    for (const auto& e : edges) {
        rows[static_cast<std::size_t>(e.first)].insert(e.second);
        rows[static_cast<std::size_t>(e.second)].insert(e.first);
    }
    SparsityPattern p;
    p.rows = n;
    p.cols = n;
    p.row_ptr.assign(static_cast<std::size_t>(n) + 1, 0);
    for (int i = 0; i < n; ++i) {
        p.row_ptr[static_cast<std::size_t>(i) + 1] =
            p.row_ptr[static_cast<std::size_t>(i)] +
            static_cast<int>(rows[static_cast<std::size_t>(i)].size());
    }
    for (int i = 0; i < n; ++i) {
        for (int c : rows[static_cast<std::size_t>(i)]) p.col_index.push_back(c);
    }
    return p;
}

/// nnz(L) by brute force: eliminate one column at a time and add the fill
/// explicitly. O(n^3) and obviously correct, which is the point -- it shares no
/// code with `elimination_tree` or the reachable-set walk.
std::size_t dense_reference_nnz(const SparsityPattern& p) {
    const int n = p.rows;
    std::vector<std::set<int>> adj(static_cast<std::size_t>(n));
    for (int r = 0; r < n; ++r) {
        for (int k = p.row_ptr[static_cast<std::size_t>(r)];
             k < p.row_ptr[static_cast<std::size_t>(r) + 1]; ++k) {
            const int c = p.col_index[static_cast<std::size_t>(k)];
            if (c == r) continue;
            adj[static_cast<std::size_t>(r)].insert(c);
            adj[static_cast<std::size_t>(c)].insert(r);
        }
    }
    std::size_t nnz = 0;
    for (int k = 0; k < n; ++k) {
        // Column k of L: k itself, plus every remaining neighbour after k.
        std::vector<int> later;
        for (int v : adj[static_cast<std::size_t>(k)]) {
            if (v > k) later.push_back(v);
        }
        nnz += 1 + later.size();
        // Eliminating k makes its later neighbours a clique.
        for (std::size_t a = 0; a < later.size(); ++a) {
            for (std::size_t b = a + 1; b < later.size(); ++b) {
                adj[static_cast<std::size_t>(later[a])].insert(later[b]);
                adj[static_cast<std::size_t>(later[b])].insert(later[a]);
            }
        }
    }
    return nnz;
}

// ---------------------------------------------------------------------------

void test_hand_computed_factors() {
    // Diagonal: no off-diagonal anything, so no fill and no tree.
    const SparsityPattern diag = pattern_from_edges(5, {});
    const std::vector<int> dparent = elimination_tree(diag);
    check(std::count(dparent.begin(), dparent.end(), -1) == 5,
          "a diagonal matrix's elimination tree is five isolated roots");
    const FactorSize dsize = predict_factor_size(diag, Permutation::identity(5));
    check(dsize.nnz == 5, "and L is just the diagonal: nnz(L) = 5");
    check(dsize.num_roots == 5, "with five components");

    // Tridiagonal: L is bidiagonal, so nnz(L) = 2n - 1, and the tree is a path.
    const int n = 12;
    std::vector<std::pair<int, int>> chain;
    for (int i = 0; i + 1 < n; ++i) chain.push_back({i, i + 1});
    const SparsityPattern tri = pattern_from_edges(n, chain);
    const std::vector<int> tparent = elimination_tree(tri);
    bool is_path = true;
    for (int j = 0; j + 1 < n; ++j) {
        if (tparent[static_cast<std::size_t>(j)] != j + 1) is_path = false;
    }
    check(is_path && tparent[static_cast<std::size_t>(n - 1)] == -1,
          "a tridiagonal matrix's elimination tree is a path 0->1->...->n-1");
    const FactorSize tsize = predict_factor_size(tri, Permutation::identity(n));
    check(tsize.nnz == static_cast<std::size_t>(2 * n - 1),
          "and L is bidiagonal: nnz(L) = 2n - 1 = " + std::to_string(2 * n - 1) + ", got " +
              std::to_string(tsize.nnz));
    check(tsize.num_roots == 1, "one component");

    // Dense: complete fill, nnz(L) = n(n+1)/2, tree a path.
    const int m = 9;
    std::vector<std::pair<int, int>> all;
    for (int i = 0; i < m; ++i) {
        for (int j = i + 1; j < m; ++j) all.push_back({i, j});
    }
    const SparsityPattern dense = pattern_from_edges(m, all);
    const FactorSize msize = predict_factor_size(dense, Permutation::identity(m));
    check(msize.nnz == static_cast<std::size_t>(m * (m + 1) / 2),
          "a dense matrix fills completely: nnz(L) = n(n+1)/2 = " +
              std::to_string(m * (m + 1) / 2) + ", got " + std::to_string(msize.nnz));
}

void test_against_the_dense_reference() {
    std::mt19937 rng(20260926);
    double worst_ratio = 0.0;
    int cases = 0, agreed = 0;

    for (int trial = 0; trial < 40; ++trial) {
        const int n = 8 + static_cast<int>(rng() % 40);
        const double density = 0.03 + 0.12 * (static_cast<double>(rng() % 100) / 100.0);
        std::vector<std::pair<int, int>> edges;
        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                if (static_cast<double>(rng() % 10000) / 10000.0 < density) edges.push_back({i, j});
            }
        }
        const SparsityPattern p = pattern_from_edges(n, edges);

        // Both orderings, because an ordering bug and a symbolic bug would
        // otherwise be able to cancel.
        for (const Ordering o : {Ordering::Natural, Ordering::ReverseCuthillMcKee}) {
            const Permutation perm = compute_ordering(p, o);
            const SparsityPattern permuted = permute_pattern(p, perm);
            const FactorSize mine = predict_factor_size(p, perm);
            const std::size_t reference = dense_reference_nnz(permuted);
            ++cases;
            if (mine.nnz == reference) ++agreed;
            if (reference > 0) {
                worst_ratio = std::max(worst_ratio,
                                       std::abs(static_cast<double>(mine.nnz) -
                                                static_cast<double>(reference)) /
                                           static_cast<double>(reference));
            }
        }
    }
    std::cout << "  dense reference: " << agreed << "/" << cases
              << " random patterns agree exactly on nnz(L)\n";
    check(agreed == cases,
          "nnz(L) matches an independent dense symbolic elimination on every random pattern (" +
              std::to_string(agreed) + " of " + std::to_string(cases) + ", worst relative "
              "disagreement " + std::to_string(worst_ratio) + ")");
}

void test_tree_and_pattern_invariants() {
    const int n = 60;
    std::vector<std::pair<int, int>> edges;
    for (int i = 0; i + 1 < n; ++i) edges.push_back({i, i + 1});
    for (int i = 0; i + 9 < n; ++i) edges.push_back({i, i + 9});
    const SparsityPattern p = pattern_from_edges(n, edges);
    const SolverAnalysis a = analyze(p, Ordering::ReverseCuthillMcKee);

    // The tree points strictly upwards, which is what makes the reachable-set
    // walk terminate at all.
    bool upward = true;
    for (int j = 0; j < n; ++j) {
        const int par = a.parent[static_cast<std::size_t>(j)];
        if (par != -1 && par <= j) upward = false;
    }
    check(upward, "parent[j] is either -1 or strictly greater than j");

    check(a.predicted_nnz == a.factor.nnz(),
          "the predicted nnz(L) equals the pattern actually built (" +
              std::to_string(a.predicted_nnz) + " vs " + std::to_string(a.factor.nnz()) +
              ") -- the counting and building passes are independent code paths");

    // L is lower triangular, each row ascending, and the diagonal is present:
    // all three are relied on by a numeric factorization.
    bool lower = true, ascending = true, has_diagonal = true;
    for (int r = 0; r < n; ++r) {
        const int begin = a.factor.row_ptr[static_cast<std::size_t>(r)];
        const int end = a.factor.row_ptr[static_cast<std::size_t>(r) + 1];
        if (end <= begin) {
            has_diagonal = false;
            continue;
        }
        if (a.factor.col_index[static_cast<std::size_t>(end - 1)] != r) has_diagonal = false;
        for (int k = begin; k < end; ++k) {
            if (a.factor.col_index[static_cast<std::size_t>(k)] > r) lower = false;
            if (k > begin && a.factor.col_index[static_cast<std::size_t>(k)] <=
                                 a.factor.col_index[static_cast<std::size_t>(k - 1)]) {
                ascending = false;
            }
        }
    }
    check(lower, "L is lower triangular");
    check(ascending, "each row of L is strictly ascending, so it can be binary-searched");
    check(has_diagonal, "and every row has its diagonal, which is the last entry");

    // Fill only adds. L must contain the lower triangle of the permuted matrix.
    bool contains_a = true;
    for (int r = 0; r < n; ++r) {
        for (int k = a.permuted.row_ptr[static_cast<std::size_t>(r)];
             k < a.permuted.row_ptr[static_cast<std::size_t>(r) + 1]; ++k) {
            const int c = a.permuted.col_index[static_cast<std::size_t>(k)];
            if (c > r) continue;
            if (a.factor.find_slot(r, c) < 0) contains_a = false;
        }
    }
    check(contains_a, "L contains every entry of the permuted matrix's lower triangle -- "
                      "factorization adds fill, it never removes an entry");

    std::size_t lower_a = 0;
    for (int r = 0; r < n; ++r) {
        for (int k = a.permuted.row_ptr[static_cast<std::size_t>(r)];
             k < a.permuted.row_ptr[static_cast<std::size_t>(r) + 1]; ++k) {
            if (a.permuted.col_index[static_cast<std::size_t>(k)] <= r) ++lower_a;
        }
    }
    check(a.predicted_nnz >= lower_a,
          "so nnz(L) is at least nnz(tril(A)) (" + std::to_string(a.predicted_nnz) + " >= " +
              std::to_string(lower_a) + ")");
}

void test_upper_only_pattern_gives_the_same_answer() {
    // An upper-triangle pattern describes the same matrix, so it must give the
    // same tree and the same nnz(L). Without this, the symmetric storage path
    // and the full one could silently disagree about fill.
    const int n = 30;
    std::vector<std::pair<int, int>> edges;
    for (int i = 0; i + 1 < n; ++i) edges.push_back({i, i + 1});
    for (int i = 0; i + 4 < n; ++i) edges.push_back({i, i + 4});
    const SparsityPattern full = pattern_from_edges(n, edges);

    SparsityPattern upper;
    upper.rows = n;
    upper.cols = n;
    upper.upper_only = true;
    upper.row_ptr.assign(static_cast<std::size_t>(n) + 1, 0);
    for (int r = 0; r < n; ++r) {
        int kept = 0;
        for (int k = full.row_ptr[static_cast<std::size_t>(r)];
             k < full.row_ptr[static_cast<std::size_t>(r) + 1]; ++k) {
            if (full.col_index[static_cast<std::size_t>(k)] >= r) ++kept;
        }
        upper.row_ptr[static_cast<std::size_t>(r) + 1] =
            upper.row_ptr[static_cast<std::size_t>(r)] + kept;
    }
    for (int r = 0; r < n; ++r) {
        for (int k = full.row_ptr[static_cast<std::size_t>(r)];
             k < full.row_ptr[static_cast<std::size_t>(r) + 1]; ++k) {
            const int c = full.col_index[static_cast<std::size_t>(k)];
            if (c >= r) upper.col_index.push_back(c);
        }
    }

    check(elimination_tree(full) == elimination_tree(upper),
          "an upper-triangle pattern gives the same elimination tree as the full one");
    const FactorSize a = predict_factor_size(full, Permutation::identity(n));
    const FactorSize b = predict_factor_size(upper, Permutation::identity(n));
    check(a.nnz == b.nnz && a.num_roots == b.num_roots,
          "and the same nnz(L) -- storing one triangle does not change the factorization");
}

void test_budget_is_checked_before_allocating() {
    const int n = 40;
    std::vector<std::pair<int, int>> all;
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) all.push_back({i, j});
    }
    const SparsityPattern dense = pattern_from_edges(n, all);  // fills completely

    bool refused = false;
    std::string message;
    try {
        analyze(dense, Ordering::Natural, 100);
    } catch (const std::invalid_argument& e) {
        refused = true;
        message = e.what();
    }
    check(refused, "a factor over the budget is refused");
    check(message.find(std::to_string(n * (n + 1) / 2)) != std::string::npos,
          "and the message names the predicted nonzero count, so the budget can be set from it");

    bool accepted = true;
    try {
        analyze(dense, Ordering::Natural, 0);
    } catch (...) {
        accepted = false;
    }
    check(accepted, "control: budget 0 means no limit, so the same call succeeds");
}

// The measurement this step exists to produce.
void test_ordering_table_on_real_meshes() {
#ifdef APHI_MESH_DIR
    const std::string dir = APHI_MESH_DIR;
    for (const char* name : {"cylinder_box.msh", "loop_cut.msh"}) {
        Mesh m;
        try {
            m = read_gmsh_msh(dir + "/" + name);
        } catch (const std::exception& e) {
            check(false, std::string(name) + " unreadable: " + e.what());
            continue;
        }
        scale_mesh_to_metres(m, LengthUnit::Millimetre);

        const bool cylinder = std::string(name)[0] == 'c';
        Problem p;
        p.type = AnalysisType::DC;
        p.length_unit = LengthUnit::Millimetre;
        Body cond;
        cond.name = "B1";
        cond.volume = cylinder ? "wire" : "ring";
        cond.sigma = 5.8e7;
        cond.line = 1;
        Body air;
        air.name = "B2";
        air.volume = "air";
        air.sigma = 0.0;
        air.line = 2;
        p.bodies = {cond, air};
        if (cylinder) {
            Port a;
            a.name = "P1";
            a.type = PortType::BoundaryCurrent;
            a.surface = {"wire_bottom"};
            a.amplitude = 1.0;
            a.line = 3;
            Port b;
            b.name = "P2";
            b.type = PortType::BoundaryVoltage;
            b.surface = {"wire_top"};
            b.amplitude = 0.0;
            b.line = 4;
            p.ports = {a, b};
        } else {
            Port a;
            a.name = "P1";
            a.type = PortType::InternalCurrent;
            a.surface = {"loop_cut"};
            a.amplitude = 1.0;
            a.current_direction = Vec3{0.0, 1.0, 0.0};
            a.line = 3;
            p.ports = {a};
        }

        const BoundProblem b = bind_to_mesh(p, m);
        const DofMap d = build_dof_map(b, m);
        const SparsityPattern pattern = build_sparsity(d, b, m);

        // Prediction only, never the pattern: with the natural ordering the
        // factor is far too large to build, which is the whole point.
        const FactorSize nat =
            predict_factor_size(pattern, compute_ordering(pattern, Ordering::Natural));
        const FactorSize rcm =
            predict_factor_size(pattern, compute_ordering(pattern, Ordering::ReverseCuthillMcKee));

        std::cout << "  " << name << ": " << pattern.rows << " unknowns, " << pattern.nnz()
                  << " nonzeros in A\n";
        std::cout << "      natural  nnz(L) = " << nat.nnz << "  (" << nat.mean_row
                  << " per row, " << nat.bytes(pattern.rows) / 1048576 << " MB)\n";
        std::cout << "      rcm      nnz(L) = " << rcm.nnz << "  (" << rcm.mean_row
                  << " per row, " << rcm.bytes(pattern.rows) / 1048576 << " MB)   "
                  << static_cast<double>(nat.nnz) / static_cast<double>(rcm.nnz) << "x less fill\n";

        check(nat.num_roots == rcm.num_roots,
              std::string(name) + ": both orderings see the same number of components");
        check(rcm.nnz < nat.nnz,
              std::string(name) + ": RCM produces LESS fill than the natural ordering (" +
                  std::to_string(rcm.nnz) + " vs " + std::to_string(nat.nnz) +
                  ") -- the control for the whole reordering step");
        check(rcm.nnz >= pattern.nnz() / 2,
              std::string(name) + ": and at least as much as A's own lower triangle");

        // Pinned from measurement, because it is what catches a missing
        // reversal -- the control that step 1 could not test at all.
        //
        // Correct RCM:            1027 (cylinder) and 1049 (loop) per row of L.
        // Plain Cuthill-McKee:    2071 and 2017 -- the reversal dropped.
        //
        // The threshold sits between, with about 1.4x margin either side. No
        // bandwidth or envelope measure can see this difference on a symmetric
        // pattern (SOLVER_PLAN Sec. 13 has the reason), so nnz(L) is the only
        // place it can be caught, and 2x of fill is far too much to leave
        // uncovered.
        check(rcm.mean_row < 1500.0,
              std::string(name) + ": RCM's factor averages under 1500 nonzeros per row "
                                  "(got " + std::to_string(rcm.mean_row) +
                  ") -- plain Cuthill-McKee, without the reversal, gives about 2050");
    }
#endif
}

}  // namespace

int main() {
    test_hand_computed_factors();
    test_against_the_dense_reference();
    test_tree_and_pattern_invariants();
    test_upper_only_pattern_gives_the_same_answer();
    test_budget_is_checked_before_allocating();
    test_ordering_table_on_real_meshes();

    std::cout << (g_checks - g_failures) << "/" << g_checks << " checks passed\n";
    return g_failures == 0 ? 0 : 1;
}
