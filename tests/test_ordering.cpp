// Tests for the reordering step -- step 1 of `docs/SOLVER_PLAN.md` §9.
//
// There is no factorization yet, so nothing here can measure fill. What it can
// measure is what RCM actually reduces (bandwidth and profile) and the
// invariants a permutation has to satisfy before it is safe to permute a matrix
// with one: bijectivity, mutual inversion, and that permuting preserves the
// nonzero count.
//
// The sharpest test is the shuffled chain: a path graph whose labels have been
// scrambled has a terrible natural bandwidth and an optimal one of 1, so RCM
// either recovers it or it does not work.

#include <algorithm>
#include <iostream>
#include <numeric>
#include <random>
#include <set>
#include <string>
#include <vector>

#include "aphi_solver/gmsh_reader.hpp"
#include "aphi_solver/ordering.hpp"
#include "aphi_solver/problem_binding.hpp"

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

/// A pattern from an explicit edge list, symmetric, with the diagonal included
/// -- the shape an assembled matrix has.
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
    p.col_index.reserve(static_cast<std::size_t>(p.row_ptr.back()));
    for (int i = 0; i < n; ++i) {
        for (int c : rows[static_cast<std::size_t>(i)]) p.col_index.push_back(c);
    }
    return p;
}

// ---------------------------------------------------------------------------

void test_keywords() {
    check(std::string(ordering_keyword(Ordering::Natural)) == "natural" &&
              std::string(ordering_keyword(Ordering::ReverseCuthillMcKee)) == "rcm" &&
              std::string(ordering_keyword(Ordering::ApproximateMinimumDegree)) == "amd",
          "every ordering has an input-file keyword");
    Ordering got = Ordering::Natural;
    check(ordering_from_keyword("rcm", got) && got == Ordering::ReverseCuthillMcKee,
          "and the keyword maps back");
    check(!ordering_from_keyword("minimum_degree", got),
          "an unrecognised keyword is refused rather than guessed at");
}

void test_permutation_invariants() {
    const Permutation id = Permutation::identity(5);
    check(id.is_valid(), "the identity is a valid permutation");
    check(id.perm == std::vector<int>({0, 1, 2, 3, 4}), "and it is the identity");
    check(id.inverse().perm == id.perm, "whose inverse is itself");

    Permutation p;
    p.perm = {2, 0, 1};
    p.iperm = {1, 2, 0};
    check(p.is_valid(), "a genuine permutation validates");
    const Permutation back = p.inverse();
    check(back.perm == p.iperm && back.iperm == p.perm, "inverse swaps the two directions");
    for (int i = 0; i < 3; ++i) {
        check(p.iperm[static_cast<std::size_t>(p.perm[static_cast<std::size_t>(i)])] == i,
              "iperm[perm[i]] == i");
    }

    // The failures that matter: a non-bijection silently loses or duplicates an
    // unknown, which downstream reads as a plausible wrong answer.
    Permutation duplicated;
    duplicated.perm = {0, 0, 2};
    duplicated.iperm = {1, 0, 2};
    check(!duplicated.is_valid(), "a repeated entry is refused");

    Permutation out_of_range;
    out_of_range.perm = {0, 3, 2};
    out_of_range.iperm = {0, 0, 2};
    check(!out_of_range.is_valid(), "an out-of-range entry is refused");

    Permutation inconsistent;
    inconsistent.perm = {2, 0, 1};
    inconsistent.iperm = {0, 1, 2};  // not the inverse
    check(!inconsistent.is_valid(), "perm and iperm that are not mutually inverse are refused");

    Permutation mismatched;
    mismatched.perm = {0, 1};
    mismatched.iperm = {0};
    check(!mismatched.is_valid(), "different lengths are refused");
}

void test_natural_is_the_identity() {
    const SparsityPattern p = pattern_from_edges(6, {{0, 3}, {1, 4}, {2, 5}, {0, 5}});
    const Permutation q = compute_ordering(p, Ordering::Natural);
    check(q.is_valid() && q.perm == Permutation::identity(6).perm,
          "the natural ordering is exactly the identity -- it is the control, so it must not "
          "quietly reorder anything");

    const SparsityPattern same = permute_pattern(p, q);
    check(same.row_ptr == p.row_ptr && same.col_index == p.col_index,
          "and permuting by it changes nothing at all");
}

void test_amd_is_refused_until_implemented() {
    const SparsityPattern p = pattern_from_edges(4, {{0, 1}, {1, 2}, {2, 3}});
    bool refused = false;
    try {
        compute_ordering(p, Ordering::ApproximateMinimumDegree);
    } catch (const std::invalid_argument&) {
        refused = true;
    }
    check(refused,
          "AMD is refused while unimplemented, rather than falling back to another ordering and "
          "reporting its fill as AMD's");
}

// The test that says whether RCM works at all.
void test_shuffled_chain() {
    const int n = 200;
    std::vector<int> label(static_cast<std::size_t>(n));
    std::iota(label.begin(), label.end(), 0);
    std::mt19937 rng(12345);  // fixed seed: this has to be reproducible
    std::shuffle(label.begin(), label.end(), rng);

    // A path graph 0-1-2-...-(n-1), relabelled by `label`. Optimal bandwidth 1.
    std::vector<std::pair<int, int>> edges;
    for (int i = 0; i + 1 < n; ++i) {
        edges.push_back({label[static_cast<std::size_t>(i)], label[static_cast<std::size_t>(i + 1)]});
    }
    const SparsityPattern shuffled = pattern_from_edges(n, edges);

    const BandwidthStats before = bandwidth_stats(shuffled);
    const Permutation p = compute_ordering(shuffled, Ordering::ReverseCuthillMcKee);
    check(p.is_valid(), "RCM returns a valid permutation");
    const SparsityPattern after = permute_pattern(shuffled, p);
    const BandwidthStats stats = bandwidth_stats(after);

    std::cout << "  shuffled chain of " << n << ": bandwidth " << before.bandwidth << " -> "
              << stats.bandwidth << ", profile " << before.profile << " -> " << stats.profile
              << "\n";

    check(before.bandwidth > 20,
          "the shuffled labelling really is bad to begin with (bandwidth " +
              std::to_string(before.bandwidth) + ")");
    check(stats.bandwidth == 1,
          "RCM recovers the optimal bandwidth of 1 on a chain (got " +
              std::to_string(stats.bandwidth) + ")");
    // n, not n - 1: the diagonal is stored, so EVERY row of a chain has a
    // furthest distance of 1, including the two ends. Worth pinning, because
    // an off-by-one here would also be consistent with a permutation that
    // left one row unordered.
    check(stats.profile == n,
          "and the optimal profile of n (got " + std::to_string(stats.profile) + ")");
}

void test_permute_preserves_structure() {
    const int n = 40;
    std::vector<std::pair<int, int>> edges;
    for (int i = 0; i + 1 < n; ++i) edges.push_back({i, i + 1});
    for (int i = 0; i + 7 < n; ++i) edges.push_back({i, i + 7});
    const SparsityPattern p = pattern_from_edges(n, edges);
    const Permutation q = compute_ordering(p, Ordering::ReverseCuthillMcKee);
    const SparsityPattern permuted = permute_pattern(p, q);

    check(permuted.nnz() == p.nnz(),
          "permuting preserves the nonzero count (" + std::to_string(permuted.nnz()) + " vs " +
              std::to_string(p.nnz()) + ")");
    check(permuted.rows == p.rows && permuted.cols == p.cols, "and the shape");

    // Every entry is still findable at its permuted position, and the diagonal
    // stays the diagonal -- a permutation is a relabelling, not a change of
    // matrix.
    bool all_found = true, diagonal_intact = true;
    for (int r = 0; r < p.rows; ++r) {
        for (int k = p.row_ptr[static_cast<std::size_t>(r)];
             k < p.row_ptr[static_cast<std::size_t>(r) + 1]; ++k) {
            const int c = p.col_index[static_cast<std::size_t>(k)];
            if (permuted.find_slot(q.iperm[static_cast<std::size_t>(r)],
                                   q.iperm[static_cast<std::size_t>(c)]) < 0) {
                all_found = false;
            }
        }
        if (permuted.find_slot(r, r) < 0) diagonal_intact = false;
    }
    check(all_found, "every entry lands at (iperm[row], iperm[col])");
    check(diagonal_intact, "and the diagonal is still the diagonal");

    // Undoing it gets the original back exactly.
    const SparsityPattern back = permute_pattern(permuted, q.inverse());
    check(back.row_ptr == p.row_ptr && back.col_index == p.col_index,
          "permuting by the inverse restores the original pattern exactly");

    // The control for the direction: using perm where iperm belongs does NOT
    // round-trip, so the round-trip above is evidence rather than a tautology.
    Permutation swapped;
    swapped.perm = q.iperm;
    swapped.iperm = q.perm;
    const SparsityPattern wrong = permute_pattern(p, swapped);
    check(!(wrong.row_ptr == permuted.row_ptr && wrong.col_index == permuted.col_index) ||
              q.perm == q.iperm,
          "swapping perm and iperm gives a different pattern -- the two directions are not "
          "interchangeable");
}

void test_isolated_and_disconnected() {
    // Row 3 has nothing but its own diagonal, which is exactly the shape of a
    // voltage port's constraint row (ASSEMBLY_PLAN §10). And {4,5} is a second
    // component. Both are normal, not edge cases to be assumed away.
    const SparsityPattern p = pattern_from_edges(6, {{0, 1}, {1, 2}, {4, 5}});
    const Permutation q = compute_ordering(p, Ordering::ReverseCuthillMcKee);
    check(q.is_valid(),
          "RCM handles an isolated unknown and several components, and still returns a "
          "permutation of everything");
    check(permute_pattern(p, q).nnz() == p.nnz(), "with the nonzero count preserved");

    // Every unknown appears exactly once -- an isolated vertex must not be
    // dropped, which would shrink the system without saying so.
    std::vector<char> seen(6, 0);
    for (int v : q.perm) seen[static_cast<std::size_t>(v)] = 1;
    check(std::count(seen.begin(), seen.end(), 1) == 6, "and nothing is left out");

    const SparsityPattern empty = pattern_from_edges(0, {});
    check(compute_ordering(empty, Ordering::ReverseCuthillMcKee).size() == 0,
          "an empty pattern is ordered without complaint");
}

void test_determinism() {
    const int n = 60;
    std::vector<std::pair<int, int>> edges;
    for (int i = 0; i + 1 < n; ++i) edges.push_back({i, i + 1});
    for (int i = 0; i + 5 < n; ++i) edges.push_back({i, i + 5});
    const SparsityPattern p = pattern_from_edges(n, edges);
    const Permutation a = compute_ordering(p, Ordering::ReverseCuthillMcKee);
    const Permutation b = compute_ordering(p, Ordering::ReverseCuthillMcKee);
    check(a.perm == b.perm,
          "RCM is deterministic -- the same pattern gives the same permutation, which is what "
          "lets a factorization built on it be compared bitwise");
}

// The real meshes. No fill to measure yet, so this reports what RCM does reduce.
void test_real_meshes() {
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

        Problem p;
        p.type = AnalysisType::DC;
        p.length_unit = LengthUnit::Millimetre;
        Body cond;
        cond.name = "B1";
        cond.volume = std::string(name)[0] == 'c' ? "wire" : "ring";
        cond.sigma = 5.8e7;
        cond.line = 1;
        Body air;
        air.name = "B2";
        air.volume = "air";
        air.sigma = 0.0;
        air.line = 2;
        p.bodies = {cond, air};
        if (std::string(name)[0] == 'c') {
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

        const BandwidthStats before = bandwidth_stats(pattern);
        const Permutation q = compute_ordering(pattern, Ordering::ReverseCuthillMcKee);
        const SparsityPattern after = permute_pattern(pattern, q);
        const BandwidthStats stats = bandwidth_stats(after);

        std::cout << "  " << name << ": " << pattern.rows << " unknowns, " << pattern.nnz()
                  << " nonzeros\n"
                  << "      bandwidth " << before.bandwidth << " -> " << stats.bandwidth
                  << ",  mean " << before.mean_bandwidth << " -> " << stats.mean_bandwidth << "\n";

        check(q.is_valid(), std::string(name) + ": RCM returns a valid permutation");
        check(after.nnz() == pattern.nnz(), std::string(name) + ": nonzero count preserved");
        check(stats.mean_bandwidth < before.mean_bandwidth,
              std::string(name) +
                  ": RCM reduces the mean bandwidth, which is the quantity that tracks fill "
                  "(from " + std::to_string(before.mean_bandwidth) + " to " +
                  std::to_string(stats.mean_bandwidth) + ")");

        // The round trip on a real pattern, where a direction error is far
        // likelier to hide than on a toy.
        const SparsityPattern restored = permute_pattern(after, q.inverse());
        check(restored.row_ptr == pattern.row_ptr && restored.col_index == pattern.col_index,
              std::string(name) + ": permuting and un-permuting restores the pattern exactly");
    }
#endif
}

}  // namespace

int main() {
    test_keywords();
    test_permutation_invariants();
    test_natural_is_the_identity();
    test_amd_is_refused_until_implemented();
    test_shuffled_chain();
    test_permute_preserves_structure();
    test_isolated_and_disconnected();
    test_determinism();
    test_real_meshes();

    std::cout << (g_checks - g_failures) << "/" << g_checks << " checks passed\n";
    return g_failures == 0 ? 0 : 1;
}
