// Minimal, dependency-free test runner -- see the comment in
// tests/CMakeLists.txt for why no external framework is fetched here.
//
// Covers docs/ROADMAP.md Phase 03 step 1 (spanning-tree search / tree-cotree
// gauge decomposition, grounded in Lee, Lee & Lee (2003) Sec. V) and its two
// supporting pieces from tree_cotree.hpp: the CSR NodeAdjacency builder and
// the UnionFind disjoint-set structure.

#include <iostream>
#include <string>

#include "aphi_solver/mesh.hpp"
#include "aphi_solver/tree_cotree.hpp"

using aphi_solver::build_tree_cotree;
using aphi_solver::Mesh;
using aphi_solver::NodeAdjacency;
using aphi_solver::TreeCotreeResult;
using aphi_solver::UnionFind;
using aphi_solver::Vec3;

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

// Same reference tet as test_incidence.cpp: nodes 0,1,2,3 at the origin and
// the three unit axes. A single tet's mesh graph is K4 (every pair of its 4
// vertices is an edge), which makes hand-checking degrees/components easy.
Mesh make_single_tet() {
    Mesh m;
    m.nodes = {Vec3(0, 0, 0), Vec3(1, 0, 0), Vec3(0, 1, 0), Vec3(0, 0, 1)};
    m.tets = {{0, 1, 2, 3}};
    m.build_topology();
    return m;
}

// Same two-tets-sharing-a-face mesh as test_incidence.cpp: 5 nodes, 9 edges.
// Nodes 0 and 4 are NOT directly connected (no edge (0,4)) -- used below to
// test that two PEC-tagged nodes with no direct PEC-PEC edge between them
// stay separate bodies rather than merging.
Mesh make_two_tets_sharing_a_face() {
    Mesh m;
    m.nodes = {Vec3(0, 0, 0), Vec3(1, 0, 0), Vec3(0, 1, 0), Vec3(0, 0, 1), Vec3(1, 1, 1)};
    m.tets = {{0, 1, 2, 3}, {1, 2, 3, 4}};
    m.build_topology();
    return m;
}

void check_invariant(const TreeCotreeResult& r, const Mesh& mesh, const std::string& label) {
    check(static_cast<int>(r.node_group.size()) == mesh.num_nodes(), label + ": node_group sized to num_nodes");
    check(static_cast<int>(r.is_tree_edge.size()) == mesh.num_edges(), label + ": is_tree_edge sized to num_edges");

    int counted_tree_edges = 0;
    for (bool b : r.is_tree_edge) {
        if (b) ++counted_tree_edges;
    }
    check(counted_tree_edges == r.tree_edge_count, label + ": tree_edge_count matches is_tree_edge's own true-count");
    check(r.tree_edge_count == r.num_groups - r.num_reference_groups,
          label + ": tree_edge_count == num_groups - num_reference_groups (Lee 2003's invariant, generalized)");

    // Every node must land in some valid group.
    for (int g : r.node_group) {
        check(g >= 0 && g < r.num_groups, label + ": every node_group entry is a valid group id");
    }
}

}  // namespace

int main() {
    // --- NodeAdjacency -------------------------------------------------
    {
        const Mesh m = make_single_tet();
        const NodeAdjacency adj = NodeAdjacency::build(m);
        check(static_cast<int>(adj.offset.size()) == m.num_nodes() + 1, "adjacency: offset sized num_nodes+1");
        check(static_cast<int>(adj.neighbor.size()) == 2 * m.num_edges(), "adjacency: neighbor sized 2*num_edges");
        // Single tet's graph is K4: every node has degree 3.
        for (int i = 0; i < m.num_nodes(); ++i) {
            const int degree = adj.offset[static_cast<std::size_t>(i) + 1] - adj.offset[static_cast<std::size_t>(i)];
            check(degree == 3, "adjacency: node " + std::to_string(i) + " has degree 3 in a single tet (K4)");
        }
    }

    // --- UnionFind -------------------------------------------------------
    {
        UnionFind uf(5);
        check(uf.find(0) != uf.find(1), "union-find: 0 and 1 start in different sets");
        uf.unite(0, 1);
        check(uf.find(0) == uf.find(1), "union-find: 0 and 1 share a set after unite(0,1)");
        uf.unite(2, 3);
        check(uf.find(0) != uf.find(2), "union-find: {0,1} and {2,3} still separate");
        uf.unite(1, 2);
        check(uf.find(0) == uf.find(3), "union-find: transitive merge via unite(1,2) joins {0,1} and {2,3}");
        check(uf.find(4) != uf.find(0), "union-find: node 4 stays in its own singleton set");
    }

    // --- build_tree_cotree: single tet, no PEC --------------------------
    {
        const Mesh m = make_single_tet();
        std::vector<bool> is_pec(static_cast<std::size_t>(m.num_nodes()), false);
        const TreeCotreeResult r = build_tree_cotree(m, is_pec);
        check_invariant(r, m, "single tet, no PEC");
        check(r.num_reference_groups == 1, "single tet, no PEC: exactly 1 fallback reference group");
        check(r.num_groups == 4, "single tet, no PEC: 4 total groups (one per node)");
        check(r.tree_edge_count == 3, "single tet, no PEC: 3 tree edges (spanning tree on 4 nodes)");
    }

    // --- build_tree_cotree: single tet, one PEC node --------------------
    {
        const Mesh m = make_single_tet();
        std::vector<bool> is_pec(static_cast<std::size_t>(m.num_nodes()), false);
        is_pec[0] = true;
        const TreeCotreeResult r = build_tree_cotree(m, is_pec);
        check_invariant(r, m, "single tet, node 0 PEC");
        check(r.num_reference_groups == 1, "single tet, node 0 PEC: exactly 1 reference group");
        check(r.node_group[0] == r.node_group[0], "single tet, node 0 PEC: sanity (node 0 has a group)");
        check(r.tree_edge_count == 3, "single tet, node 0 PEC: 3 tree edges");
    }

    // --- build_tree_cotree: single tet, two PEC nodes joined by a
    //     PEC-tagged edge -- must merge into ONE body, not two. -----------
    {
        const Mesh m = make_single_tet();
        std::vector<bool> is_pec(static_cast<std::size_t>(m.num_nodes()), false);
        is_pec[0] = true;
        is_pec[1] = true;  // edge (0,1) exists in a single tet (K4)
        const TreeCotreeResult r = build_tree_cotree(m, is_pec);
        check_invariant(r, m, "single tet, nodes 0&1 PEC (connected)");
        check(r.num_reference_groups == 1, "single tet, nodes 0&1 PEC connected: merge into 1 body");
        check(r.node_group[0] == r.node_group[1], "single tet, nodes 0&1 PEC connected: same group id");
        check(!r.is_tree_edge[static_cast<std::size_t>(m.find_edge(0, 1))],
              "single tet, nodes 0&1 PEC connected: the internal PEC-PEC edge (0,1) is not a tree edge");
        check(r.num_groups == 3, "single tet, nodes 0&1 PEC connected: 3 total groups (1 body + nodes 2,3)");
        check(r.tree_edge_count == 2, "single tet, nodes 0&1 PEC connected: 2 tree edges");
    }

    // --- build_tree_cotree: two tets sharing a face, no PEC -------------
    {
        const Mesh m = make_two_tets_sharing_a_face();
        std::vector<bool> is_pec(static_cast<std::size_t>(m.num_nodes()), false);
        const TreeCotreeResult r = build_tree_cotree(m, is_pec);
        check_invariant(r, m, "two tets, no PEC");
        check(r.num_reference_groups == 1, "two tets, no PEC: 1 fallback reference group");
        check(r.num_groups == 5, "two tets, no PEC: 5 total groups");
        check(r.tree_edge_count == 4, "two tets, no PEC: 4 tree edges (spanning tree on 5 nodes)");
    }

    // --- build_tree_cotree: two tets sharing a face, two PEC nodes with
    //     NO direct PEC-PEC edge between them -- must stay two bodies. ---
    {
        const Mesh m = make_two_tets_sharing_a_face();
        check(m.find_edge(0, 4) == -1, "two tets: sanity check, (0,4) is not an edge of this mesh");
        std::vector<bool> is_pec(static_cast<std::size_t>(m.num_nodes()), false);
        is_pec[0] = true;
        is_pec[4] = true;
        const TreeCotreeResult r = build_tree_cotree(m, is_pec);
        check_invariant(r, m, "two tets, nodes 0&4 PEC (disconnected)");
        check(r.num_reference_groups == 2, "two tets, nodes 0&4 PEC disconnected: 2 separate bodies");
        check(r.node_group[0] != r.node_group[4],
              "two tets, nodes 0&4 PEC disconnected: different group ids (not merged)");
        check(r.num_groups == 5, "two tets, nodes 0&4 PEC disconnected: 5 total groups");
        check(r.tree_edge_count == 3, "two tets, nodes 0&4 PEC disconnected: 3 tree edges");
    }

    // --- build_tree_cotree: dimension-mismatch throws -------------------
    {
        const Mesh m = make_single_tet();
        std::vector<bool> wrong_size(2, false);
        bool threw = false;
        try {
            (void)build_tree_cotree(m, wrong_size);
        } catch (const std::invalid_argument&) {
            threw = true;
        }
        check(threw, "build_tree_cotree: throws std::invalid_argument on is_pec size mismatch");
    }

    std::cout << g_checks - g_failures << "/" << g_checks << " checks passed.\n";
    return g_failures == 0 ? 0 : 1;
}
