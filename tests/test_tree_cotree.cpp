// Minimal, dependency-free test runner -- see the comment in
// tests/CMakeLists.txt for why no external framework is fetched here.
//
// Covers the boundary-first tree-cotree construction (tree_cotree.hpp;
// `docs/TREE_COTREE_BOUNDARY_FIRST.md`) and its two
// supporting pieces: the CSR NodeAdjacency builder and the UnionFind
// disjoint-set structure.
//
// These are STRUCTURAL checks -- that the tree is a genuine spanning tree
// over every node, that each Dirichlet surface is spanned by its own edges
// and entered exactly once, that thin layers keep surfaces separate. Whether
// the resulting gauge gives correct fields (B recovered exactly; no null
// space) is checked against the curl-curl operator in test_gauge_variants.cpp.

#include <array>
#include <iostream>
#include <string>
#include <vector>

#include "aphi_solver/gmsh_reader.hpp"
#include "aphi_solver/mesh.hpp"
#include "aphi_solver/tree_cotree.hpp"

using aphi_solver::boundary_edge_mask;
using aphi_solver::build_tree_cotree;
using aphi_solver::Mesh;
using aphi_solver::NodeAdjacency;
using aphi_solver::tagged_face_edge_mask;
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

Mesh make_single_tet() {
    Mesh m;
    m.nodes = {Vec3(0, 0, 0), Vec3(1, 0, 0), Vec3(0, 1, 0), Vec3(0, 0, 1)};
    m.tets = {{0, 1, 2, 3}};
    m.build_topology();
    return m;
}

Mesh make_two_tets_sharing_a_face() {
    Mesh m;
    m.nodes = {Vec3(0, 0, 0), Vec3(1, 0, 0), Vec3(0, 1, 0), Vec3(0, 0, 1), Vec3(1, 1, 1)};
    m.tets = {{0, 1, 2, 3}, {1, 2, 3, 4}};
    m.build_topology();
    return m;
}

// A structured n x n x n cube of unit cells, each split into 6 tets by the
// Kuhn triangulation -- the same construction as tools/generate_cube_mesh.py,
// built in code so these tests need no mesh file.
Mesh make_kuhn_cube(int n) {
    Mesh m;
    const int s = n + 1;
    auto id = [s](int i, int j, int k) { return i + s * j + s * s * k; };
    for (int k = 0; k <= n; ++k)
        for (int j = 0; j <= n; ++j)
            for (int i = 0; i <= n; ++i) m.nodes.emplace_back(i, j, k);
    const int kuhn[6][4][3] = {{{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {1, 1, 1}}, {{0, 0, 0}, {1, 1, 0}, {0, 1, 0}, {1, 1, 1}},
                               {{0, 0, 0}, {0, 1, 0}, {0, 1, 1}, {1, 1, 1}}, {{0, 0, 0}, {0, 1, 1}, {0, 0, 1}, {1, 1, 1}},
                               {{0, 0, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}}, {{0, 0, 0}, {1, 0, 1}, {1, 0, 0}, {1, 1, 1}}};
    for (int k = 0; k < n; ++k)
        for (int j = 0; j < n; ++j)
            for (int i = 0; i < n; ++i)
                for (const auto& t : kuhn) {
                    m.tets.push_back({id(i + t[0][0], j + t[0][1], k + t[0][2]), id(i + t[1][0], j + t[1][1], k + t[1][2]),
                                      id(i + t[2][0], j + t[2][1], k + t[2][2]), id(i + t[3][0], j + t[3][1], k + t[3][2])});
                }
    m.build_topology();
    return m;
}

double coord(const Vec3& p, int axis) { return axis == 0 ? p.x : (axis == 1 ? p.y : p.z); }

// Edges of the boundary faces lying in the plane coord[axis] == value.
std::vector<bool> plane_face_mask(const Mesh& m, int axis, double value) {
    std::vector<bool> mask(static_cast<std::size_t>(m.num_edges()), false);
    for (int f = 0; f < m.num_faces(); ++f) {
        if (!m.is_boundary_face(f)) continue;
        const auto& v = m.faces[static_cast<std::size_t>(f)];
        bool in_plane = true;
        for (int a = 0; a < 3; ++a) in_plane = in_plane && coord(m.nodes[static_cast<std::size_t>(v[a])], axis) == value;
        if (!in_plane) continue;
        mask[static_cast<std::size_t>(m.find_edge(v[0], v[1]))] = true;
        mask[static_cast<std::size_t>(m.find_edge(v[1], v[2]))] = true;
        mask[static_cast<std::size_t>(m.find_edge(v[0], v[2]))] = true;
    }
    return mask;
}

std::vector<bool> mask_or(std::vector<bool> a, const std::vector<bool>& b) {
    for (std::size_t e = 0; e < a.size(); ++e) a[e] = a[e] || b[e];
    return a;
}

int count(const std::vector<bool>& v) {
    int c = 0;
    for (bool b : v) c += b ? 1 : 0;
    return c;
}

int surface_node_count(const TreeCotreeResult& r) {
    int c = 0;
    for (int comp : r.dirichlet_component) c += (comp >= 0) ? 1 : 0;
    return c;
}

// Invariants that must hold for every result, whatever the mesh and mask.
void check_invariants(const TreeCotreeResult& r, const Mesh& m, const std::string& label) {
    const int N = m.num_nodes();
    check(static_cast<int>(r.node_group.size()) == N, label + ": node_group sized to nodes");
    check(static_cast<int>(r.is_tree_edge.size()) == m.num_edges(), label + ": is_tree_edge sized to edges");
    check(count(r.is_tree_edge) == r.tree_edge_count, label + ": tree_edge_count matches is_tree_edge");
    check(r.tree_edge_count == N - r.num_reference_groups,
          label + ": a genuine spanning tree over every node (N - roots edges)");
    check(r.surface_tree_edge_count + r.interior_tree_edge_count == r.tree_edge_count,
          label + ": surface + interior tree edges == all tree edges");
    check(r.interior_tree_edge_count == r.num_groups - r.num_reference_groups,
          label + ": one interior tree edge per non-root group");
    check(r.surface_tree_edge_count == surface_node_count(r) - r.num_dirichlet_components,
          label + ": each Dirichlet surface spanned by its own edges");

    bool groups_valid = true;
    for (int g : r.node_group) groups_valid = groups_valid && g >= 0 && g < r.num_groups;
    check(groups_valid, label + ": every node belongs to a valid group");

    // Surface tree edges are exactly the tree edges on Dirichlet surfaces;
    // every group-discovering edge is an interior (non-Dirichlet) tree edge.
    int tree_on_surface = 0;
    for (int e = 0; e < m.num_edges(); ++e) {
        if (r.is_tree_edge[static_cast<std::size_t>(e)] && r.is_dirichlet_edge[static_cast<std::size_t>(e)]) {
            ++tree_on_surface;
        }
    }
    check(tree_on_surface == r.surface_tree_edge_count, label + ": surface tree edges are the Dirichlet tree edges");
    bool discovering_ok = true;
    for (int g = 0; g < r.num_groups; ++g) {
        const int e = r.discovering_edge[static_cast<std::size_t>(g)];
        if (r.parent_group[static_cast<std::size_t>(g)] == -1) {
            discovering_ok = discovering_ok && e == -1;
        } else {
            discovering_ok = discovering_ok && e >= 0 && r.is_tree_edge[static_cast<std::size_t>(e)] &&
                             !r.is_dirichlet_edge[static_cast<std::size_t>(e)];
        }
    }
    check(discovering_ok, label + ": every non-root group is entered by exactly one interior tree edge");

    // All nodes of one surface share a group; nodes off any surface are
    // singletons.
    std::vector<int> group_size(static_cast<std::size_t>(r.num_groups), 0);
    for (int g : r.node_group) group_size[static_cast<std::size_t>(g)] += 1;
    bool grouping_ok = true;
    for (int i = 0; i < N; ++i) {
        if (r.dirichlet_component[static_cast<std::size_t>(i)] < 0) {
            grouping_ok = grouping_ok && group_size[static_cast<std::size_t>(r.node_group[static_cast<std::size_t>(i)])] == 1;
        }
    }
    check(grouping_ok, label + ": nodes off every Dirichlet surface are singleton groups (nothing collapsed)");
}

}  // namespace

int main() {
    // --- NodeAdjacency -------------------------------------------------
    {
        const Mesh m = make_single_tet();
        const NodeAdjacency adj = NodeAdjacency::build(m);
        check(static_cast<int>(adj.offset.size()) == m.num_nodes() + 1, "adjacency: offset sized num_nodes+1");
        check(static_cast<int>(adj.neighbor.size()) == 2 * m.num_edges(), "adjacency: neighbor sized 2*num_edges");
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

    // --- No Dirichlet surface: the classical construction ---------------
    {
        const Mesh m = make_single_tet();
        const TreeCotreeResult r = build_tree_cotree(m, std::vector<bool>(static_cast<std::size_t>(m.num_edges()), false));
        check_invariants(r, m, "single tet, no Dirichlet");
        check(r.num_reference_groups == 1 && r.num_groups == 4, "single tet, no Dirichlet: 1 root, 4 singleton groups");
        check(r.tree_edge_count == 3 && r.surface_tree_edge_count == 0, "single tet, no Dirichlet: 3 interior tree edges");
        check(r.node_group[0] == 0, "single tet, no Dirichlet: rooted at node 0, as before");
    }
    {
        const Mesh m = make_two_tets_sharing_a_face();
        const TreeCotreeResult r = build_tree_cotree(m, std::vector<bool>(static_cast<std::size_t>(m.num_edges()), false));
        check_invariants(r, m, "two tets, no Dirichlet");
        check(r.num_groups == 5 && r.tree_edge_count == 4, "two tets, no Dirichlet: 5 groups, 4 tree edges");
    }

    // --- One Dirichlet face on a single tet ------------------------------
    // Face (0,1,2) carries n x A = 0: those three nodes form one surface,
    // spanned by two of its three edges; node 3 is reached by one interior
    // edge.
    {
        const Mesh m = make_single_tet();
        std::vector<bool> mask(static_cast<std::size_t>(m.num_edges()), false);
        mask[static_cast<std::size_t>(m.find_edge(0, 1))] = true;
        mask[static_cast<std::size_t>(m.find_edge(1, 2))] = true;
        mask[static_cast<std::size_t>(m.find_edge(0, 2))] = true;
        const TreeCotreeResult r = build_tree_cotree(m, mask);
        check_invariants(r, m, "single tet, face (0,1,2) Dirichlet");
        check(r.num_dirichlet_components == 1, "single tet, one face: one surface");
        check(r.surface_tree_edge_count == 2 && r.interior_tree_edge_count == 1,
              "single tet, one face: 2 surface tree edges + 1 interior (N - V_b = 1)");
        check(r.node_group[0] == r.node_group[1] && r.node_group[1] == r.node_group[2] &&
                  r.node_group[3] != r.node_group[0],
              "single tet, one face: surface nodes share a group; node 3 does not");
    }

    // --- Whole outer boundary of a cube ----------------------------------
    // n = 3: 64 nodes, 8 of them interior. One surface, so the interior tree
    // must add exactly N - V_b = 8 edges.
    {
        const Mesh m = make_kuhn_cube(3);
        const TreeCotreeResult r = build_tree_cotree(m, boundary_edge_mask(m));
        check_invariants(r, m, "cube n=3, whole boundary");
        check(r.num_dirichlet_components == 1 && surface_node_count(r) == 56, "cube n=3, whole boundary: one surface, 56 nodes");
        check(r.interior_tree_edge_count == 8, "cube n=3, whole boundary: interior tree edges == N - V_b == 8");
        check(r.surface_tree_edge_count == 55, "cube n=3, whole boundary: surface tree edges == V_b - 1 == 55");
        check(r.num_reference_groups == 1, "cube n=3, whole boundary: one root");
    }

    // --- Two separate Dirichlet surfaces: ONE root, entered once ---------
    // Faces x = 0 and x = 3 do not touch. The previous construction made
    // each its own root -- one null direction too many. Here there must be
    // one root, the other surface entered by exactly one interior edge, and
    // interior tree edges == (N - V_b) + (k - 1) = 32 + 1 = 33.
    {
        const Mesh m = make_kuhn_cube(3);
        const auto mask = mask_or(plane_face_mask(m, 0, 0.0), plane_face_mask(m, 0, 3.0));
        const TreeCotreeResult r = build_tree_cotree(m, mask);
        check_invariants(r, m, "cube n=3, faces x=0 and x=3");
        check(r.num_dirichlet_components == 2, "two faces: two separate surfaces");
        check(r.num_reference_groups == 1, "two faces: ONE root, not one per surface");
        check(r.interior_tree_edge_count == 33, "two faces: interior tree edges == (N - V_b) + (k - 1) == 33");
        int entered_surfaces = 0;
        for (int g = 0; g < r.num_groups; ++g) {
            const int e = r.discovering_edge[static_cast<std::size_t>(g)];
            if (e < 0) continue;
            const auto& [i, j] = m.edges[static_cast<std::size_t>(e)];
            const int gi = r.node_group[static_cast<std::size_t>(i)];
            const int gj = r.node_group[static_cast<std::size_t>(j)];
            const int into = (gi == g) ? i : j;
            if (gi != gj && r.dirichlet_component[static_cast<std::size_t>(into)] >= 0) ++entered_surfaces;
        }
        check(entered_surfaces == 1, "two faces: the non-root surface is entered by exactly one edge");
    }

    // --- Thin layer: surfaces joined by interior edges stay separate -----
    // One cell thick: every node lies on face z = 0 or z = 1, and interior
    // edges run straight from one to the other. Unioning tagged nodes along
    // ANY edge (the previous rule) merges the two into one body -- a short
    // circuit across the layer. Unioning only along edges ON the faces keeps
    // them two surfaces.
    {
        const Mesh m = make_kuhn_cube(1);
        const auto mask = mask_or(plane_face_mask(m, 2, 0.0), plane_face_mask(m, 2, 1.0));
        bool has_crossing_edge = false;
        for (int e = 0; e < m.num_edges(); ++e) {
            const auto& [i, j] = m.edges[static_cast<std::size_t>(e)];
            if (m.nodes[static_cast<std::size_t>(i)].z != m.nodes[static_cast<std::size_t>(j)].z) has_crossing_edge = true;
        }
        check(has_crossing_edge, "thin layer: sanity, some edge joins the two faces directly");
        const TreeCotreeResult r = build_tree_cotree(m, mask);
        check_invariants(r, m, "thin layer, faces z=0 and z=1");
        check(r.num_dirichlet_components == 2, "thin layer: the two faces stay two separate surfaces");
        check(r.interior_tree_edge_count == 1, "thin layer: one interior edge, entering the second surface");
    }

    // --- Mask helpers ----------------------------------------------------
    {
        const Mesh m = make_kuhn_cube(2);
        const auto mask = boundary_edge_mask(m);
        // Closed genus-0 boundary: 12 n^2 = 48 triangles, so 72 edges.
        check(count(mask) == 72, "boundary_edge_mask, cube n=2: 72 boundary edges");
        int centre = -1;
        for (int i = 0; i < m.num_nodes(); ++i) {
            const Vec3& p = m.nodes[static_cast<std::size_t>(i)];
            if (p.x == 1 && p.y == 1 && p.z == 1) centre = i;
        }
        bool centre_edges_free = true;
        for (int e = 0; e < m.num_edges(); ++e) {
            const auto& [i, j] = m.edges[static_cast<std::size_t>(e)];
            if ((i == centre || j == centre) && mask[static_cast<std::size_t>(e)]) centre_edges_free = false;
        }
        check(centre_edges_free, "boundary_edge_mask: edges at the interior node are not boundary edges");
    }
    {
        Mesh m = make_two_tets_sharing_a_face();
        m.tagged_boundary_faces.push_back({{0, 1, 2}, 5});
        const auto mask = tagged_face_edge_mask(m, {5});
        check(count(mask) == 3 && mask[static_cast<std::size_t>(m.find_edge(0, 1))],
              "tagged_face_edge_mask: the three edges of the tagged triangle");
        check(count(tagged_face_edge_mask(m, {6})) == 0, "tagged_face_edge_mask: other tags select nothing");

        m.tagged_boundary_faces.push_back({{0, 1, 4}, 7});  // not a face of either tet
        bool threw = false;
        try {
            (void)tagged_face_edge_mask(m, {7});
        } catch (const std::invalid_argument&) {
            threw = true;
        }
        check(threw, "tagged_face_edge_mask: throws for a tagged triangle that is not a mesh face");
    }
    {
        const Mesh m = make_single_tet();
        bool threw = false;
        try {
            (void)build_tree_cotree(m, std::vector<bool>(2, false));
        } catch (const std::invalid_argument&) {
            threw = true;
        }
        check(threw, "build_tree_cotree: throws on a mask of the wrong length");
    }

    // --- Completeness on the two shipped meshes, in global terms ----------
    //
    // The gauge is complete exactly when the tree spans every node that
    // belongs to a tet:
    //
    //     tree edges == (nodes in tets) - (components among them)
    //
    // One edge short leaves a gradient in the null space and a singular
    // curl-curl block; one too many would have eliminated a real unknown.
    // Every other check in this file is about HOW the tree is built; this is
    // the one that says it did enough.
    //
    // The loop mesh is here because it is a TORUS, the case where a gauge can
    // legitimately be incomplete: on a region with first Betti number b1 > 0,
    // ker(curl) exceeds the nodal gradients by b1 and a spanning tree cannot
    // reach that part. It does not bite here, because A lives on the whole
    // BOX, which is contractible -- the cuts-for-multiply-connected problem
    // belongs to formulations that confine the potential to the air with the
    // conductor as a hole. Asserting it rather than arguing it.
    //
    // It also pins something that looks alarming and is not: loop_cut.msh
    // carries one node belonging to no tet (a point gmsh inserted and did not
    // use), so `num_reference_groups` comes out 2 rather than 1. That node has
    // no edges and no DOFs, so the tree over the real component is exactly the
    // right size, which is what this measures.
#ifdef APHI_MESH_DIR
    for (const char* name : {"cylinder_box.msh", "loop_cut.msh"}) {
        const Mesh m = aphi_solver::read_gmsh_msh(std::string(APHI_MESH_DIR) + "/" + name);
        const auto mask = boundary_edge_mask(m);
        const TreeCotreeResult r = build_tree_cotree(m, mask);
        check_invariants(r, m, std::string(name) + ", whole boundary");

        std::vector<char> used(static_cast<std::size_t>(m.num_nodes()), 0);
        for (const auto& t : m.tets) {
            for (int v : t) used[static_cast<std::size_t>(v)] = 1;
        }
        int used_nodes = 0;
        for (int i = 0; i < m.num_nodes(); ++i) used_nodes += used[static_cast<std::size_t>(i)];

        // Components of the node graph, over tet edges, counting only those
        // that contain a tet.
        UnionFind uf(m.num_nodes());
        for (const auto& t : m.tets) {
            for (int a = 0; a < 4; ++a) {
                for (int b = a + 1; b < 4; ++b) {
                    uf.unite(t[static_cast<std::size_t>(a)], t[static_cast<std::size_t>(b)]);
                }
            }
        }
        int used_components = 0;
        for (int i = 0; i < m.num_nodes(); ++i) {
            if (uf.find(i) == i && used[static_cast<std::size_t>(i)]) ++used_components;
        }

        int tree_edges = 0;
        for (int e = 0; e < m.num_edges(); ++e) {
            if (r.is_tree_edge[static_cast<std::size_t>(e)]) ++tree_edges;
        }

        check(used_components == 1, std::string(name) + ": the tets form one connected piece");
        check(tree_edges == used_nodes - used_components,
              std::string(name) + ": the tree spans every tet node -- " +
                  std::to_string(tree_edges) + " edges for " + std::to_string(used_nodes) +
                  " nodes in " + std::to_string(used_components) + " piece(s), so the gauge is "
                  "complete and the curl-curl block has no gradient left in its null space");
    }
#endif

    // --- The coaxial via -------------------------------------------------
    // n x A = 0 on the whole outer boundary. The boundary is a closed genus-0
    // surface of 1,710 triangles, so it has 2,565 edges and 857 nodes; the
    // interior tree must add 5,099 - 857 = 4,242 edges, leaving
    // 34,539 - 2,565 - 4,242 = 27,732 free A unknowns -- the count the port
    // and DOF plan predicted before any of this existed.
#ifdef APHI_MESH_DIR
    {
        const Mesh m = aphi_solver::read_gmsh_msh(std::string(APHI_MESH_DIR) + "/coax_via.msh");
        const auto mask = boundary_edge_mask(m);
        const TreeCotreeResult r = build_tree_cotree(m, mask);
        check_invariants(r, m, "coax, whole boundary");
        check(count(mask) == 2565, "coax: 2,565 boundary (Dirichlet) edges");
        check(surface_node_count(r) == 857 && r.num_dirichlet_components == 1, "coax: one surface of 857 nodes");
        check(r.interior_tree_edge_count == 4242, "coax: 4,242 interior tree edges");
        int free_edges = 0;
        for (int e = 0; e < m.num_edges(); ++e) {
            if (!r.is_tree_edge[static_cast<std::size_t>(e)] && !mask[static_cast<std::size_t>(e)]) ++free_edges;
        }
        check(free_edges == 27732, "coax: exactly 27,732 free A unknowns, as predicted");

        // The via is an ordinary region to the tree: its interior nodes are
        // singleton groups, and interior tree edges run through it.
        std::vector<bool> via_node(static_cast<std::size_t>(m.num_nodes()), false);
        for (int t = 0; t < m.num_tets(); ++t) {
            if (m.tet_tags[static_cast<std::size_t>(t)] != 1) continue;
            for (int v : m.tets[static_cast<std::size_t>(t)]) via_node[static_cast<std::size_t>(v)] = true;
        }
        int via_interior_tree_edges = 0;
        for (int e = 0; e < m.num_edges(); ++e) {
            if (!r.is_tree_edge[static_cast<std::size_t>(e)] || mask[static_cast<std::size_t>(e)]) continue;
            const auto& [i, j] = m.edges[static_cast<std::size_t>(e)];
            const bool i_in = via_node[static_cast<std::size_t>(i)] && r.dirichlet_component[static_cast<std::size_t>(i)] < 0;
            const bool j_in = via_node[static_cast<std::size_t>(j)] && r.dirichlet_component[static_cast<std::size_t>(j)] < 0;
            if (i_in && j_in) ++via_interior_tree_edges;
        }
        check(via_interior_tree_edges > 0, "coax: interior tree edges run through the via (conductor not collapsed)");
    }
#endif

    std::cout << g_checks - g_failures << "/" << g_checks << " checks passed.\n";
    return g_failures == 0 ? 0 : 1;
}
