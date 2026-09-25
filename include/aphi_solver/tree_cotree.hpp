#pragma once

#include <vector>

#include "aphi_solver/mesh.hpp"

namespace aphi_solver {

/// Compact CSR (compressed sparse row) adjacency for the mesh's node-edge
/// graph. Built once, O(V+E), directly from Mesh::edges -- a flat,
/// cache-friendly alternative to a std::map/adjacency-list representation.
/// Mesh already keeps a std::map for find_edge (mesh.hpp), which is the
/// right tool for the occasional point query build_curl_matrix makes (a
/// handful of lookups per face); it is the wrong tool for the
/// all-neighbors-of-every-node traversal a spanning-tree search performs
/// once per node, which is why this is built as its own flat structure
/// rather than reusing that map. See docs/ENGINEERING_STANDARDS.md.
struct NodeAdjacency {
    std::vector<int> offset;      // size num_nodes + 1
    std::vector<int> neighbor;    // size 2 * num_edges
    std::vector<int> edge_index;  // size 2 * num_edges, parallel to `neighbor`

    /// Neighbors of node `i` are neighbor[offset[i] .. offset[i+1]), with
    /// edge_index[k] giving the global edge index that connects node i to
    /// neighbor[k], for each k in that range.
    static NodeAdjacency build(const Mesh& mesh);
};

/// Disjoint-set (Union-Find) with union by rank and path compression:
/// near-O(1) amortized find()/unite(), used below to find the Dirichlet
/// surfaces: the connected components of the graph made of edges that lie on
/// an n x A = 0 face. The standard, asymptotically optimal tool for
/// discovering connected components from a list of edges.
class UnionFind {
public:
    explicit UnionFind(int n);

    /// Returns the representative (root) of the set containing `x`, with
    /// path compression (nodes visited along the way are re-parented
    /// directly to the root, flattening future lookups).
    int find(int x);

    /// Merges the sets containing `a` and `b` (no-op if already the same
    /// set), attaching the lower-rank tree under the higher-rank one to
    /// keep the resulting trees shallow.
    void unite(int a, int b);

private:
    std::vector<int> parent_;
    std::vector<int> rank_;
};

/// Result of the boundary-first tree-cotree decomposition. See
/// build_tree_cotree for the construction and
/// `docs/TREE_COTREE_BOUNDARY_FIRST.md` for why it is
/// built this way.
///
/// Two views of one spanning tree are recorded:
///
///  - **Edge view** (`is_tree_edge`): a genuine spanning tree over every mesh
///    node. Every node is its own vertex -- nothing is merged -- and the
///    tree's edges are of two kinds: *surface tree edges*, which span each
///    Dirichlet surface using only edges lying on it, and *interior tree
///    edges*, which connect everything else.
///  - **Group view** (`node_group`, `parent_group`, `discovering_edge`): the
///    same tree with each Dirichlet surface seen as one group. A surface's
///    nodes are held together by its surface tree edges, so to the rest of
///    the tree the surface behaves as one vertex; this is the view the
///    Munteanu gauge's fundamental-cycle walk needs. Every node that is not
///    on a Dirichlet surface is its own singleton group.
///
/// The group view is bookkeeping, not collapsing: no node loses its identity,
/// and conductivity plays no part in it -- a conductor's interior nodes are
/// ordinary singleton groups.
struct TreeCotreeResult {
    /// node_group[i]: the group mesh node i belongs to. Every node of one
    /// Dirichlet surface shares one group; every other node is a singleton.
    std::vector<int> node_group;

    /// Number of root groups: one per connected piece of the mesh. A root is
    /// a Dirichlet surface when the piece has one, otherwise a single node.
    /// Exactly one root per piece is what keeps the gauge complete -- one
    /// root per Dirichlet surface leaves a singular matrix.
    int num_reference_groups = 0;

    /// Total number of groups: roots, plus every entered Dirichlet surface,
    /// plus every singleton node reached by the search.
    int num_groups = 0;

    /// parent_group[g]: for a non-root group, the group it was reached from
    /// (-1 for a root). Together with discovering_edge, this is the group
    /// tree that gauge_variants walks to express tree edges through cotree
    /// edges, without re-deriving it from is_tree_edge.
    std::vector<int> parent_group;

    /// discovering_edge[g]: for a non-root group, the interior tree edge that
    /// first reached it (-1 for a root). For a Dirichlet surface this is the
    /// single edge through which the tree enters it.
    std::vector<int> discovering_edge;

    /// is_tree_edge[e]: true for every spanning-tree edge, surface and
    /// interior alike.
    std::vector<bool> is_tree_edge;

    /// is_dirichlet_edge[e]: true if edge e lies on an n x A = 0 surface --
    /// the mask passed to build_tree_cotree, kept here so a gauge reduction
    /// can eliminate these edges alongside the tree edges. Every surface tree
    /// edge is also a Dirichlet edge.
    std::vector<bool> is_dirichlet_edge;

    /// dirichlet_component[i]: which Dirichlet surface node i lies on (0 ..
    /// num_dirichlet_components - 1), or -1 if it lies on none.
    std::vector<int> dirichlet_component;
    int num_dirichlet_components = 0;

    /// All tree edges: N - num_reference_groups, i.e. N - 1 for a connected
    /// mesh.
    int tree_edge_count = 0;

    /// Tree edges lying on Dirichlet surfaces: (surface nodes) - (surfaces).
    /// Already zero by n x A = 0; recorded so the tree stays a genuine
    /// spanning tree and its invariants can be checked.
    int surface_tree_edge_count = 0;

    /// Tree edges NOT on any Dirichlet surface -- the ones the gauge
    /// actually adds as constraints. Equal to num_groups -
    /// num_reference_groups, and to (N - V_b) + (k - 1) for one connected
    /// piece with V_b Dirichlet nodes on k surfaces.
    int interior_tree_edge_count = 0;
};

/// Builds the tree-cotree gauge decomposition, boundary-first.
///
/// `dirichlet_edge[e]` is true iff edge e lies ON a surface carrying
/// n x A = 0 -- an edge of such a face, never merely an edge whose two
/// endpoints both touch one. That distinction matters: an interior edge can
/// join two surface nodes (across a corner tet, or through a thin layer),
/// and in a thin substrate a single tet edge can run from a trace to the
/// ground plane below it. Build the mask with boundary_edge_mask or
/// tagged_face_edge_mask below. An all-false mask means no Dirichlet
/// surface, and reproduces the classical construction exactly.
///
/// Construction:
///  1. Dirichlet surfaces are the connected components of the graph that
///     uses Dirichlet edges only.
///  2. Each surface gets its own spanning tree, from Dirichlet edges only.
///  3. The interior tree grows breadth-first from ONE root per connected
///     mesh piece (the surface containing the lowest-numbered Dirichlet
///     node, else node 0). When the search first reaches another surface,
///     it enters through that one edge and takes the whole surface at once,
///     so no surface is entered twice.
///
/// Why each step is needed (measured on cube_4, M = C^T C):
///  - A plain spanning tree ignoring the surfaces reaches a Dirichlet
///    surface from inside at many points. Any two such points are also
///    joined along the surface, where A is already zero, forming a closed
///    loop of zero-A edges that forces zero magnetic flux through it: a
///    physical constraint, not a gauge. The matrix stays non-singular and
///    the answer is simply wrong (26-36 % error in B). Step 2 prevents it.
///  - One root per surface leaves one null direction per extra surface (a
///    different constant potential on each root's tree): a singular matrix.
///    Step 3 prevents it.
///
/// Throws std::invalid_argument if dirichlet_edge has the wrong length.
TreeCotreeResult build_tree_cotree(const Mesh& mesh, const std::vector<bool>& dirichlet_edge);

/// Edges lying on the outer boundary of the mesh -- every edge of every face
/// with exactly one adjacent tet. The mask for "n x A = 0 on the whole
/// domain boundary", the default boundary condition.
std::vector<bool> boundary_edge_mask(const Mesh& mesh);

/// Edges lying on the tagged surface triangles whose tag is in `tags`
/// (Mesh::tagged_boundary_faces). Throws std::invalid_argument if a tagged
/// triangle is not a face of the tet mesh.
std::vector<bool> tagged_face_edge_mask(const Mesh& mesh, const std::vector<int>& tags);

}  // namespace aphi_solver
