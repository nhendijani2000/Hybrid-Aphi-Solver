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
/// near-O(1) amortized find()/unite(), used below to group PEC-tagged nodes
/// into connected "PEC bodies" without a separate BFS pass per body. This
/// is the standard, asymptotically optimal tool for discovering connected
/// components incrementally from a list of edges, and (unlike a one-off BFS
/// per body) it stays cheap if PEC tagging is later extended or updated
/// incrementally in a future phase.
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

/// Result of the tree-cotree gauge decomposition (docs/ROADMAP.md Phase 03
/// step 1; grounded in S.-C. Lee, J.-F. Lee, R. Lee, "Hierarchical Vector
/// Finite Elements for Analyzing Waveguiding Structures," IEEE Trans.
/// Microwave Theory Tech., vol. 51, no. 8, 2003, Sec. V).
struct TreeCotreeResult {
    /// node_group[i]: the DOF-numbering group index mesh node i belongs to.
    /// Every node physically on the same PEC body shares one group (held at
    /// the same potential -- Lee, Lee & Lee (2003)'s stated design
    /// principle for a single body); per docs/ROADMAP.md Phase 03 step 4,
    /// *different* PEC bodies each get their own separate group rather than
    /// being merged into one global ground. Every other node gets its own
    /// singleton group, discovered by the spanning-forest search below.
    std::vector<int> node_group;

    /// Number of distinct reference groups formed directly from PEC-tagged
    /// nodes (one per electrically distinct PEC body; at least 1 -- an
    /// arbitrary single node is used as the sole reference group if the
    /// mesh has no PEC-tagged nodes at all, so the gauge always has at
    /// least one fixed reference point).
    int num_reference_groups = 0;

    /// Total number of distinct DOF-numbering groups (num_reference_groups
    /// plus one new singleton group per node discovered during the
    /// spanning-forest search).
    int num_groups = 0;

    /// parent_group[g]: for a non-reference group g, the group id of the
    /// node it was first reached FROM during the spanning-forest search
    /// (-1 for a reference group, which has no parent). Together with
    /// discovering_edge below, this records the tree structure itself, not
    /// just which edges belong to it -- needed to walk the tree (e.g. to
    /// express a tree edge's value as a cumulative sum of cotree
    /// contributions down from the root, in the gauge-construction code
    /// this feeds -- see gauge_variants.hpp) without re-deriving it from
    /// is_tree_edge via a second graph traversal.
    std::vector<int> parent_group;

    /// discovering_edge[g]: for a non-reference group g, the global edge
    /// index of the tree edge that first reached it (-1 for a reference
    /// group). Every entry here is, by construction, an edge with
    /// is_tree_edge[discovering_edge[g]] == true.
    std::vector<int> discovering_edge;

    /// is_tree_edge[e]: true if global edge e is a spanning-tree edge (the
    /// edge a node was first reached through during the search -- per Lee,
    /// Lee & Lee (2003), its A-field DOF gets replaced by the potential
    /// difference of its two endpoint groups); false means it is a cotree
    /// edge, keeping its own independent DOF.
    std::vector<bool> is_tree_edge;

    /// Self-check invariant, always true on return (asserted in the .cpp):
    /// for a connected mesh, the number of tree edges found equals
    /// num_groups - num_reference_groups. This generalizes Lee, Lee & Lee
    /// (2003)'s single-ground invariant ("the number of tree edges is
    /// exactly the same as the number of unknowns... numbered") to this
    /// project's multiple-PEC-body requirement.
    int tree_edge_count = 0;
};

/// Builds the tree-cotree gauge decomposition for `mesh`.
///
/// `is_pec[i]` is true iff mesh node i lies on some PEC boundary (the
/// caller's own boundary-condition tagging -- this project does not yet
/// have a boundary-condition subsystem, so tests supply this directly). It
/// deliberately does *not* require the caller to have already partitioned
/// PEC nodes into separate bodies: any two PEC-tagged nodes connected by a
/// mesh edge whose both endpoints are PEC-tagged are unioned into the same
/// body automatically via UnionFind, so multiple disjoint PEC surfaces in
/// the same mesh are discovered rather than assumed.
///
/// Implementation note: Algorithm 1 (node numbering) and Algorithm 2
/// (tree/cotree edge marking) in Lee, Lee & Lee (2003) Sec. V are
/// implemented here as a *single* combined breadth-first traversal rather
/// than two separate passes over the mesh graph -- a node is assigned its
/// DOF group in the same step that the edge it was discovered through is
/// marked as a tree edge, halving the number of full graph traversals
/// compared to running the two algorithms as literally separate passes,
/// with no change in the result.
TreeCotreeResult build_tree_cotree(const Mesh& mesh, const std::vector<bool>& is_pec);

}  // namespace aphi_solver
