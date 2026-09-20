#include "aphi_solver/tree_cotree.hpp"

#include <stdexcept>
#include <unordered_map>
#include <utility>

namespace aphi_solver {

NodeAdjacency NodeAdjacency::build(const Mesh& mesh) {
    NodeAdjacency adj;
    const int n = mesh.num_nodes();
    const int m = mesh.num_edges();

    // Bucket sizes via a counting pass, then prefix-sum into `offset` -- the
    // standard CSR construction, O(V+E) with no per-edge allocation. Not
    // parallelized: this is a single linear scan over `mesh.edges` (already
    // memory-bandwidth-bound at any mesh size this project targets), and a
    // parallel bucket-fill of the `neighbor`/`edge_index` arrays below would
    // need atomic cursor increments (or a more elaborate parallel-prefix-sum
    // scheme) to stay race-free -- not worth the complexity for a pass this
    // cheap. See docs/ENGINEERING_STANDARDS.md.
    adj.offset.assign(static_cast<std::size_t>(n) + 1, 0);
    for (const auto& e : mesh.edges) {
        adj.offset[static_cast<std::size_t>(e.first) + 1] += 1;
        adj.offset[static_cast<std::size_t>(e.second) + 1] += 1;
    }
    for (int i = 0; i < n; ++i) {
        adj.offset[static_cast<std::size_t>(i) + 1] += adj.offset[static_cast<std::size_t>(i)];
    }

    adj.neighbor.resize(static_cast<std::size_t>(2 * m));
    adj.edge_index.resize(static_cast<std::size_t>(2 * m));

    // Running write cursor per node, seeded from `offset` (copy, not a
    // reference into it -- `offset` itself must stay as the fixed
    // range boundaries callers rely on).
    std::vector<int> cursor(adj.offset.begin(), adj.offset.end() - 1);
    for (int e = 0; e < m; ++e) {
        const int i = mesh.edges[static_cast<std::size_t>(e)].first;
        const int j = mesh.edges[static_cast<std::size_t>(e)].second;

        int& ci = cursor[static_cast<std::size_t>(i)];
        adj.neighbor[static_cast<std::size_t>(ci)] = j;
        adj.edge_index[static_cast<std::size_t>(ci)] = e;
        ++ci;

        int& cj = cursor[static_cast<std::size_t>(j)];
        adj.neighbor[static_cast<std::size_t>(cj)] = i;
        adj.edge_index[static_cast<std::size_t>(cj)] = e;
        ++cj;
    }
    return adj;
}

UnionFind::UnionFind(int n) : parent_(static_cast<std::size_t>(n)), rank_(static_cast<std::size_t>(n), 0) {
    for (int i = 0; i < n; ++i) {
        parent_[static_cast<std::size_t>(i)] = i;
    }
}

int UnionFind::find(int x) {
    // Path halving: every visited node is re-parented to its grandparent on
    // the way up, which flattens the tree over repeated calls without the
    // extra recursion/second-pass a textbook full path-compression
    // implementation uses -- same O(alpha(n)) amortized bound, cheaper
    // per-call constant factor.
    while (parent_[static_cast<std::size_t>(x)] != x) {
        parent_[static_cast<std::size_t>(x)] = parent_[static_cast<std::size_t>(parent_[static_cast<std::size_t>(x)])];
        x = parent_[static_cast<std::size_t>(x)];
    }
    return x;
}

void UnionFind::unite(int a, int b) {
    int ra = find(a);
    int rb = find(b);
    if (ra == rb) return;
    if (rank_[static_cast<std::size_t>(ra)] < rank_[static_cast<std::size_t>(rb)]) {
        std::swap(ra, rb);
    }
    parent_[static_cast<std::size_t>(rb)] = ra;
    if (rank_[static_cast<std::size_t>(ra)] == rank_[static_cast<std::size_t>(rb)]) {
        rank_[static_cast<std::size_t>(ra)] += 1;
    }
}

TreeCotreeResult build_tree_cotree(const Mesh& mesh, const std::vector<bool>& is_pec) {
    const int n = mesh.num_nodes();
    if (static_cast<int>(is_pec.size()) != n) {
        throw std::invalid_argument("build_tree_cotree: is_pec must have one entry per mesh node");
    }

    TreeCotreeResult result;
    result.node_group.assign(static_cast<std::size_t>(n), -1);
    result.is_tree_edge.assign(static_cast<std::size_t>(mesh.num_edges()), false);
    // parent_group / discovering_edge are indexed by GROUP id, so they grow
    // as groups are discovered below rather than being pre-sized to `n`
    // (the final size is result.num_groups, not the mesh node count).
    result.parent_group.clear();
    result.discovering_edge.clear();
    if (n == 0) {
        return result;
    }

    const NodeAdjacency adj = NodeAdjacency::build(mesh);
    std::vector<bool> visited(static_cast<std::size_t>(n), false);
    std::vector<int> queue;
    queue.reserve(static_cast<std::size_t>(n));

    int next_group_id = 0;
    int num_reference_groups = 0;

    // Step 1 (Lee, Lee & Lee 2003 Sec. V, Algorithm 1's PEC handling):
    // union any two PEC-tagged nodes joined by an edge whose both endpoints
    // are PEC-tagged, so electrically-connected PEC surfaces are discovered
    // as a single body automatically rather than assumed by the caller.
    UnionFind uf(n);
    for (const auto& e : mesh.edges) {
        if (is_pec[static_cast<std::size_t>(e.first)] && is_pec[static_cast<std::size_t>(e.second)]) {
            uf.unite(e.first, e.second);
        }
    }

    // Step 2: one DOF group per distinct PEC body (per docs/ROADMAP.md Phase
    // 03 step 4 -- separate bodies stay separate, never merged into one
    // global ground), seeding the BFS frontier with every PEC-tagged node.
    std::unordered_map<int, int> body_group;
    bool any_pec = false;
    for (int i = 0; i < n; ++i) {
        if (!is_pec[static_cast<std::size_t>(i)]) continue;
        any_pec = true;
        const int root = uf.find(i);
        auto it = body_group.find(root);
        int gid;
        if (it == body_group.end()) {
            gid = next_group_id++;
            body_group.emplace(root, gid);
            ++num_reference_groups;
            result.parent_group.push_back(-1);
            result.discovering_edge.push_back(-1);
        } else {
            gid = it->second;
        }
        result.node_group[static_cast<std::size_t>(i)] = gid;
        visited[static_cast<std::size_t>(i)] = true;
        queue.push_back(i);
    }

    // Step 3 (Algorithm 1's fallback, steps 11-12 in the paper): if the mesh
    // has no PEC-tagged nodes at all, an arbitrary single node still needs
    // to be fixed as a reference, or the whole system stays gauge-free (an
    // all-cotree assignment with every node in its own ungrounded group,
    // which leaves the eventual A-Phi system singular).
    if (!any_pec) {
        result.node_group[0] = next_group_id++;
        visited[0] = true;
        queue.push_back(0);
        ++num_reference_groups;
        result.parent_group.push_back(-1);
        result.discovering_edge.push_back(-1);
    }

    // Step 4: Algorithm 1 (node numbering) and Algorithm 2 (tree/cotree edge
    // marking) fused into one breadth-first traversal -- a node's group id
    // is assigned in the same step that discovers the tree edge it was
    // first reached through, rather than running the two as separate graph
    // passes (see the header comment on build_tree_cotree). `queue` is
    // grown in place and read by index rather than popped from the front,
    // so this is a plain O(V+E) scan with no per-node deque overhead.
    std::size_t head = 0;
    auto drain_queue = [&]() {
        while (head < queue.size()) {
            const int u = queue[head++];
            for (int k = adj.offset[static_cast<std::size_t>(u)]; k < adj.offset[static_cast<std::size_t>(u) + 1]; ++k) {
                const int v = adj.neighbor[static_cast<std::size_t>(k)];
                if (visited[static_cast<std::size_t>(v)]) continue;
                visited[static_cast<std::size_t>(v)] = true;
                result.node_group[static_cast<std::size_t>(v)] = next_group_id++;
                result.is_tree_edge[static_cast<std::size_t>(adj.edge_index[static_cast<std::size_t>(k)])] = true;
                result.parent_group.push_back(result.node_group[static_cast<std::size_t>(u)]);
                result.discovering_edge.push_back(adj.edge_index[static_cast<std::size_t>(k)]);
                queue.push_back(v);
            }
        }
    };
    drain_queue();

    // Step 5: any mesh nodes still unvisited belong to a connected component
    // with no PEC-tagged node of its own -- not handled by the literal
    // per-node catch-up in the paper's Algorithm 1 steps 11-13, but the
    // correct generalization for an arbitrary (possibly disconnected) mesh
    // graph: each such component needs its own fallback reference root and
    // its own spanning tree, exactly like step 3 above.
    for (int start = 0; start < n; ++start) {
        if (visited[static_cast<std::size_t>(start)]) continue;
        result.node_group[static_cast<std::size_t>(start)] = next_group_id++;
        visited[static_cast<std::size_t>(start)] = true;
        ++num_reference_groups;
        result.parent_group.push_back(-1);
        result.discovering_edge.push_back(-1);
        queue.push_back(start);
        drain_queue();
    }

    result.num_reference_groups = num_reference_groups;
    result.num_groups = next_group_id;

    int tree_edge_count = 0;
    for (bool b : result.is_tree_edge) {
        if (b) ++tree_edge_count;
    }
    result.tree_edge_count = tree_edge_count;

    // Self-check: generalizes Lee, Lee & Lee (2003)'s stated invariant
    // ("the number of tree edges is exactly the same as the number of
    // unknowns... numbered in Algorithm 1") to this project's
    // multiple-reference-group construction -- one new tree edge is created
    // for every node discovered beyond the initial reference roots.
    if (result.tree_edge_count != result.num_groups - result.num_reference_groups) {
        throw std::logic_error("build_tree_cotree: tree edge count / group count invariant violated");
    }

    return result;
}

}  // namespace aphi_solver
