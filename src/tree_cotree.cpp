#include "aphi_solver/tree_cotree.hpp"

#include <algorithm>
#include <stdexcept>
#include <string>
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


TreeCotreeResult build_tree_cotree(const Mesh& mesh, const std::vector<bool>& dirichlet_edge) {
    const int n = mesh.num_nodes();
    const int m = mesh.num_edges();
    if (static_cast<int>(dirichlet_edge.size()) != m) {
        throw std::invalid_argument("build_tree_cotree: dirichlet_edge must have one entry per mesh edge");
    }

    TreeCotreeResult result;
    result.node_group.assign(static_cast<std::size_t>(n), -1);
    result.is_tree_edge.assign(static_cast<std::size_t>(m), false);
    result.is_dirichlet_edge = dirichlet_edge;
    result.dirichlet_component.assign(static_cast<std::size_t>(n), -1);
    if (n == 0) {
        return result;
    }

    const NodeAdjacency adj = NodeAdjacency::build(mesh);

    // Step 1: Dirichlet surfaces are the connected components of the graph
    // made of Dirichlet edges only. Unioning along edges that lie ON a
    // surface -- never along an edge merely because both endpoints touch
    // one -- is what keeps two surfaces separated by a thin layer separate.
    UnionFind uf(n);
    std::vector<bool> on_surface(static_cast<std::size_t>(n), false);
    for (int e = 0; e < m; ++e) {
        if (!dirichlet_edge[static_cast<std::size_t>(e)]) continue;
        const auto& [i, j] = mesh.edges[static_cast<std::size_t>(e)];
        uf.unite(i, j);
        on_surface[static_cast<std::size_t>(i)] = true;
        on_surface[static_cast<std::size_t>(j)] = true;
    }
    // Component ids in order of each surface's lowest-numbered node, so
    // component 0 always holds the lowest Dirichlet node -- a deterministic
    // root choice.
    std::vector<int> component_of_root(static_cast<std::size_t>(n), -1);
    int num_components = 0;
    int surface_nodes = 0;
    for (int i = 0; i < n; ++i) {
        if (!on_surface[static_cast<std::size_t>(i)]) continue;
        ++surface_nodes;
        int& c = component_of_root[static_cast<std::size_t>(uf.find(i))];
        if (c == -1) c = num_components++;
        result.dirichlet_component[static_cast<std::size_t>(i)] = c;
    }
    result.num_dirichlet_components = num_components;
    std::vector<std::vector<int>> members(static_cast<std::size_t>(num_components));
    for (int i = 0; i < n; ++i) {
        const int c = result.dirichlet_component[static_cast<std::size_t>(i)];
        if (c >= 0) members[static_cast<std::size_t>(c)].push_back(i);
    }

    // Step 2: a spanning tree of each surface, from Dirichlet edges only.
    // These edges are zero anyway by n x A = 0; building the surface tree
    // FIRST is what guarantees the interior tree can never close a loop of
    // zero-A edges through the surface (which would pin the magnetic flux
    // through that loop to zero -- a physical constraint, not a gauge).
    int surface_tree_edges = 0;
    {
        std::vector<bool> seen(static_cast<std::size_t>(n), false);
        std::vector<int> queue;
        for (int c = 0; c < num_components; ++c) {
            const int start = members[static_cast<std::size_t>(c)].front();
            seen[static_cast<std::size_t>(start)] = true;
            queue.assign(1, start);
            for (std::size_t head = 0; head < queue.size(); ++head) {
                const int u = queue[head];
                for (int k = adj.offset[static_cast<std::size_t>(u)]; k < adj.offset[static_cast<std::size_t>(u) + 1]; ++k) {
                    const int e = adj.edge_index[static_cast<std::size_t>(k)];
                    const int v = adj.neighbor[static_cast<std::size_t>(k)];
                    if (!dirichlet_edge[static_cast<std::size_t>(e)] || seen[static_cast<std::size_t>(v)]) continue;
                    seen[static_cast<std::size_t>(v)] = true;
                    result.is_tree_edge[static_cast<std::size_t>(e)] = true;
                    ++surface_tree_edges;
                    queue.push_back(v);
                }
            }
        }
    }

    // Step 3: the interior tree, grown breadth-first from ONE root per
    // connected mesh piece. Reaching a node on a not-yet-entered surface
    // takes that whole surface as one group, through that one edge -- so no
    // surface is entered twice, and no surface is a second root (which would
    // leave a null direction). With no Dirichlet surface at all this is the
    // classical construction, rooted at node 0, visiting nodes in the same
    // order as before.
    std::vector<bool> visited(static_cast<std::size_t>(n), false);
    std::vector<int> queue;
    queue.reserve(static_cast<std::size_t>(n));
    int next_group = 0;
    int num_roots = 0;
    int interior_tree_edges = 0;

    auto take_group = [&](int node, int parent, int edge) {
        const int g = next_group++;
        result.parent_group.push_back(parent);
        result.discovering_edge.push_back(edge);
        const int c = result.dirichlet_component[static_cast<std::size_t>(node)];
        if (c >= 0) {
            for (int w : members[static_cast<std::size_t>(c)]) {
                result.node_group[static_cast<std::size_t>(w)] = g;
                visited[static_cast<std::size_t>(w)] = true;
                queue.push_back(w);
            }
        } else {
            result.node_group[static_cast<std::size_t>(node)] = g;
            visited[static_cast<std::size_t>(node)] = true;
            queue.push_back(node);
        }
    };

    std::size_t head = 0;
    auto drain = [&]() {
        while (head < queue.size()) {
            const int u = queue[head++];
            for (int k = adj.offset[static_cast<std::size_t>(u)]; k < adj.offset[static_cast<std::size_t>(u) + 1]; ++k) {
                const int v = adj.neighbor[static_cast<std::size_t>(k)];
                if (visited[static_cast<std::size_t>(v)]) continue;
                const int e = adj.edge_index[static_cast<std::size_t>(k)];
                // An edge into an unvisited node cannot lie on a surface: a
                // surface is always taken whole, so both ends of any of its
                // edges are visited together.
                if (dirichlet_edge[static_cast<std::size_t>(e)]) {
                    throw std::logic_error("build_tree_cotree: interior tree tried to use a Dirichlet edge");
                }
                result.is_tree_edge[static_cast<std::size_t>(e)] = true;
                ++interior_tree_edges;
                take_group(v, result.node_group[static_cast<std::size_t>(u)], e);
            }
        }
    };

    take_group(num_components > 0 ? members.front().front() : 0, -1, -1);
    ++num_roots;
    drain();
    for (int start = 0; start < n; ++start) {
        if (visited[static_cast<std::size_t>(start)]) continue;
        take_group(start, -1, -1);  // a further connected piece of the mesh
        ++num_roots;
        drain();
    }

    result.num_reference_groups = num_roots;
    result.num_groups = next_group;
    result.surface_tree_edge_count = surface_tree_edges;
    result.interior_tree_edge_count = interior_tree_edges;
    result.tree_edge_count = surface_tree_edges + interior_tree_edges;

    // Self-checks. Together they say: this is a genuine spanning tree over
    // every node (one per mesh piece), each surface is spanned by its own
    // edges, and the interior part adds exactly one edge per group beyond
    // the roots.
    if (surface_tree_edges != surface_nodes - num_components) {
        throw std::logic_error("build_tree_cotree: surface tree does not span each surface");
    }
    if (interior_tree_edges != result.num_groups - result.num_reference_groups) {
        throw std::logic_error("build_tree_cotree: interior tree edge / group count invariant violated");
    }
    if (result.tree_edge_count != n - result.num_reference_groups) {
        throw std::logic_error("build_tree_cotree: tree is not a spanning tree over all nodes");
    }
    return result;
}

std::vector<bool> boundary_edge_mask(const Mesh& mesh) {
    std::vector<bool> mask(static_cast<std::size_t>(mesh.num_edges()), false);
    for (int f = 0; f < mesh.num_faces(); ++f) {
        if (!mesh.is_boundary_face(f)) continue;
        const auto& v = mesh.faces[static_cast<std::size_t>(f)];
        mask[static_cast<std::size_t>(mesh.find_edge(v[0], v[1]))] = true;
        mask[static_cast<std::size_t>(mesh.find_edge(v[1], v[2]))] = true;
        mask[static_cast<std::size_t>(mesh.find_edge(v[0], v[2]))] = true;
    }
    return mask;
}

std::vector<bool> tagged_face_edge_mask(const Mesh& mesh, const std::vector<int>& tags) {
    std::vector<bool> mask(static_cast<std::size_t>(mesh.num_edges()), false);
    for (const TaggedFace& tf : mesh.tagged_boundary_faces) {
        if (std::find(tags.begin(), tags.end(), tf.tag) == tags.end()) continue;
        const auto& v = tf.nodes;
        if (mesh.find_face(v[0], v[1], v[2]) < 0) {
            throw std::invalid_argument("tagged_face_edge_mask: a triangle tagged " + std::to_string(tf.tag) +
                                        " is not a face of the tet mesh");
        }
        mask[static_cast<std::size_t>(mesh.find_edge(v[0], v[1]))] = true;
        mask[static_cast<std::size_t>(mesh.find_edge(v[1], v[2]))] = true;
        mask[static_cast<std::size_t>(mesh.find_edge(v[0], v[2]))] = true;
    }
    return mask;
}

}  // namespace aphi_solver
