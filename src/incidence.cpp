#include "aphi_solver/incidence.hpp"

namespace aphi_solver {

SparseMatrix build_gradient_matrix(const Mesh& mesh) {
    SparseMatrix g(mesh.num_edges(), mesh.num_nodes());
    g.reserve(static_cast<std::size_t>(mesh.num_edges()) * 2);
    for (int e = 0; e < mesh.num_edges(); ++e) {
        const auto& [i, j] = mesh.edges[static_cast<std::size_t>(e)];  // i < j by construction
        g.add(e, i, -1.0);
        g.add(e, j, +1.0);
    }
    g.compress();
    return g;
}

SparseMatrix build_curl_matrix(const Mesh& mesh) {
    SparseMatrix c(mesh.num_faces(), mesh.num_edges());
    c.reserve(static_cast<std::size_t>(mesh.num_faces()) * 3);
    for (int f = 0; f < mesh.num_faces(); ++f) {
        const auto& face = mesh.faces[static_cast<std::size_t>(f)];  // a < b < c by construction
        const int a = face[0], b = face[1], cc = face[2];
        const int e_ab = mesh.find_edge(a, b);
        const int e_bc = mesh.find_edge(b, cc);
        const int e_ac = mesh.find_edge(a, cc);
        // These must exist: every face's 3 edges were registered as edges of
        // the same tet(s) that registered the face, in Mesh::build_topology.
        c.add(f, e_ab, +1.0);
        c.add(f, e_bc, +1.0);
        c.add(f, e_ac, -1.0);
    }
    c.compress();
    return c;
}

}  // namespace aphi_solver
