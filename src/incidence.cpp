#include "aphi_solver/incidence.hpp"

#include <cmath>
#include <map>

namespace aphi_solver {

void SparseMatrix::coalesce() {
    std::map<std::pair<int, int>, double> acc;
    for (const auto& e : entries) {
        acc[{e.row, e.col}] += e.value;
    }
    entries.clear();
    entries.reserve(acc.size());
    for (const auto& [rc, v] : acc) {
        entries.push_back({rc.first, rc.second, v});
    }
}

std::vector<std::vector<double>> SparseMatrix::to_dense() const {
    std::vector<std::vector<double>> d(static_cast<std::size_t>(rows),
                                        std::vector<double>(static_cast<std::size_t>(cols), 0.0));
    for (const auto& e : entries) {
        d[static_cast<std::size_t>(e.row)][static_cast<std::size_t>(e.col)] += e.value;
    }
    return d;
}

SparseMatrix multiply(const SparseMatrix& a, const SparseMatrix& b) {
    // a is rows x a.cols, b is b.rows x cols; requires a.cols == b.rows.
    SparseMatrix result;
    result.rows = a.rows;
    result.cols = b.cols;

    // Index b's entries by row for a simple row-wise gather. Fine for the
    // small, correctness-focused meshes this module targets (see the header
    // comment on SparseMatrix).
    std::map<int, std::vector<SparseMatrix::Entry>> b_by_row;
    for (const auto& be : b.entries) {
        b_by_row[be.row].push_back(be);
    }

    for (const auto& ae : a.entries) {
        auto it = b_by_row.find(ae.col);
        if (it == b_by_row.end()) continue;
        for (const auto& be : it->second) {
            result.add(ae.row, be.col, ae.value * be.value);
        }
    }
    result.coalesce();
    return result;
}

SparseMatrix transpose(const SparseMatrix& m) {
    SparseMatrix t;
    t.rows = m.cols;
    t.cols = m.rows;
    t.entries.reserve(m.entries.size());
    for (const auto& e : m.entries) {
        t.add(e.col, e.row, e.value);
    }
    t.coalesce();
    return t;
}

bool is_zero_matrix(const SparseMatrix& m, double tol) {
    for (const auto& e : m.entries) {
        if (std::fabs(e.value) > tol) return false;
    }
    return true;
}

SparseMatrix build_gradient_matrix(const Mesh& mesh) {
    SparseMatrix g;
    g.rows = mesh.num_edges();
    g.cols = mesh.num_nodes();
    for (int e = 0; e < mesh.num_edges(); ++e) {
        const auto& [i, j] = mesh.edges[static_cast<std::size_t>(e)];  // i < j by construction
        g.add(e, i, -1.0);
        g.add(e, j, +1.0);
    }
    return g;
}

SparseMatrix build_curl_matrix(const Mesh& mesh) {
    SparseMatrix c;
    c.rows = mesh.num_faces();
    c.cols = mesh.num_edges();
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
    return c;
}

}  // namespace aphi_solver
