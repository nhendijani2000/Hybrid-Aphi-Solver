// Minimal, dependency-free test runner -- see the comment in
// tests/CMakeLists.txt for why no external framework is fetched here.
//
// This is Phase 02's own exit criterion from docs/ROADMAP.md: "C and G can
// be assembled for a test mesh and CG = 0 holds numerically (the discrete
// curl.grad = 0 identity)". That is checked directly below, on two
// programmatically-built test meshes (no file I/O, so this test is robust
// to whatever working directory ctest runs it from).

#include <iostream>
#include <string>

#include "aphi_solver/incidence.hpp"
#include "aphi_solver/mesh.hpp"

using aphi_solver::Mesh;
using aphi_solver::SparseMatrix;
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

// A single reference tetrahedron: nodes 0,1,2,3 at the origin and the three
// unit axes. Smallest possible test of the incidence construction.
Mesh make_single_tet() {
    Mesh m;
    m.nodes = {Vec3(0, 0, 0), Vec3(1, 0, 0), Vec3(0, 1, 0), Vec3(0, 0, 1)};
    m.tets = {{0, 1, 2, 3}};
    m.build_topology();
    return m;
}

// Two tets sharing a face (1,2,3): tet A = (0,1,2,3), tet B = (1,2,3,4).
// Exercises global edge/face deduplication across tets, which the
// single-tet case above cannot -- the shared face's 3 edges and the face
// itself must each collapse to one global index, not be double-counted.
Mesh make_two_tets_sharing_a_face() {
    Mesh m;
    m.nodes = {
        Vec3(0, 0, 0),   // 0
        Vec3(1, 0, 0),   // 1
        Vec3(0, 1, 0),   // 2
        Vec3(0, 0, 1),   // 3
        Vec3(1, 1, 1),   // 4 -- apex of the second tet, on the other side of face (1,2,3)
    };
    m.tets = {{0, 1, 2, 3}, {1, 2, 3, 4}};
    m.build_topology();
    return m;
}

void check_cg_zero(const Mesh& mesh, const std::string& label) {
    const SparseMatrix g = aphi_solver::build_gradient_matrix(mesh);
    const SparseMatrix c = aphi_solver::build_curl_matrix(mesh);

    check(g.rows == mesh.num_edges() && g.cols == mesh.num_nodes(), label + ": G has the expected shape");
    check(c.rows == mesh.num_faces() && c.cols == mesh.num_edges(), label + ": C has the expected shape");

    const SparseMatrix cg = aphi_solver::multiply(c, g);
    check(cg.rows == mesh.num_faces() && cg.cols == mesh.num_nodes(), label + ": C*G has the expected shape");
    check(aphi_solver::is_zero_matrix(cg), label + ": C*G == 0 (discrete curl.grad = 0 identity)");
}

}  // namespace

int main() {
    {
        const Mesh m = make_single_tet();
        check(m.num_nodes() == 4, "single tet: 4 nodes");
        check(m.num_edges() == 6, "single tet: 6 edges");
        check(m.num_faces() == 4, "single tet: 4 faces");
        check_cg_zero(m, "single tet");
    }

    {
        const Mesh m = make_two_tets_sharing_a_face();
        check(m.num_nodes() == 5, "two tets: 5 nodes");
        // 2 tets x 6 edges = 12, minus the 3 edges of the shared face counted twice = 9.
        check(m.num_edges() == 9, "two tets: 9 unique edges (shared face's 3 edges deduplicated)");
        // 2 tets x 4 faces = 8, minus the 1 shared face counted twice = 7.
        check(m.num_faces() == 7, "two tets: 7 unique faces (shared face deduplicated)");
        check_cg_zero(m, "two tets sharing a face");

        // Sanity check on the shared face itself: it must resolve to exactly
        // one global face index, and each tet's tet_faces must reference it.
        const auto shared = std::array<int, 3>{1, 2, 3};
        bool found_in_a = false, found_in_b = false;
        for (int lf = 0; lf < 4; ++lf) {
            if (m.faces[static_cast<std::size_t>(m.tet_faces[0][static_cast<std::size_t>(lf)])] == shared)
                found_in_a = true;
            if (m.faces[static_cast<std::size_t>(m.tet_faces[1][static_cast<std::size_t>(lf)])] == shared)
                found_in_b = true;
        }
        check(found_in_a && found_in_b, "two tets: shared face (1,2,3) appears in both tets' tet_faces");
    }

    std::cout << g_checks - g_failures << "/" << g_checks << " checks passed.\n";
    return g_failures == 0 ? 0 : 1;
}
