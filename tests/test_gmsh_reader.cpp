// Minimal, dependency-free test runner -- see tests/CMakeLists.txt.
//
// Writes a small Gmsh 2.2 ASCII mesh to a temp file, reads it back with
// read_gmsh_msh, and checks the result -- including that CG = 0 still holds
// on a mesh that came through the file-parsing path, not just the
// programmatically-built meshes in test_incidence.cpp.

#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>

#include "aphi_solver/gmsh_reader.hpp"
#include "aphi_solver/incidence.hpp"

using aphi_solver::GmshReadError;
using aphi_solver::Mesh;

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

// Two tets sharing a face (1,2,3), same geometry as
// test_incidence.cpp's make_two_tets_sharing_a_face, but written out as a
// Gmsh 2.2 ASCII file -- deliberately using non-contiguous, non-1-based-
// from-the-start Gmsh node ids (11..15) to check the id remapping in
// gmsh_reader.cpp actually remaps rather than assuming ids == 1..N. A
// triangle element (type 2) is included on the shared face to check that
// non-tet elements are correctly skipped rather than breaking the parse.
const char* kSampleMsh =
    "$MeshFormat\n"
    "2.2 0 8\n"
    "$EndMeshFormat\n"
    "$Nodes\n"
    "5\n"
    "11 0 0 0\n"
    "12 1 0 0\n"
    "13 0 1 0\n"
    "14 0 0 1\n"
    "15 1 1 1\n"
    "$EndNodes\n"
    "$Elements\n"
    "3\n"
    "1 2 2 0 1 12 13 14\n"
    "2 4 2 0 1 11 12 13 14\n"
    "3 4 2 0 1 12 13 14 15\n"
    "$EndElements\n";

}  // namespace

int main() {
    const std::string path = "aphi_test_gmsh_reader_tmp.msh";
    {
        std::ofstream out(path);
        out << kSampleMsh;
    }

    Mesh m;
    bool threw = false;
    try {
        m = aphi_solver::read_gmsh_msh(path);
    } catch (const GmshReadError& e) {
        threw = true;
        std::cerr << "unexpected GmshReadError: " << e.what() << "\n";
    }
    std::remove(path.c_str());

    check(!threw, "read_gmsh_msh parses the sample file without throwing");
    check(m.num_nodes() == 5, "read_gmsh_msh: 5 nodes read (non-contiguous Gmsh ids remapped)");
    check(m.num_tets() == 2, "read_gmsh_msh: 2 tets read (the type-2 triangle element was skipped)");
    check(m.num_edges() == 9, "read_gmsh_msh: 9 unique edges after build_topology");
    check(m.num_faces() == 7, "read_gmsh_msh: 7 unique faces after build_topology");

    const auto g = aphi_solver::build_gradient_matrix(m);
    const auto c = aphi_solver::build_curl_matrix(m);
    const auto cg = c.multiply(g);
    check(cg.is_zero(), "read_gmsh_msh: CG = 0 holds on a mesh read from a file");

    // A missing file must raise GmshReadError, not crash or return silently.
    bool missing_file_threw = false;
    try {
        aphi_solver::read_gmsh_msh("aphi_this_file_does_not_exist.msh");
    } catch (const GmshReadError&) {
        missing_file_threw = true;
    }
    check(missing_file_threw, "read_gmsh_msh throws GmshReadError on a missing file");

    std::cout << g_checks - g_failures << "/" << g_checks << " checks passed.\n";
    return g_failures == 0 ? 0 : 1;
}
