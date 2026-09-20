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
// gmsh_reader.cpp actually remaps rather than assuming ids == 1..N.
//
// Tagging, exercising what docs/ROADMAP.md Phase 03.5 step 3 added:
//  - the two tets carry DIFFERENT physical-group tags (7 and 8), so a
//    single shared tag could not pass by accident;
//  - two type-2 triangles carry surface tags (20 and 21) -- one on the
//    shared interior face, one on an outer face -- to check they are
//    captured rather than skipped, and that their node ids are remapped and
//    sorted the same way Mesh::faces is;
//  - a type-1 line element is present to check that element types nothing
//    consumes are still skipped without breaking the parse.
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
    "5\n"
    "1 2 2 20 1 12 13 14\n"
    "2 4 2 7 1 11 12 13 14\n"
    "3 4 2 8 1 12 13 14 15\n"
    "4 2 2 21 1 11 12 13\n"
    "5 1 2 30 1 11 12\n"
    "$EndElements\n";

// A mesh whose elements carry no tags at all (num-tags == 0) -- legal Gmsh,
// and the case where tet_tags must still come back the right length rather
// than empty, since everything downstream indexes it in parallel with tets.
const char* kUntaggedMsh =
    "$MeshFormat\n"
    "2.2 0 8\n"
    "$EndMeshFormat\n"
    "$Nodes\n"
    "4\n"
    "1 0 0 0\n"
    "2 1 0 0\n"
    "3 0 1 0\n"
    "4 0 0 1\n"
    "$EndNodes\n"
    "$Elements\n"
    "1\n"
    "1 4 0 1 2 3 4\n"
    "$EndElements\n";

Mesh read_from_string(const std::string& contents, const std::string& path) {
    {
        std::ofstream out(path);
        out << contents;
    }
    Mesh m = aphi_solver::read_gmsh_msh(path);
    std::remove(path.c_str());
    return m;
}

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
    check(m.num_tets() == 2, "read_gmsh_msh: 2 tets read (triangle and line elements are not tets)");
    check(m.num_edges() == 9, "read_gmsh_msh: 9 unique edges after build_topology");
    check(m.num_faces() == 7, "read_gmsh_msh: 7 unique faces after build_topology");

    // --- Region tags (Phase 03.5, step 3) -----------------------------------
    check(m.tet_tags.size() == m.tets.size(), "tet_tags is parallel to tets");
    check(m.tet_tags.size() == 2 && m.tet_tags[0] == 7 && m.tet_tags[1] == 8,
          "tet_tags keeps each tet's own physical-group tag, in file order");

    // --- Tagged boundary faces ----------------------------------------------
    check(m.tagged_boundary_faces.size() == 2,
          "tagged_boundary_faces captures both type-2 triangles (and not the type-1 line)");

    if (m.tagged_boundary_faces.size() == 2) {
        const auto& f0 = m.tagged_boundary_faces[0];
        const auto& f1 = m.tagged_boundary_faces[1];
        check(f0.tag == 20 && f1.tag == 21, "tagged_boundary_faces keeps each triangle's surface tag");

        // Gmsh ids 12,13,14 remap to local 1,2,3; 11,12,13 to 0,1,2. Sorted
        // ascending, matching Mesh::faces' canonical orientation.
        check(f0.nodes[0] == 1 && f0.nodes[1] == 2 && f0.nodes[2] == 3,
              "tagged face nodes are remapped to 0-based local indices and sorted");
        check(f1.nodes[0] == 0 && f1.nodes[1] == 1 && f1.nodes[2] == 2,
              "second tagged face is remapped and sorted the same way");

        // Every tagged triangle in this mesh is a real face of a tet, so it
        // must resolve against the topology derived from the tets. This is
        // what makes the tags usable: Phase 04 derives is_pec from the nodes
        // of faces tagged "pec", and Phase 07's ABC needs the face indices
        // themselves.
        const int i0 = m.find_face(f0.nodes[0], f0.nodes[1], f0.nodes[2]);
        const int i1 = m.find_face(f1.nodes[0], f1.nodes[1], f1.nodes[2]);
        check(i0 >= 0 && i0 < m.num_faces(), "tagged interior face resolves via find_face");
        check(i1 >= 0 && i1 < m.num_faces(), "tagged outer face resolves via find_face");
        check(i0 != i1, "the two tagged faces resolve to distinct global faces");

        // find_face is order-independent, like find_edge.
        check(m.find_face(f0.nodes[2], f0.nodes[0], f0.nodes[1]) == i0,
              "find_face does not depend on the order its three vertices are given in");
    }

    // A vertex triple that is not a face of any tet must report -1 rather
    // than a plausible-looking index.
    check(m.find_face(0, 1, 4) == -1, "find_face returns -1 for a triple that is not a face");

    // --- Untagged mesh ------------------------------------------------------
    {
        const Mesh u = read_from_string(kUntaggedMsh, "aphi_test_gmsh_untagged_tmp.msh");
        check(u.num_tets() == 1, "untagged mesh: 1 tet read");
        check(u.tet_tags.size() == u.tets.size() && u.tet_tags[0] == -1,
              "untagged mesh: tet_tags stays parallel to tets, filled with -1");
        check(u.tagged_boundary_faces.empty(), "untagged mesh: no boundary faces captured");
    }

    // --- Hand-built mesh ----------------------------------------------------
    // build_topology must keep the invariant for meshes built in code, which
    // never go through the reader -- every other test fixture in this project
    // is one of these.
    {
        Mesh h;
        h.nodes = {aphi_solver::Vec3(0, 0, 0), aphi_solver::Vec3(1, 0, 0), aphi_solver::Vec3(0, 1, 0),
                   aphi_solver::Vec3(0, 0, 1)};
        h.tets = {{0, 1, 2, 3}};
        h.build_topology();
        check(h.tet_tags.size() == h.tets.size() && h.tet_tags[0] == -1,
              "hand-built mesh: build_topology fills tet_tags with -1");
    }

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
