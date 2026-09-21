// Minimal, dependency-free test runner -- see tests/CMakeLists.txt.
//
// Writes a small Gmsh 2.2 ASCII mesh to a temp file, reads it back with
// read_gmsh_msh, and checks the result -- including that CG = 0 still holds
// on a mesh that came through the file-parsing path, not just the
// programmatically-built meshes in test_incidence.cpp.

#include <cmath>
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

// The SAME two-tet mesh as kSampleMsh, in Gmsh 4.1 format. Same geometry,
// same physical groups (volumes 7 and 8, surfaces 20 and 21), so the two
// readers can be checked to agree rather than each merely "not throwing".
//
// The structural differences 4.1 forces, and the reason this is a real
// parser rather than a reformat:
//  - element lines carry NO physical tag; they carry a geometric entity,
//    and $Entities maps entity -> physical group;
//  - $Nodes is blocked per entity, with node tags listed first and
//    coordinates after (this is what 4.0 does differently, hence 4.0 being
//    rejected rather than guessed at);
//  - $Elements is blocked per entity too, with the element type on the
//    block header rather than on each line.
const char* kSampleMsh41 =
    "$MeshFormat\n"
    "4.1 0 8\n"
    "$EndMeshFormat\n"
    "$PhysicalNames\n"
    "4\n"
    "3 7 \"via_conductor\"\n"
    "3 8 \"fr4 dielectric\"\n"
    "2 20 \"pec_wall\"\n"
    "2 21 \"port 1\"\n"
    "$EndPhysicalNames\n"
    "$Entities\n"
    "0 0 2 2\n"
    "1 0 0 0 1 1 1 1 20 0\n"
    "2 0 0 0 1 1 1 1 21 0\n"
    "1 0 0 0 1 1 1 1 7 0\n"
    "2 0 0 0 1 1 1 1 8 0\n"
    "$EndEntities\n"
    "$Nodes\n"
    "1 5 11 15\n"
    "3 1 0 5\n"
    "11\n"
    "12\n"
    "13\n"
    "14\n"
    "15\n"
    "0 0 0\n"
    "1 0 0\n"
    "0 1 0\n"
    "0 0 1\n"
    "1 1 1\n"
    "$EndNodes\n"
    "$Elements\n"
    "4 4 1 4\n"
    "2 1 2 1\n"
    "1 12 13 14\n"
    "2 2 2 1\n"
    "2 11 12 13\n"
    "3 1 4 1\n"
    "3 11 12 13 14\n"
    "3 2 4 1\n"
    "4 12 13 14 15\n"
    "$EndElements\n";

// Two tets sharing a face where one is flat: nodes 0,1,2,4 are coplanar
// (all at z = 0), so that tet has zero volume. A mesher can emit this from
// a sliver, and it parses perfectly well -- it only fails later, deep in
// compute_tet_geometry, unless the reader checks.
const char* kDegenerateMsh =
    "$MeshFormat\n"
    "2.2 0 8\n"
    "$EndMeshFormat\n"
    "$Nodes\n"
    "5\n"
    "1 0 0 0\n"
    "2 1 0 0\n"
    "3 0 1 0\n"
    "4 0 0 1\n"
    "5 1 1 0\n"
    "$EndNodes\n"
    "$Elements\n"
    "2\n"
    "1 4 2 7 1 1 2 3 4\n"
    "2 4 2 7 1 1 2 3 5\n"
    "$EndElements\n";

// A tet whose node list repeats a node -- also parseable, also nonsense.
const char* kRepeatedNodeMsh =
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
    "1 4 2 7 1 1 2 3 3\n"
    "$EndElements\n";

const char* kDuplicateNodeIdMsh =
    "$MeshFormat\n"
    "2.2 0 8\n"
    "$EndMeshFormat\n"
    "$Nodes\n"
    "4\n"
    "1 0 0 0\n"
    "2 1 0 0\n"
    "2 0 1 0\n"
    "4 0 0 1\n"
    "$EndNodes\n"
    "$Elements\n"
    "1\n"
    "1 4 2 7 1 1 2 2 4\n"
    "$EndElements\n";

// The same 5-node, 2-tet geometry again, but written the way Gmsh actually
// emits 4.1 rather than in the simplest form that satisfies the spec. This
// exercises the paths kSampleMsh41 does not:
//
//  - $Nodes split across FOUR entity blocks (a point, a curve, a surface, a
//    volume) instead of one. Real Gmsh always blocks per geometric entity;
//    a single-block fixture never tests the block loop at all.
//  - parametric == 1 on the curve and surface blocks, so those coordinate
//    lines carry extra u / u,v values after x y z that must be ignored.
//  - $Entities with real bounding-entity lists trailing the physical tags,
//    including negative tags (Gmsh's orientation convention), all of which
//    must be read past.
//  - An entity with ZERO physical tags (the common case when no physical
//    group was defined) -> tag -1.
//  - An entity with TWO physical tags (surface 2: 21 and 22) -> the first
//    is taken.
//  - Coordinates in scientific notation.
//
// Geometry and expected tags are identical to kSampleMsh41, so the two can
// be compared rather than each just "not throwing".
const char* kRealisticMsh41 =
    "$MeshFormat\n"
    "4.1 0 8\n"
    "$EndMeshFormat\n"
    "$PhysicalNames\n"
    "4\n"
    "3 7 \"via_conductor\"\n"
    "3 8 \"fr4 dielectric\"\n"
    "2 20 \"pec_wall\"\n"
    "2 21 \"port 1\"\n"
    "$EndPhysicalNames\n"
    "$Entities\n"
    "1 1 2 2\n"
    // point 1: tag, x y z, numPhysicalTags = 0  (no physical group)
    "1 0 0 0 0\n"
    // curve 5: tag, bbox(6), numPhys = 0, numBoundingPoints = 2, tags (one negative)
    "5 0 0 0 1 1 1 0 2 1 -1\n"
    // surface 1: bbox, numPhys = 1 (20), numBoundingCurves = 3
    "1 0 0 0 1 1 1 1 20 3 5 -5 5\n"
    // surface 2: bbox, numPhys = 2 (21, 22 -> first wins), numBoundingCurves = 1
    "2 0 0 0 1 1 1 2 21 22 1 5\n"
    // volume 1: bbox, numPhys = 1 (7), numBoundingSurfaces = 2
    "1 0 0 0 1 1 1 1 7 2 1 2\n"
    // volume 2: bbox, numPhys = 1 (8), numBoundingSurfaces = 1
    "2 0 0 0 1 1 1 1 8 1 1\n"
    "$EndEntities\n"
    "$Nodes\n"
    "4 5 11 15\n"
    // block 1: dim 0, entity 1, parametric 0, 1 node
    "0 1 0 1\n"
    "11\n"
    "0 0 0\n"
    // block 2: dim 1, entity 5, PARAMETRIC, 2 nodes -> x y z u
    "1 5 1 2\n"
    "12\n"
    "13\n"
    "1e+00 0 0 0.5\n"
    "0 1.0e+00 0 0.75\n"
    // block 3: dim 2, entity 1, PARAMETRIC, 1 node -> x y z u v
    "2 1 1 1\n"
    "14\n"
    "0 0 1 0.25 0.5\n"
    // block 4: dim 3, entity 1, parametric 0, 1 node
    "3 1 0 1\n"
    "15\n"
    "1 1 1\n"
    "$EndNodes\n"
    "$Elements\n"
    "4 4 1 4\n"
    "2 1 2 1\n"
    "1 12 13 14\n"
    "2 2 2 1\n"
    "2 11 12 13\n"
    "3 1 4 1\n"
    "3 11 12 13 14\n"
    "3 2 4 1\n"
    "4 12 13 14 15\n"
    "$EndElements\n";

const char* kMsh40 =
    "$MeshFormat\n"
    "4.0 0 8\n"
    "$EndMeshFormat\n"
    "$Nodes\n"
    "0 0\n"
    "$EndNodes\n"
    "$Elements\n"
    "0 0\n"
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

// Writes the same content with CRLF line endings, as a .msh exported from
// Gmsh on Windows would have.
Mesh read_from_string_crlf(const std::string& contents, const std::string& path) {
    {
        std::ofstream out(path, std::ios::binary);
        for (char ch : contents) {
            if (ch == '\n') out << '\r';
            out << ch;
        }
    }
    Mesh m = aphi_solver::read_gmsh_msh(path);
    std::remove(path.c_str());
    return m;
}

// Returns the GmshReadError message, or "" if the parse unexpectedly
// succeeded -- for the cases that are supposed to be rejected.
std::string expect_read_error(const std::string& contents, const std::string& path) {
    {
        std::ofstream out(path);
        out << contents;
    }
    std::string message;
    try {
        aphi_solver::read_gmsh_msh(path);
    } catch (const GmshReadError& e) {
        message = e.what();
    }
    std::remove(path.c_str());
    return message;
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

    // --- face -> tet adjacency ----------------------------------------------
    // Derived from the tets, so it is independent of what the file tagged.
    // Two tets sharing one face: that face has 2 neighbours, the other 6
    // have 1 each and are therefore the domain boundary.
    check(static_cast<int>(m.face_tets.size()) == m.num_faces(), "face_tets is parallel to faces");
    {
        int interior = 0, boundary = 0, bad = 0;
        for (int f = 0; f < m.num_faces(); ++f) {
            const int count = m.face_tets[static_cast<std::size_t>(f)].count;
            if (count == 2) ++interior;
            else if (count == 1) ++boundary;
            else ++bad;
        }
        check(interior == 1, "exactly one interior face (the shared one)");
        check(boundary == 6, "the other six faces are boundary faces");
        check(bad == 0, "no face has 0 or >2 adjacent tets");

        // The shared face must name both tets; a boundary face only one.
        const int shared = m.find_face(1, 2, 3);
        check(shared >= 0 && !m.is_boundary_face(shared), "the shared face is not a boundary face");
        if (shared >= 0) {
            const auto& ft = m.face_tets[static_cast<std::size_t>(shared)];
            check((ft.tets[0] == 0 && ft.tets[1] == 1) || (ft.tets[0] == 1 && ft.tets[1] == 0),
                  "the shared face names both tets");
        }
        const int outer = m.find_face(0, 1, 2);
        check(outer >= 0 && m.is_boundary_face(outer), "an outer face is reported as boundary");
    }

    // --- Gmsh 4.1 -----------------------------------------------------------
    // The same mesh in the format current Gmsh actually writes. Checked
    // against the 2.2 result rather than only for self-consistency.
    {
        const Mesh m41 = read_from_string(kSampleMsh41, "aphi_test_gmsh_v41_tmp.msh");
        check(m41.num_nodes() == 5 && m41.num_tets() == 2, "4.1: same node and tet counts as 2.2");
        check(m41.num_edges() == m.num_edges() && m41.num_faces() == m.num_faces(),
              "4.1: same derived topology as the 2.2 file");
        check(m41.tet_tags.size() == 2 && m41.tet_tags[0] == 7 && m41.tet_tags[1] == 8,
              "4.1: tet physical tags resolved through $Entities");
        check(m41.tagged_boundary_faces.size() == 2, "4.1: both tagged surface triangles captured");
        if (m41.tagged_boundary_faces.size() == 2) {
            check(m41.tagged_boundary_faces[0].tag == 20 && m41.tagged_boundary_faces[1].tag == 21,
                  "4.1: surface physical tags resolved through $Entities");
        }
        check(m41.nodes[0].x == 0.0 && m41.nodes[4].x == 1.0,
              "4.1: blocked $Nodes (tags then coordinates) read in the right order");

        // $PhysicalNames, including one with a space in it.
        check(m41.physical_names.size() == 4, "4.1: all four $PhysicalNames entries read");
        check(m41.physical_name(3, 7) == "via_conductor", "physical_name resolves a volume group");
        check(m41.physical_name(3, 8) == "fr4 dielectric", "physical_name keeps spaces inside the quotes");
        check(m41.physical_name(2, 20) == "pec_wall", "physical_name resolves a surface group");
        check(m41.physical_name(2, 21) == "port 1", "physical_name distinguishes surface tag 21");
        check(m41.physical_name(3, 20).empty(),
              "physical_name is keyed by (dimension, tag) -- surface tag 20 is not volume tag 20");
        check(m41.physical_name(3, 999).empty(), "physical_name returns empty for an unknown group");

        const auto g41 = aphi_solver::build_gradient_matrix(m41);
        const auto c41 = aphi_solver::build_curl_matrix(m41);
        check(c41.multiply(g41).is_zero(), "4.1: CG = 0 holds on the 4.1-parsed mesh");

        // --- 4.1 as Gmsh actually writes it ---------------------------------
        // Same geometry and tags, but multi-block $Nodes, parametric
        // coordinates, real bounding-entity lists, an entity with no
        // physical tags and one with two. Compared against the simple 4.1
        // fixture rather than only checked for self-consistency.
        const Mesh mr = read_from_string(kRealisticMsh41, "aphi_test_gmsh_v41r_tmp.msh");
        check(mr.num_nodes() == 5 && mr.num_tets() == 2,
              "4.1 realistic: multi-block $Nodes yields the same node and tet counts");
        check(mr.num_edges() == m41.num_edges() && mr.num_faces() == m41.num_faces(),
              "4.1 realistic: same derived topology as the single-block fixture");
        check(mr.tet_tags.size() == 2 && mr.tet_tags[0] == 7 && mr.tet_tags[1] == 8,
              "4.1 realistic: volume physical tags resolved past bounding-entity lists");
        check(mr.tagged_boundary_faces.size() == 2 && mr.tagged_boundary_faces[0].tag == 20 &&
                  mr.tagged_boundary_faces[1].tag == 21,
              "4.1 realistic: surface with TWO physical tags resolves to the first (21, not 22)");

        // The parametric blocks are the real trap: those coordinate lines
        // carry u / u,v after x y z, so a reader that consumed the whole
        // line would misplace every node on a curve or surface.
        bool coords_match = true;
        for (int i = 0; i < mr.num_nodes(); ++i) {
            const aphi_solver::Vec3 a = mr.nodes[static_cast<std::size_t>(i)];
            const aphi_solver::Vec3 b = m41.nodes[static_cast<std::size_t>(i)];
            if (std::abs(a.x - b.x) > 1e-12 || std::abs(a.y - b.y) > 1e-12 ||
                std::abs(a.z - b.z) > 1e-12) {
                coords_match = false;
            }
        }
        check(coords_match,
              "4.1 realistic: parametric blocks and scientific notation give the same coordinates");

        const auto gr = aphi_solver::build_gradient_matrix(mr);
        const auto cr = aphi_solver::build_curl_matrix(mr);
        check(cr.multiply(gr).is_zero(), "4.1 realistic: CG = 0");

        // --- CRLF -----------------------------------------------------------
        // A .msh exported from Gmsh on Windows has CRLF endings.
        const Mesh mc = read_from_string_crlf(kRealisticMsh41, "aphi_test_gmsh_crlf_tmp.msh");
        check(mc.num_nodes() == 5 && mc.num_tets() == 2 && mc.num_faces() == m41.num_faces(),
              "CRLF line endings parse identically");
        check(mc.tet_tags == mr.tet_tags, "CRLF: tags unaffected by line endings");
        check(mc.physical_name(2, 21) == "port 1",
              "CRLF: a quoted $PhysicalNames value keeps no trailing carriage return");
    }

    // A 2.2 file with no $PhysicalNames has none -- and that is not an error.
    check(m.physical_names.empty(), "2.2 sample without $PhysicalNames yields no names");

    // --- Validation ---------------------------------------------------------
    {
        const std::string degenerate = expect_read_error(kDegenerateMsh, "aphi_test_gmsh_degen_tmp.msh");
        check(degenerate.find("degenerate") != std::string::npos,
              "a flat (zero-volume) tet is rejected at read time, not at first use");

        const std::string repeated = expect_read_error(kRepeatedNodeMsh, "aphi_test_gmsh_repeat_tmp.msh");
        check(repeated.find("same node twice") != std::string::npos,
              "a tet that repeats a node is rejected");

        const std::string duplicate = expect_read_error(kDuplicateNodeIdMsh, "aphi_test_gmsh_dupid_tmp.msh");
        check(duplicate.find("duplicate node id") != std::string::npos,
              "a repeated node id is rejected rather than silently overwriting");

        const std::string v40 = expect_read_error(kMsh40, "aphi_test_gmsh_v40_tmp.msh");
        check(v40.find("4.0") != std::string::npos || v40.find("unsupported") != std::string::npos,
              "MSH 4.0 is rejected with an explanatory message rather than misread");
    }

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
