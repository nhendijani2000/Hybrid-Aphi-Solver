// Tests for binding a parsed Problem to a mesh.
//
// Step 4 of `Claude outputs/input_file_plan.md` Sec. 6. Two kinds of case:
// the real cylinder fixture, where the expected numbers are known from the
// geometry; and tiny hand-built meshes, one per rejection, because a
// deliberate mistake is far easier to construct in six tetrahedra than in
// sixteen thousand.

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "aphi_solver/gmsh_reader.hpp"
#include "aphi_solver/problem_binding.hpp"

using namespace aphi_solver;

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

bool near(double a, double b, double tol = 1e-12) { return std::abs(a - b) <= tol; }

std::string expect_error(const Problem& p, const Mesh& m, const std::string& what) {
    try {
        bind_to_mesh(p, m);
    } catch (const InputError& e) {
        return e.message();
    } catch (const std::exception& e) {
        check(false, what + " -- threw the wrong exception type: " + e.what());
        return {};
    }
    check(false, what + " -- was accepted, but should have been rejected");
    return {};
}

BoundProblem expect_ok(const Problem& p, const Mesh& m, const std::string& what) {
    try {
        return bind_to_mesh(p, m);
    } catch (const std::exception& e) {
        check(false, what + " -- was rejected: " + e.what());
        return {};
    }
}

// --- a minimal hand-built mesh -------------------------------------------
//
// A unit cube split into 6 tets (Kuhn's triangulation), the same
// construction tests/test_tree_cotree.cpp uses. Tets 0..2 are tagged 1 and
// 3..5 tagged 2, so there are two "bodies" to bind, and every face of the
// cube is on the boundary.
Mesh make_cube() {
    Mesh m;
    m.nodes = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
               {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};
    m.tets = {{0, 1, 2, 6}, {0, 2, 3, 6}, {0, 3, 7, 6},
              {0, 7, 4, 6}, {0, 4, 5, 6}, {0, 5, 1, 6}};
    m.build_topology();
    m.tet_tags = {1, 1, 1, 2, 2, 2};
    m.physical_names = {{3, 1, "left"}, {3, 2, "right"}, {2, 10, "face_a"}, {2, 11, "face_b"}};
    return m;
}

// Tags a boundary face by its three nodes.
void tag_face(Mesh& m, int a, int b, int c, int tag) {
    TaggedFace f;
    f.nodes = {a, b, c};
    std::sort(f.nodes.begin(), f.nodes.end());
    f.tag = tag;
    m.tagged_boundary_faces.push_back(f);
}

Body make_body(const std::string& name, const std::string& volume, double sigma, int line = 1) {
    Body b;
    b.name = name;
    b.volume = volume;
    b.sigma = sigma;
    b.line = line;
    return b;
}

Port make_port(const std::string& name, PortType type, const std::string& surface, double amplitude,
               int line = 2) {
    Port p;
    p.name = name;
    p.type = type;
    p.surface = {surface};
    p.amplitude = amplitude;
    p.line = line;
    return p;
}

// A cube problem that binds cleanly: both bodies claimed, one current port
// and one voltage port on opposite faces. Each rejection below is this,
// with one thing changed.
Problem make_cube_problem(Mesh& m) {
    tag_face(m, 0, 1, 2, 10);  // z = 0
    tag_face(m, 0, 2, 3, 10);
    tag_face(m, 4, 5, 6, 11);  // z = 1
    tag_face(m, 4, 6, 7, 11);

    Problem p;
    p.type = AnalysisType::DC;
    p.bodies = {make_body("B1", "left", 5.8e7, 1), make_body("B2", "right", 5.8e7, 5)};
    p.ports = {make_port("P1", PortType::BoundaryCurrent, "face_a", 1.0, 9),
               make_port("P2", PortType::BoundaryVoltage, "face_b", 0.0, 13)};
    return p;
}

// ---------------------------------------------------------------------------

void test_scaling() {
    Mesh m = make_cube();
    scale_mesh_to_metres(m, LengthUnit::Millimetre);
    check(near(m.nodes[1].x, 1e-3), "mm scales a unit coordinate to 1e-3 m");
    check(near(m.nodes[6].z, 1e-3), "every coordinate is scaled");

    Mesh n = make_cube();
    scale_mesh_to_metres(n, LengthUnit::Metre);
    check(near(n.nodes[1].x, 1.0), "metres are left alone");
}

void test_cube_binds() {
    Mesh m = make_cube();
    const Problem p = make_cube_problem(m);
    const BoundProblem b = expect_ok(p, m, "the cube problem");

    check(b.bodies.size() == 2u, "two bodies");
    check(b.bodies[0].tag == 1 && b.bodies[1].tag == 2, "names resolved to tags");
    check(b.bodies[0].tets.size() == 3u && b.bodies[1].tets.size() == 3u,
          "each body owns three tets");
    check(b.body_of_tet.size() == 6u, "body_of_tet covers every tet");
    check(b.body_of_tet[0] == 0 && b.body_of_tet[5] == 1, "tets map to the right body");

    check(b.ports.size() == 2u, "two ports");
    check(b.ports[0].surface.faces.size() == 2u, "P1 resolved to two faces");
    check(b.ports[0].surface.nodes.size() == 4u, "P1's four corner nodes");
    check(b.ports[0].surface.edges.size() == 5u, "P1's five edges (two triangles share one)");

    // Positive current enters the domain, so d points inward: +z on the
    // z = 0 face, -z on the z = 1 face.
    check(near(b.ports[0].direction.z, 1.0), "the z=0 terminal's inward normal is +z");
    check(near(b.ports[1].direction.z, -1.0), "the z=1 terminal's inward normal is -z");
    check(!b.ports[0].direction_from_hint, "a boundary port's direction is always derived");

    check(b.phi_tet.size() == 6u && b.phi_tet[0], "at DC Phi lives on the conductors");
    check(b.warnings.empty(), "no warnings");
}

void test_unknown_names_are_rejected() {
    Mesh m = make_cube();
    {
        Problem p = make_cube_problem(m);
        p.bodies[0].volume = "lefft";
        const std::string msg = expect_error(p, m, "a misspelled volume name");
        check(msg.find("'left'") != std::string::npos,
              "the message lists the volume names the mesh does have");
    }
    {
        Problem p = make_cube_problem(m);
        p.ports[0].surface = {"face_zzz"};
        const std::string msg = expect_error(p, m, "a misspelled surface name");
        check(msg.find("Physical Surface") != std::string::npos,
              "the message says which dimension was searched");
    }
    {
        // A tag with no elements is different from a name that does not
        // exist, and says so.
        Problem p = make_cube_problem(m);
        p.bodies[0].volume = "77";
        const std::string msg = expect_error(p, m, "a numeric tag with no tetrahedra");
        check(msg.find("no tetrahedra") != std::string::npos,
              "an empty tag is reported as empty, not as missing");
    }
}

void test_unclaimed_and_doubly_claimed_volumes() {
    Mesh m = make_cube();
    {
        Problem p = make_cube_problem(m);
        p.bodies.pop_back();  // leave tag 2 unclaimed
        const std::string msg = expect_error(p, m, "a region no body claims");
        check(msg.find("no material") != std::string::npos ||
                  msg.find("assemble as nothing") != std::string::npos,
              "the message explains why an unclaimed region matters");
    }
    {
        Problem p = make_cube_problem(m);
        p.bodies[1].volume = "left";  // both bodies claim tag 1
        expect_error(p, m, "two bodies claiming one volume");
    }
}

void test_boundary_vs_interior_is_checked() {
    Mesh m = make_cube();
    Problem p = make_cube_problem(m);
    p.ports[0].type = PortType::InternalCurrent;
    p.ports[0].current_direction = Vec3(0, 0, 1);
    const std::string msg = expect_error(p, m, "a boundary surface declared internal");
    check(msg.find("internal") != std::string::npos,
          "the message says the declaration and the mesh disagree");
}

void test_shared_terminal_nodes_are_rejected() {
    Mesh m = make_cube();
    Problem p = make_cube_problem(m);
    // Tag a second face with tag 11 that touches tag 10's nodes.
    tag_face(m, 0, 1, 5, 11);
    const std::string msg = expect_error(p, m, "two terminals sharing nodes");
    check(msg.find("short") != std::string::npos, "the message names it a short circuit");
}

void test_dc_reference_rules() {
    {
        Mesh m = make_cube();
        Problem p = make_cube_problem(m);
        p.ports[1].type = PortType::BoundaryCurrent;  // both are now current sources
        p.ports[1].amplitude = -1.0;
        const std::string msg = expect_error(p, m, "a conductor with no potential reference");
        check(msg.find("up to a constant") != std::string::npos,
              "the message explains the singularity rather than just refusing");
    }
    {
        Mesh m = make_cube();
        Problem p = make_cube_problem(m);
        p.ports[0].type = PortType::BoundaryVoltage;  // two references on one conductor
        p.ports[0].amplitude = 1.0;
        const std::string msg = expect_error(p, m, "two references on one conductor");
        check(msg.find("over-constrained") != std::string::npos,
              "the message says the problem is over-constrained");
    }
    {
        // An insulator carrying a port can carry no DC current.
        Mesh m = make_cube();
        Problem p = make_cube_problem(m);
        p.bodies[0].sigma = 0.0;
        p.bodies[1].sigma = 0.0;
        expect_error(p, m, "ports on an all-insulator mesh");
    }
}

// The real fixture. Numbers here come from the geometry, not from a previous
// run of this code -- which is what makes them a check rather than a record.
void test_cylinder_mesh() {
#ifdef APHI_MESH_DIR
    const std::string dir = APHI_MESH_DIR;
    Mesh m;
    try {
        m = read_gmsh_msh(dir + "/cylinder_box.msh");
    } catch (const std::exception& e) {
        check(false, std::string("cylinder_box.msh unreadable: ") + e.what());
        return;
    }
    scale_mesh_to_metres(m, LengthUnit::Millimetre);

    Problem p;
    p.type = AnalysisType::DC;
    p.length_unit = LengthUnit::Millimetre;
    p.bodies = {make_body("B1", "wire", 5.8e7, 10), make_body("B2", "air", 0.0, 15)};
    p.ports = {make_port("P1", PortType::BoundaryCurrent, "wire_bottom", 1.0, 20),
               make_port("P2", PortType::BoundaryVoltage, "wire_top", 0.0, 25)};

    const BoundProblem b = expect_ok(p, m, "the cylinder problem");
    if (b.bodies.empty()) return;

    check(b.bodies[0].tets.size() == 7535u, "the wire owns 7535 tets");
    check(b.bodies[1].tets.size() == 8505u, "the air owns 8505 tets");
    check(b.bodies[0].is_conductor() && !b.bodies[1].is_conductor(),
          "the wire conducts and the air does not");

    check(b.ports[0].surface.faces.size() == 122u, "wire_bottom has 122 faces");
    check(b.ports[1].surface.faces.size() == 122u, "wire_top has 122 faces");
    check(b.ports[0].surface.nodes.size() == 74u, "wire_bottom has 74 vertices");

    // The caps are the two ends of a z-aligned prism, so the inward normals
    // are +z and -z. Getting these backwards would flip the sign of every
    // reported current.
    check(near(b.ports[0].direction.z, 1.0, 1e-9), "wire_bottom's inward normal is +z");
    check(near(b.ports[1].direction.z, -1.0, 1e-9), "wire_top's inward normal is -z");

    // At DC Phi lives only in the wire.
    int phi_tets = 0;
    for (bool v : b.phi_tet) {
        if (v) ++phi_tets;
    }
    check(phi_tets == 7535, "at DC, Phi lives on the wire's tets and no others");

    // The scaling reached the geometry: the wire's volume comes out in
    // cubic metres, not cubic millimetres -- a factor of 1e9, which is the
    // kind of error that would otherwise surface much later as a resistance
    // off by nine orders of magnitude.
    //
    // The expected value is computed from the polygon, not written as a
    // literal: a decimal literal truncated at ten digits is already 4.9e-20
    // adrift, which is larger than the tolerance this deserves. (Found by
    // this check failing on its first run.)
    const double pi = 3.14159265358979323846;
    const double a_mm = 0.2, height_mm = 1.0;
    const int N = 24;
    const double expected_m3 =
        0.5 * N * a_mm * a_mm * std::sin(2.0 * pi / N) * height_mm * 1e-9;

    double wire_volume = 0.0;
    for (int t : b.bodies[0].tets) wire_volume += std::abs(m.signed_tet_volume(t));
    check(std::abs(wire_volume / expected_m3 - 1.0) < 1e-12,
          "the wire's volume is the exact 24-gon prism, expressed in cubic metres");
    check(wire_volume < 1e-9, "the mesh really was scaled to metres, not left in mm");

    check(b.warnings.empty(), "the cylinder problem produces no warnings");
#endif
}

}  // namespace

int main() {
    test_scaling();
    test_cube_binds();
    test_unknown_names_are_rejected();
    test_unclaimed_and_doubly_claimed_volumes();
    test_boundary_vs_interior_is_checked();
    test_shared_terminal_nodes_are_rejected();
    test_dc_reference_rules();
    test_cylinder_mesh();

    std::cout << (g_checks - g_failures) << "/" << g_checks << " checks passed\n";
    return g_failures == 0 ? 0 : 1;
}
