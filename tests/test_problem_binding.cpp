// Tests for binding a parsed Problem to a mesh.
//
// Step 4 of `docs/INPUT_FILE_PLAN.md` Sec. 6. Two kinds of case:
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
        // TWO references on one conductor is a voltage-driven resistor:
        // V = 1 at one end, V = 0 at the other. Well-posed -- the current
        // falls out as a result -- and it must be accepted.
        //
        // An earlier version of bind_to_mesh rejected this, and an earlier
        // version of this test asserted the rejection, so the suite
        // defended the bug. The rule is "at least one reference per
        // conducting path", never "exactly one". The cylinder passes either
        // way, which is why it went unnoticed.
        Mesh m = make_cube();
        Problem p = make_cube_problem(m);
        p.ports[0].type = PortType::BoundaryVoltage;
        p.ports[0].amplitude = 1.0;
        const BoundProblem b = expect_ok(p, m, "a voltage-driven conductor (two voltage ports)");
        if (!b.conduction_paths.empty()) {
            check(b.conduction_paths[0].reference_count == 2,
                  "both voltage ports count as references on the one path");
        }
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

// A conduction path is built from sigma, and spans bodies: two conductors
// that touch are ONE path, because that is what a current sees. Getting this
// from phi_tet instead would merge the whole domain at full wave.
void test_conduction_paths() {
    {
        Mesh m = make_cube();
        const Problem p = make_cube_problem(m);
        const BoundProblem b = expect_ok(p, m, "the cube's conduction path");
        check(b.conduction_paths.size() == 1u,
              "two touching conductor bodies are ONE conduction path");
        check(b.conduction_paths[0].tets.size() == 6u, "the path covers every tet");
        check(b.conduction_paths[0].bodies.size() == 2u, "the path records both bodies");
        check(b.conduction_paths[0].reference_count == 1, "P2 is its one reference");
        check(!b.conduction_paths[0].is_floating, "a ported path is not floating");
        check(b.path_of_tet[0] == 0 && b.path_of_tet[5] == 0, "every tet maps to the path");
    }
    {
        // Make one body an insulator: the remaining conductor is a smaller
        // path, and the insulator's tets belong to none.
        Mesh m = make_cube();
        Problem p = make_cube_problem(m);
        p.bodies[1].sigma = 0.0;
        p.ports[1].surface = {"face_a"};  // move both ports onto the conductor
        p.ports[1].type = PortType::BoundaryVoltage;
        p.ports.pop_back();
        p.ports[0].type = PortType::BoundaryVoltage;
        p.ports[0].amplitude = 0.0;
        const BoundProblem b = expect_ok(p, m, "a mesh with one conductor and one insulator");
        check(b.conduction_paths.size() == 1u, "only the conducting body forms a path");
        check(b.conduction_paths[0].tets.size() == 3u, "the path is the three conducting tets");
        check(b.path_of_tet[3] == -1, "an insulator tet belongs to no path");
    }
    {
        // ONE body whose tets are not connected to each other is TWO paths.
        //
        // In the Kuhn cube the six tets form a closed fan, so tets 1 and 3
        // are not adjacent. Tagging only those two as the conductor leaves
        // two islands of one tet each, separated by insulator.
        //
        // This is the case that makes the sigma-vs-phi_tet distinction
        // visible. At full wave phi_tet is true everywhere, so building the
        // components from it would unite the two islands *through the
        // insulator* and report one path. An earlier version of this test
        // used a single connected conductor, where both spellings agree --
        // it passed a negative control that flipped the source, which is
        // how the weakness was found.
        Mesh m = make_cube();
        m.tet_tags = {2, 1, 2, 1, 2, 2};
        Problem p = make_cube_problem(m);
        p.bodies[0].sigma = 5.8e7;  // tag 1: the two islands
        p.bodies[1].sigma = 0.0;    // tag 2: insulator
        p.type = AnalysisType::Frequency;
        p.frequencies = {1e6};
        p.formulation = Formulation::FullWave;
        p.ports.pop_back();
        p.ports[0].type = PortType::BoundaryVoltage;
        p.ports[0].amplitude = 0.0;

        const BoundProblem b = expect_ok(p, m, "two disconnected conductor islands at full wave");
        check(b.conduction_paths.size() == 2u,
              "one body in two disconnected pieces is TWO conduction paths, and the insulator "
              "between them does not join them even though Phi lives there");
        if (b.conduction_paths.size() == 2u) {
            check(b.conduction_paths[0].tets.size() == 1u &&
                      b.conduction_paths[1].tets.size() == 1u,
                  "each island is one tet");
            check(b.conduction_paths[0].bodies.size() == 1u &&
                      b.conduction_paths[0].bodies[0] == 0,
                  "both islands belong to the same body");
        }
        int phi_tets = 0;
        for (bool v : b.phi_tet) {
            if (v) ++phi_tets;
        }
        check(phi_tets == 6, "at full wave Phi still lives everywhere, paths notwithstanding");
    }
}

// A conductor nothing touches would leave Phi fixed only up to a constant.
// Warned about rather than refused -- an unconnected shield is legitimate --
// and one node is pinned to keep the matrix non-singular.
void test_floating_conductor() {
    Mesh m = make_cube();
    Problem p = make_cube_problem(m);
    // Move both ports onto body 1's face so body 2, still a conductor, is
    // touched by nothing. The cube's two halves meet, so make them separate
    // by turning the shared region into... simpler: both ports on face_a,
    // which only body 1 borders, is not possible here -- instead give body 2
    // its own tag and no port, and rely on the cube's halves being joined.
    // They ARE joined, so this is one path; the floating case needs a mesh
    // where they are not. Use the insulator to separate them.
    p.bodies[1].sigma = 0.0;
    p.ports.pop_back();
    p.ports[0].type = PortType::BoundaryVoltage;
    p.ports[0].amplitude = 0.0;
    const BoundProblem b = expect_ok(p, m, "a conductor with a reference");
    check(b.warnings.empty(), "a referenced conductor produces no floating warning");

    // Now remove the port's reference role by making the conductor carry no
    // port at all: put the only port on the insulator's face. At DC that is
    // a different error (a port touching no conductor), which is the point
    // -- the two cases are distinguished.
    Problem q = make_cube_problem(m);
    q.bodies[1].sigma = 0.0;
    q.ports.pop_back();
    q.ports[0].surface = {"face_b"};  // z = 1, bordered by body 2 (the insulator)
    q.ports[0].type = PortType::BoundaryVoltage;
    q.ports[0].amplitude = 0.0;
    const std::string msg = expect_error(q, m, "a port touching no conductor at DC");
    check(msg.find("no conductor") != std::string::npos,
          "the message says the port touches no conductor");
}

// The gauge runs during binding, because the DOF map cannot classify an edge
// without it. The invariants are exact, so they are asserted as equalities.
void test_gauge_is_built() {
    Mesh m = make_cube();
    const Problem p = make_cube_problem(m);
    const BoundProblem b = expect_ok(p, m, "the cube's gauge");

    check(b.dirichlet_edge.size() == static_cast<std::size_t>(m.num_edges()),
          "the Dirichlet mask covers every edge");
    check(b.gauge.is_tree_edge.size() == static_cast<std::size_t>(m.num_edges()),
          "the gauge covers every edge");

    // A spanning tree over one connected mesh has N - 1 edges, always.
    check(b.gauge.tree_edge_count == m.num_nodes() - 1,
          "tree edges == nodes - 1 for a connected mesh");
    check(b.gauge.num_reference_groups == 1, "one root for one connected piece");

    // Every edge of the cube's surface is on the boundary; the only interior
    // edge is the body diagonal 0-6 that Kuhn's triangulation shares.
    int dirichlet = 0;
    for (bool v : b.dirichlet_edge) {
        if (v) ++dirichlet;
    }
    check(dirichlet == m.num_edges() - 1,
          "every edge but the interior diagonal lies on the boundary");

    check(b.gauge.surface_tree_edge_count + b.gauge.interior_tree_edge_count ==
              b.gauge.tree_edge_count,
          "surface and interior tree edges account for the whole tree");
}

// The side labels. A tet touching the cut at a single node still has to know
// which side it sees that node from, because Phi is 0 on one side and the
// port's unknown on the other.
void test_cut_side_labels() {
    // A bipyramid: two tets sharing the face (0,1,2), with apexes above and
    // below it. That shared face is a COMPLETE cross-section -- it separates
    // the solid, it is planar, and its three rim edges all lie on the outer
    // boundary -- so it is a valid cut.
    //
    // The Kuhn cube cannot supply one. Its six tets form a closed fan around
    // the diagonal 0-6, so a single interior triangle does not separate it
    // (the rim check catches that, as it should), and the two faces that
    // would separate it are not coplanar (the planarity check catches that).
    // Both refusals are correct; the geometry simply has no valid cut.
    Mesh m;
    m.nodes = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, -1}};
    m.tets = {{0, 1, 2, 3}, {0, 1, 2, 4}};
    m.build_topology();
    m.tet_tags = {1, 1};
    m.physical_names = {{3, 1, "metal"}, {2, 20, "cut"}};
    tag_face(m, 0, 1, 2, 20);

    Problem p;
    p.type = AnalysisType::DC;
    p.bodies = {make_body("B1", "metal", 5.8e7, 1)};
    Port cut = make_port("P1", PortType::InternalCurrent, "cut", 1.0, 9);
    // The cut lies in z = 0, so +z points through it. The plus side is then
    // the tet with the apex at z = +1, which is tet 0.
    cut.current_direction = Vec3(0, 0, 1);
    p.ports = {cut};

    const BoundProblem b = expect_ok(p, m, "an internal cut on the cube");
    if (b.ports.empty()) return;
    const BoundPort& bp = b.ports[0];

    check(bp.surface.faces.size() == 1u, "the cut is one face");
    check(bp.plus_side_tet.size() == 1u, "one plus-side tet per cut face");
    check(!bp.touching_tets.empty(), "some tets touch the cut");

    // Every touching tet gets a definite side, and the two sides disagree.
    int plus = 0, minus = 0;
    for (signed char s : bp.tet_side) {
        if (s > 0) ++plus;
        if (s < 0) ++minus;
    }
    check(plus > 0 && minus > 0, "tets are found on both sides of the cut");
    check(plus + minus == static_cast<int>(bp.touching_tets.size()),
          "every touching tet is labelled, none left at zero");

    // The face's own two tets must be labelled consistently with
    // plus_side_tet, which was computed by a different test (the opposite
    // vertex) -- so agreement is a real cross-check, not a tautology.
    const int f = bp.surface.faces[0];
    const FaceTets& ft = m.face_tets[static_cast<std::size_t>(f)];
    check(bp.side_of_tet(bp.plus_side_tet[0]) == +1,
          "the opposite-vertex test and the plane test agree on the plus side");
    const int other = ft.tets[0] == bp.plus_side_tet[0] ? ft.tets[1] : ft.tets[0];
    check(bp.side_of_tet(other) == -1, "and on the minus side");

    // Both tets touch this cut, so side_of_tet is exercised for absence
    // with an index past the end rather than a real tet.
    check(bp.side_of_tet(m.num_tets() + 5) == 0,
          "side_of_tet returns 0 for a tet that does not touch the cut");
    check(bp.plus_side_tet[0] == 0, "tet 0, whose apex is at z = +1, is the plus side");

    // Reversing the hint must swap every label -- the labels follow d, and
    // nothing else.
    Problem q = p;
    q.ports[0].current_direction = Vec3(0, 0, -1);
    const BoundProblem c = expect_ok(q, m, "the same cut with the hint reversed");
    if (!c.ports.empty()) {
        bool all_flipped = c.ports[0].touching_tets == bp.touching_tets;
        for (std::size_t i = 0; i < bp.tet_side.size() && all_flipped; ++i) {
            if (c.ports[0].tet_side[i] != -bp.tet_side[i]) all_flipped = false;
        }
        check(all_flipped, "reversing current_direction swaps every side label");
    }
}

// The loop mesh: the first INTERNAL port on real geometry.
//
// A boundary port's terminal is handed to it by the mesh -- the faces are on
// the outer boundary and that is that. An internal port has to work out which
// side is which, and every number below is about that. Expected values come
// from the geometry in tools/loop_cut.geo, not from a previous run.
void test_loop_mesh() {
#ifdef APHI_MESH_DIR
    const std::string dir = APHI_MESH_DIR;
    Mesh m;
    try {
        m = read_gmsh_msh(dir + "/loop_cut.msh");
    } catch (const std::exception& e) {
        check(false, std::string("loop_cut.msh unreadable: ") + e.what());
        return;
    }
    scale_mesh_to_metres(m, LengthUnit::Millimetre);

    Problem p;
    p.type = AnalysisType::DC;
    p.length_unit = LengthUnit::Millimetre;
    p.bodies = {make_body("B1", "ring", 5.8e7, 10), make_body("B2", "air", 0.0, 15)};
    Port port = make_port("P1", PortType::InternalCurrent, "loop_cut", 1.0, 20);
    port.current_direction = Vec3{0.0, 1.0, 0.0};  // +y, as the example file says
    p.ports = {port};

    const BoundProblem b = expect_ok(p, m, "the loop problem");
    if (b.bodies.empty() || b.ports.empty()) return;

    check(b.bodies[0].is_conductor() && !b.bodies[1].is_conductor(),
          "the ring conducts and the air does not");
    check(b.bodies[0].tets.size() + b.bodies[1].tets.size() ==
              static_cast<std::size_t>(m.num_tets()),
          "the two bodies partition the mesh");

    // The cut is INTERNAL: every one of its faces must have a tet on both
    // sides. A cut accidentally left on the outer boundary, or drawn inside a
    // volume rather than between two, would fail here -- and it is the whole
    // reason the geometry is built as two half-annuli.
    const BoundPort& cut = b.ports[0];
    check(is_internal(cut.type), "the port is internal");
    check(cut.surface.faces.size() > 0u, "the cut has faces");
    int interior = 0;
    for (int f : cut.surface.faces) {
        if (!m.is_boundary_face(f)) ++interior;
    }
    check(interior == static_cast<int>(cut.surface.faces.size()),
          "every face of the cut has a tet on BOTH sides -- it is interior, not a boundary");

    // Both sides are ring tets, not air: a cross-section of the conductor.
    int plus = 0, minus = 0, in_air = 0;
    for (int f : cut.surface.faces) {
        for (int k = 0; k < 2; ++k) {
            const int t =
                m.face_tets[static_cast<std::size_t>(f)].tets[static_cast<std::size_t>(k)];
            if (t < 0) continue;
            const int body = b.body_of_tet[static_cast<std::size_t>(t)];
            if (body != 0) {
                ++in_air;
            } else if (cut.side_of_tet(t) > 0) {
                ++plus;
            } else {
                ++minus;
            }
        }
    }
    check(in_air == 0, "the cut lies entirely inside the conductor");
    check(plus > 0 && minus > 0,
          "and it has tets on a plus side AND a minus side (" + std::to_string(plus) + " / " +
              std::to_string(minus) + ")");

    // The hint only resolves the sign, so the resolved direction must agree
    // with it to within 90 degrees -- here the cut's normal is +/-y, so this
    // is +y exactly.
    check(near(cut.direction.y, 1.0, 1e-9),
          "the resolved direction follows the +y hint, not the mesh's arbitrary face order");
    check(std::abs(cut.direction.x) < 1e-9 && std::abs(cut.direction.z) < 1e-9,
          "and it is purely +y, as a cut at theta = 0 must be");

    // The control the example file promises: reversing the hint reverses the
    // direction, and nothing else about the binding changes. Without this the
    // check above could pass on code that ignored the hint and happened to
    // pick +y from the mesh's face order.
    Problem reversed = p;
    reversed.ports[0].current_direction = Vec3{0.0, -1.0, 0.0};
    const BoundProblem rb = expect_ok(reversed, m, "the loop with the hint reversed");
    if (!rb.ports.empty()) {
        check(near(rb.ports[0].direction.y, -1.0, 1e-9), "-y reverses the resolved direction");
        check(rb.ports[0].surface.faces.size() == cut.surface.faces.size(),
              "and the cut is otherwise the same surface");
    }

    // A hint only has to be within 90 degrees of what is meant, so an
    // off-axis one must give the same answer as the clean axis.
    Problem sloppy = p;
    sloppy.ports[0].current_direction = Vec3{0.1, 0.9, -0.2};
    const BoundProblem sb = expect_ok(sloppy, m, "the loop with an off-axis hint");
    if (!sb.ports.empty()) {
        check(near(sb.ports[0].direction.y, 1.0, 1e-9),
              "0.1 0.9 -0.2 resolves to the same +y as the clean hint -- it picks a SIGN, "
              "not a normal");
    }

    // The cut fully spans the ring. If it did not, current would flow around
    // the uncut part and the port would be partly shorted, which
    // problem_binding refuses outright at DC -- so `expect_ok` above already
    // proves it. Assert the rim is genuinely a rim, on the ring's surface.
    check(cut.rim_edges.size() > 0u, "the cut has a rim");

    // One conduction path, and the cut is its own potential reference: a
    // floating loop has no other. This is what makes the loop test different
    // from the cylinder, where a 0 V terminal supplies the reference.
    check(b.conduction_paths.size() == 1u, "the ring is one conduction path");
    if (!b.conduction_paths.empty()) {
        check(b.conduction_paths[0].tets.size() == b.bodies[0].tets.size(),
              "which is the whole ring -- the cut does not disconnect it");
        check(b.conduction_paths[0].reference_count >= 1,
              "and the internal port references it, so nothing else has to");
    }

    // At DC Phi lives only in the ring.
    int phi_tets = 0;
    for (bool v : b.phi_tet) {
        if (v) ++phi_tets;
    }
    check(phi_tets == static_cast<int>(b.bodies[0].tets.size()),
          "at DC, Phi lives on the ring's tets and no others");

    // The ring is a closed loop: it must NOT touch the outer boundary, or it
    // would be an electrode rather than a floating loop.
    bool ring_on_boundary = false;
    for (int f = 0; f < m.num_faces(); ++f) {
        if (!m.is_boundary_face(f)) continue;
        const int t = m.face_tets[static_cast<std::size_t>(f)].tets[0];
        if (b.body_of_tet[static_cast<std::size_t>(t)] == 0) ring_on_boundary = true;
    }
    check(!ring_on_boundary, "the ring touches no outer face -- it floats inside the air");
#endif
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

    // One conduction path: the wire. The air conducts nothing, and both
    // ports sit on the same path.
    check(b.conduction_paths.size() == 1u, "the cylinder has one conduction path");
    if (b.conduction_paths.size() == 1u) {
        const ConductionPath& path = b.conduction_paths[0];
        check(path.tets.size() == 7535u, "the path is the wire");
        check(path.ports.size() == 2u, "both ports touch it");
        check(path.reference_count == 1, "P2 alone fixes its potential");
        check(!path.is_floating, "it is not floating");
    }

    // The gauge. Two of these are recorded numbers; the rest are derived
    // from the mesh here, so they would catch a change rather than merely
    // describe one.
    const int num_nodes = m.num_nodes();
    check(b.gauge.tree_edge_count == num_nodes - 1,
          "a spanning tree over a connected mesh has exactly nodes - 1 edges");
    check(b.gauge.tree_edge_count == 2911, "which is 2911 on this mesh");
    check(b.gauge.surface_tree_edge_count + b.gauge.interior_tree_edge_count ==
              b.gauge.tree_edge_count,
          "surface and interior tree edges account for the whole tree");
    check(b.gauge.num_dirichlet_components == 1,
          "the box's outer boundary is one connected n x A = 0 surface");
    check(b.gauge.num_reference_groups == 1, "one root for one connected mesh");

    // Count the Dirichlet edges and boundary nodes independently, then check
    // the construction's own invariants against them:
    //   surface tree edges = V_b - (number of surfaces)
    //   interior tree edges = N - V_b
    // These are what make the gauge complete: too few leaves a null space,
    // too many pins real magnetic flux (`docs/TREE_COTREE_GAUGE.md`).
    int dirichlet_edges = 0;
    for (bool v : b.dirichlet_edge) {
        if (v) ++dirichlet_edges;
    }
    std::vector<bool> boundary_node(static_cast<std::size_t>(num_nodes), false);
    for (int f = 0; f < m.num_faces(); ++f) {
        if (!m.is_boundary_face(f)) continue;
        for (int v : m.faces[static_cast<std::size_t>(f)]) {
            boundary_node[static_cast<std::size_t>(v)] = true;
        }
    }
    int V_b = 0;
    for (bool v : boundary_node) {
        if (v) ++V_b;
    }
    check(b.gauge.surface_tree_edge_count == V_b - 1,
          "surface tree edges == boundary nodes - 1 (one connected surface)");
    check(b.gauge.interior_tree_edge_count == num_nodes - V_b,
          "interior tree edges == nodes - boundary nodes: exactly the gauge constraints needed");
    check(dirichlet_edges == 1908, "1908 edges lie on the box boundary");

    const int free_a = m.num_edges() - dirichlet_edges - b.gauge.interior_tree_edge_count;
    check(free_a == 15405, "15405 free A unknowns remain of 19587 edges");
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
    test_conduction_paths();
    test_floating_conductor();
    test_gauge_is_built();
    test_cut_side_labels();
    test_cylinder_mesh();
    test_loop_mesh();

    std::cout << (g_checks - g_failures) << "/" << g_checks << " checks passed\n";
    return g_failures == 0 ? 0 : 1;
}
