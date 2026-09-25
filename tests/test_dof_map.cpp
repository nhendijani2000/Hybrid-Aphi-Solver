// Tests for the degree-of-freedom map.
//
// Step 6 of `docs/INPUT_FILE_PLAN.md` Sec. 6 -- the last thing
// before element matrices. The map is where a mistake is cheapest to catch:
// a mis-numbered unknown or a dropped edge sign produces a plausible matrix
// that solves to a wrong field, with nothing to distinguish it from a bad
// weak form afterwards.
//
// The properties worth asserting are exhaustiveness (every edge and every
// P2 node is classified exactly once, and the classes account for the mesh)
// and the local -> global map, checked against facts derived here rather
// than read back from the code under test.

#include <algorithm>
#include <cmath>
#include <iostream>
#include <set>
#include <string>
#include <vector>

#include "aphi_solver/dof_map.hpp"
#include "aphi_solver/gmsh_reader.hpp"

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

Mesh make_cube() {
    Mesh m;
    m.nodes = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
               {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};
    m.tets = {{0, 1, 2, 6}, {0, 2, 3, 6}, {0, 3, 7, 6},
              {0, 7, 4, 6}, {0, 4, 5, 6}, {0, 5, 1, 6}};
    m.build_topology();
    m.tet_tags = {1, 1, 1, 1, 1, 1};
    m.physical_names = {{3, 1, "metal"}, {2, 10, "bottom"}, {2, 11, "top"}};
    return m;
}

void tag_face(Mesh& m, int a, int b, int c, int tag) {
    TaggedFace f;
    f.nodes = {a, b, c};
    std::sort(f.nodes.begin(), f.nodes.end());
    f.tag = tag;
    m.tagged_boundary_faces.push_back(f);
}

Body make_body(const std::string& name, const std::string& volume, double sigma) {
    Body b;
    b.name = name;
    b.volume = volume;
    b.sigma = sigma;
    b.line = 1;
    return b;
}

Port make_port(const std::string& name, PortType type, const std::string& surface,
               double amplitude) {
    Port p;
    p.name = name;
    p.type = type;
    p.surface = {surface};
    p.amplitude = amplitude;
    p.line = 2;
    return p;
}

// One conductor cube, current in at the bottom, 0 V at the top.
BoundProblem bind_cube(Mesh& m) {
    tag_face(m, 0, 1, 2, 10);
    tag_face(m, 0, 2, 3, 10);
    tag_face(m, 4, 5, 6, 11);
    tag_face(m, 4, 6, 7, 11);

    Problem p;
    p.type = AnalysisType::DC;
    p.bodies = {make_body("B1", "metal", 5.8e7)};
    p.ports = {make_port("P1", PortType::BoundaryCurrent, "bottom", 1.0),
               make_port("P2", PortType::BoundaryVoltage, "top", 0.0)};
    return bind_to_mesh(p, m);
}

// ---------------------------------------------------------------------------

// Every edge falls into exactly one class, and the classes add up. This is
// the check that would catch an edge counted twice or missed entirely.
void test_edges_are_exhaustively_classified() {
    Mesh m = make_cube();
    const BoundProblem b = bind_cube(m);
    const DofMap d = build_dof_map(b, m);

    int free_edges = 0, dirichlet = 0, tree = 0;
    for (int e = 0; e < m.num_edges(); ++e) {
        switch (d.edge_state[static_cast<std::size_t>(e)]) {
            case EdgeDof::Free: ++free_edges; break;
            case EdgeDof::Dirichlet: ++dirichlet; break;
            case EdgeDof::Tree: ++tree; break;
        }
    }
    check(free_edges + dirichlet + tree == m.num_edges(),
          "every edge is Free, Dirichlet or Tree -- no edge counted twice or missed");
    check(d.num_a == free_edges, "num_a counts exactly the free edges");

    // Free edges carry indices 0..num_a-1, each exactly once.
    std::set<int> indices;
    for (int e = 0; e < m.num_edges(); ++e) {
        const int idx = d.edge_index[static_cast<std::size_t>(e)];
        if (d.edge_state[static_cast<std::size_t>(e)] == EdgeDof::Free) {
            check(idx >= 0 && idx < d.num_a, "a free edge has an index inside [0, num_a)");
            indices.insert(idx);
        } else {
            check(idx == -1, "an eliminated edge has no index");
        }
    }
    check(static_cast<int>(indices.size()) == d.num_a, "free edge indices are a bijection");

    // A Dirichlet edge is one the boundary condition zeroes; a Tree edge is
    // one the gauge zeroes. An edge that is both must be recorded as
    // Dirichlet, because that is the reason that survives a change of tree.
    for (int e = 0; e < m.num_edges(); ++e) {
        const std::size_t i = static_cast<std::size_t>(e);
        if (b.dirichlet_edge[i]) {
            check(d.edge_state[i] == EdgeDof::Dirichlet,
                  "an edge on an n x A = 0 surface is Dirichlet even when the tree also took it");
        }
    }
}

// Likewise for Phi, which is P2: one unknown per vertex AND one per edge
// midpoint.
void test_phi_nodes_are_exhaustively_classified() {
    Mesh m = make_cube();
    const BoundProblem b = bind_cube(m);
    const DofMap d = build_dof_map(b, m);

    check(d.num_p2_nodes == m.num_nodes() + m.num_edges(),
          "P2 gives one node per vertex and one per edge midpoint");

    int free_phi = 0, absent = 0, fixed = 0, port = 0, cut = 0;
    for (int i = 0; i < d.num_p2_nodes; ++i) {
        switch (d.phi_state[static_cast<std::size_t>(i)]) {
            case PhiDof::Free: ++free_phi; break;
            case PhiDof::Absent: ++absent; break;
            case PhiDof::Fixed: ++fixed; break;
            case PhiDof::Port: ++port; break;
            case PhiDof::Cut: ++cut; break;
        }
    }
    check(free_phi + absent + fixed + port + cut == d.num_p2_nodes,
          "every P2 node is classified exactly once");
    check(d.num_phi == free_phi, "num_phi counts exactly the free P2 nodes");
    check(absent == 0, "the cube is all conductor, so no P2 node is absent");

    // Both terminals are equipotentials, so their P2 nodes read a port
    // unknown rather than carrying one each.
    const std::size_t expected_port_nodes =
        b.ports[0].surface.nodes.size() + b.ports[0].surface.edges.size() +
        b.ports[1].surface.nodes.size() + b.ports[1].surface.edges.size();
    check(static_cast<std::size_t>(port) == expected_port_nodes,
          "a terminal's P2 nodes -- vertices AND edge midpoints -- all map to its port");
}

// The global numbering must be a partition: a, then phi, then V, with no
// index used twice and none skipped.
void test_global_numbering_is_a_partition() {
    Mesh m = make_cube();
    const BoundProblem b = bind_cube(m);
    const DofMap d = build_dof_map(b, m);

    check(d.num_total == d.num_a + d.num_phi + d.num_ports, "the three blocks sum to the total");
    check(d.num_ports == 2, "two ports, two unknowns");

    std::vector<int> seen(static_cast<std::size_t>(d.num_total), 0);
    for (int e = 0; e < m.num_edges(); ++e) {
        const int idx = d.edge_index[static_cast<std::size_t>(e)];
        if (idx >= 0) ++seen[static_cast<std::size_t>(idx)];
    }
    for (int i = 0; i < d.num_p2_nodes; ++i) {
        const int idx = d.phi_index[static_cast<std::size_t>(i)];
        if (idx >= 0) {
            check(idx >= d.num_a && idx < d.num_a + d.num_phi,
                  "a free Phi index lands in the Phi block");
            ++seen[static_cast<std::size_t>(idx)];
        }
    }
    for (int k = 0; k < d.num_ports; ++k) {
        const int idx = d.port_index[static_cast<std::size_t>(k)];
        check(idx >= d.num_a + d.num_phi && idx < d.num_total,
              "a port index lands last, after a and phi");
        ++seen[static_cast<std::size_t>(idx)];
    }
    bool all_once = true;
    for (int c : seen) {
        if (c != 1) all_once = false;
    }
    check(all_once, "every global index in [0, num_total) is claimed exactly once");
}

// The local -> global map. The edge sign is the part most likely to be
// dropped, and it is checked against the mesh directly rather than against
// the map's own record of it.
void test_local_dof_map() {
    Mesh m = make_cube();
    const BoundProblem b = bind_cube(m);
    const DofMap d = build_dof_map(b, m);

    for (int t = 0; t < m.num_tets(); ++t) {
        const TetDofs td = d.local_dofs(t, m, b);

        for (int le = 0; le < 6; ++le) {
            const int e = m.tet_edges[static_cast<std::size_t>(t)][static_cast<std::size_t>(le)];
            const double expected_sign = static_cast<double>(
                m.tet_edge_signs[static_cast<std::size_t>(t)][static_cast<std::size_t>(le)]);
            const DofEntry& entry = td.edge[static_cast<std::size_t>(le)];
            check(entry.index == d.edge_index[static_cast<std::size_t>(e)],
                  "a local edge maps to its global edge's index");
            check(entry.coeff == expected_sign,
                  "and carries Mesh::tet_edge_signs, without which most contributions flip");
        }

        // Local Phi node 4 + le must sit on local edge le -- the two
        // orderings line up by construction (mesh.hpp), and assembly relies
        // on it.
        for (int le = 0; le < 6; ++le) {
            const int e = m.tet_edges[static_cast<std::size_t>(t)][static_cast<std::size_t>(le)];
            const int p2 = d.edge_p2(e);
            const DofEntry& entry = td.phi[static_cast<std::size_t>(4 + le)];
            const PhiDof state = d.phi_state[static_cast<std::size_t>(p2)];
            const int expected = state == PhiDof::Free
                                     ? d.phi_index[static_cast<std::size_t>(p2)]
                                     : (state == PhiDof::Port
                                            ? d.port_index[static_cast<std::size_t>(
                                                  d.phi_port[static_cast<std::size_t>(p2)])]
                                            : -1);
            check(entry.index == expected,
                  "local Phi node 4+le is the midpoint of local edge le, mapped consistently");
        }
    }

    // At least one edge must actually be reversed, or the sign check above
    // is vacuous on this mesh.
    int reversed = 0;
    for (int t = 0; t < m.num_tets(); ++t) {
        for (int le = 0; le < 6; ++le) {
            if (m.tet_edge_signs[static_cast<std::size_t>(t)][static_cast<std::size_t>(le)] < 0) {
                ++reversed;
            }
        }
    }
    check(reversed > 0, "the cube really does have reversed local edges, so the sign check bites");
}

// An internal cut: the same node maps to the port's unknown from one side
// and to nothing (a prescribed zero) from the other. This is the only
// mapping that depends on the tet rather than only on the node.
void test_cut_nodes_depend_on_the_side() {
    Mesh m;
    m.nodes = {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, -1}};
    m.tets = {{0, 1, 2, 3}, {0, 1, 2, 4}};
    m.build_topology();
    m.tet_tags = {1, 1};
    m.physical_names = {{3, 1, "metal"}, {2, 20, "cut"}};
    tag_face(m, 0, 1, 2, 20);

    Problem p;
    p.type = AnalysisType::DC;
    p.bodies = {make_body("B1", "metal", 5.8e7)};
    Port cut = make_port("P1", PortType::InternalCurrent, "cut", 1.0);
    cut.current_direction = Vec3(0, 0, 1);  // plus side is the z > 0 tet
    p.ports = {cut};

    const BoundProblem b = bind_to_mesh(p, m);
    const DofMap d = build_dof_map(b, m);

    const int plus_tet = b.ports[0].plus_side_tet[0];
    const int minus_tet = plus_tet == 0 ? 1 : 0;
    const int port_dof = d.port_index[0];

    const TetDofs from_plus = d.local_dofs(plus_tet, m, b);
    const TetDofs from_minus = d.local_dofs(minus_tet, m, b);

    // Find a local Phi node that is a cut node in both tets: the cut's
    // vertices 0, 1, 2 are shared.
    int checked = 0;
    for (int ln = 0; ln < 4; ++ln) {
        const int v = m.tets[static_cast<std::size_t>(plus_tet)][static_cast<std::size_t>(ln)];
        if (v != 0 && v != 1 && v != 2) continue;
        const int p2 = d.vertex_p2(v);
        check(d.phi_state[static_cast<std::size_t>(p2)] == PhiDof::Cut,
              "a cut vertex is classified Cut, not Free");
        check(from_plus.phi[static_cast<std::size_t>(ln)].index == port_dof,
              "seen from the plus side, a cut node reads the port's unknown");
        ++checked;
    }
    check(checked == 3, "all three cut vertices were checked");

    for (int ln = 0; ln < 4; ++ln) {
        const int v = m.tets[static_cast<std::size_t>(minus_tet)][static_cast<std::size_t>(ln)];
        if (v != 0 && v != 1 && v != 2) continue;
        check(from_minus.phi[static_cast<std::size_t>(ln)].index == -1,
              "seen from the grounded minus side, the same node reads a prescribed zero");
    }

    // The apex of each tet is not on the cut, so it is an ordinary unknown.
    check(d.phi_state[static_cast<std::size_t>(d.vertex_p2(3))] == PhiDof::Free,
          "a node off the cut is free");
    check(d.phi_state[static_cast<std::size_t>(d.vertex_p2(4))] == PhiDof::Free,
          "on both sides");
}

// At DC, Phi lives only on the conductors, so P2 nodes touched by no
// conducting tet must be Absent -- not free unknowns with empty rows.
void test_phi_absent_outside_conductors() {
    Mesh m = make_cube();
    tag_face(m, 0, 1, 2, 10);
    tag_face(m, 0, 2, 3, 10);
    tag_face(m, 4, 5, 6, 11);
    tag_face(m, 4, 6, 7, 11);
    m.tet_tags = {1, 1, 1, 2, 2, 2};
    m.physical_names.push_back({3, 2, "air"});

    Problem p;
    p.type = AnalysisType::DC;
    p.bodies = {make_body("B1", "metal", 5.8e7), make_body("B2", "air", 0.0)};
    p.ports = {make_port("P1", PortType::BoundaryVoltage, "bottom", 0.0)};

    const BoundProblem b = bind_to_mesh(p, m);
    const DofMap d = build_dof_map(b, m);

    int absent = 0;
    for (int i = 0; i < d.num_p2_nodes; ++i) {
        if (d.phi_state[static_cast<std::size_t>(i)] == PhiDof::Absent) ++absent;
    }
    check(absent > 0, "at DC some P2 nodes lie outside the conductor and are Absent");

    // Derive the expected present-set independently.
    std::set<int> expected;
    for (int t = 0; t < m.num_tets(); ++t) {
        if (!b.phi_tet[static_cast<std::size_t>(t)]) continue;
        for (int v : m.tets[static_cast<std::size_t>(t)]) expected.insert(d.vertex_p2(v));
        for (int e : m.tet_edges[static_cast<std::size_t>(t)]) expected.insert(d.edge_p2(e));
    }
    int present = 0;
    for (int i = 0; i < d.num_p2_nodes; ++i) {
        if (d.phi_state[static_cast<std::size_t>(i)] != PhiDof::Absent) ++present;
    }
    check(static_cast<std::size_t>(present) == expected.size(),
          "exactly the P2 nodes of conducting tets are present");
}

void test_ports_are_free_or_fixed() {
    Mesh m = make_cube();
    const BoundProblem b = bind_cube(m);
    const DofMap d = build_dof_map(b, m);

    check(!d.port_is_fixed[0], "a current port's V is an unknown -- the solve returns it");
    check(d.port_is_fixed[1], "a voltage port's V is prescribed");
    check(d.port_value[1] == std::complex<double>(0.0, 0.0), "P2 is prescribed 0 V");

    // Both still get an index: a voltage port keeps its row, because that
    // row IS the port current, recovered as a residual rather than by
    // differentiating the field.
    check(d.port_index[0] >= 0 && d.port_index[1] >= 0,
          "both ports have a global index, fixed or not");
}

// The real fixture, where the A counts are already known from the gauge.
void test_cylinder() {
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
    p.bodies = {make_body("B1", "wire", 5.8e7), make_body("B2", "air", 0.0)};
    p.ports = {make_port("P1", PortType::BoundaryCurrent, "wire_bottom", 1.0),
               make_port("P2", PortType::BoundaryVoltage, "wire_top", 0.0)};

    const BoundProblem b = bind_to_mesh(p, m);
    const DofMap d = build_dof_map(b, m);

    // A: the same 15405 the gauge predicted, reached by a different route
    // -- counting classified edges rather than subtracting counts.
    check(d.num_a == 15405, "15405 free A unknowns, agreeing with the gauge's own arithmetic");

    int dirichlet = 0, tree = 0;
    for (int e = 0; e < m.num_edges(); ++e) {
        if (d.edge_state[static_cast<std::size_t>(e)] == EdgeDof::Dirichlet) ++dirichlet;
        if (d.edge_state[static_cast<std::size_t>(e)] == EdgeDof::Tree) ++tree;
    }
    check(dirichlet == 1908, "1908 Dirichlet edges");
    check(tree == b.gauge.interior_tree_edge_count,
          "the Tree class holds exactly the interior tree edges -- the surface ones are "
          "already Dirichlet and must not be counted twice");
    check(d.num_a + dirichlet + tree == m.num_edges(), "the three classes cover every edge");

    // Phi at DC lives on the wire only. Both terminals are equipotentials.
    const std::size_t terminal_nodes =
        b.ports[0].surface.nodes.size() + b.ports[0].surface.edges.size() +
        b.ports[1].surface.nodes.size() + b.ports[1].surface.edges.size();
    check(terminal_nodes == 2u * (74u + 195u), "each terminal has 74 + 195 = 269 P2 nodes");

    int port_nodes = 0, absent = 0;
    for (int i = 0; i < d.num_p2_nodes; ++i) {
        if (d.phi_state[static_cast<std::size_t>(i)] == PhiDof::Port) ++port_nodes;
        if (d.phi_state[static_cast<std::size_t>(i)] == PhiDof::Absent) ++absent;
    }
    check(static_cast<std::size_t>(port_nodes) == terminal_nodes,
          "every terminal P2 node is slaved to its port");
    check(absent > 0, "the air's P2 nodes are absent at DC");

    check(d.num_total == d.num_a + d.num_phi + 2, "two port unknowns on top of a and phi");
    std::cout << "  cylinder DOFs: a = " << d.num_a << ", phi = " << d.num_phi
              << ", ports = " << d.num_ports << ", total = " << d.num_total << "\n";
#endif
}

}  // namespace

int main() {
    test_edges_are_exhaustively_classified();
    test_phi_nodes_are_exhaustively_classified();
    test_global_numbering_is_a_partition();
    test_local_dof_map();
    test_cut_nodes_depend_on_the_side();
    test_phi_absent_outside_conductors();
    test_ports_are_free_or_fixed();
    test_cylinder();

    std::cout << (g_checks - g_failures) << "/" << g_checks << " checks passed\n";
    return g_failures == 0 ? 0 : 1;
}
