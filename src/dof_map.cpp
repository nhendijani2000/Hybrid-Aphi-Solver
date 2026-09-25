#include "aphi_solver/dof_map.hpp"

#include <stdexcept>
#include <string>

namespace aphi_solver {

DofMap build_dof_map(const BoundProblem& bound, const Mesh& mesh) {
    DofMap map;
    map.num_nodes = mesh.num_nodes();
    map.num_p2_nodes = mesh.num_nodes() + mesh.num_edges();

    // ---- A: one unknown per edge that is neither prescribed nor gauged ----
    //
    // An edge on an n x A = 0 surface that the tree also chose is BOTH. It
    // is recorded as Dirichlet, because that is the reason that would still
    // apply if the gauge changed: the boundary condition is physics, the
    // tree is a choice.
    map.edge_state.assign(static_cast<std::size_t>(mesh.num_edges()), EdgeDof::Free);
    map.edge_index.assign(static_cast<std::size_t>(mesh.num_edges()), -1);
    for (int e = 0; e < mesh.num_edges(); ++e) {
        const std::size_t i = static_cast<std::size_t>(e);
        if (bound.dirichlet_edge[i]) {
            map.edge_state[i] = EdgeDof::Dirichlet;
        } else if (bound.gauge.is_tree_edge[i]) {
            map.edge_state[i] = EdgeDof::Tree;
        } else {
            map.edge_index[i] = map.num_a++;
        }
    }

    // ---- Phi: P2, so one unknown per vertex and one per edge midpoint ----
    map.phi_state.assign(static_cast<std::size_t>(map.num_p2_nodes), PhiDof::Absent);
    map.phi_index.assign(static_cast<std::size_t>(map.num_p2_nodes), -1);
    map.phi_port.assign(static_cast<std::size_t>(map.num_p2_nodes), -1);
    map.phi_fixed_value.assign(static_cast<std::size_t>(map.num_p2_nodes),
                               std::complex<double>(0.0, 0.0));

    // Present wherever some tet of Phi's support touches it. At full wave
    // that is everything; at DC only the conductors, where Phi's equation is
    // not identically empty.
    std::vector<bool> present(static_cast<std::size_t>(map.num_p2_nodes), false);
    for (int t = 0; t < mesh.num_tets(); ++t) {
        if (!bound.phi_tet[static_cast<std::size_t>(t)]) continue;
        for (int v : mesh.tets[static_cast<std::size_t>(t)]) {
            present[static_cast<std::size_t>(map.vertex_p2(v))] = true;
        }
        for (int e : mesh.tet_edges[static_cast<std::size_t>(t)]) {
            present[static_cast<std::size_t>(map.edge_p2(e))] = true;
        }
    }

    // Port terminals and cuts claim their nodes before anything else does.
    // A boundary terminal is an equipotential reading the port's unknown; a
    // cut reads 0 or that unknown depending on the side, which only the tet
    // asking can say (local_dofs below).
    for (std::size_t k = 0; k < bound.ports.size(); ++k) {
        const BoundPort& port = bound.ports[k];
        const PhiDof state = port.is_internal() ? PhiDof::Cut : PhiDof::Port;
        for (int v : port.surface.nodes) {
            const std::size_t i = static_cast<std::size_t>(map.vertex_p2(v));
            map.phi_state[i] = state;
            map.phi_port[i] = static_cast<int>(k);
        }
        for (int e : port.surface.edges) {
            const std::size_t i = static_cast<std::size_t>(map.edge_p2(e));
            map.phi_state[i] = state;
            map.phi_port[i] = static_cast<int>(k);
        }
    }

    // A floating conductor's pinned node, if it was not already claimed.
    for (const ConductionPath& path : bound.conduction_paths) {
        if (!path.is_floating || path.pin_node < 0) continue;
        const std::size_t i = static_cast<std::size_t>(map.vertex_p2(path.pin_node));
        if (map.phi_state[i] == PhiDof::Port || map.phi_state[i] == PhiDof::Cut) continue;
        map.phi_state[i] = PhiDof::Fixed;
        map.phi_fixed_value[i] = std::complex<double>(0.0, 0.0);
    }

    for (int i = 0; i < map.num_p2_nodes; ++i) {
        const std::size_t u = static_cast<std::size_t>(i);
        if (!present[u]) {
            map.phi_state[u] = PhiDof::Absent;
            continue;
        }
        if (map.phi_state[u] == PhiDof::Absent) map.phi_state[u] = PhiDof::Free;
        if (map.phi_state[u] == PhiDof::Free) map.phi_index[u] = map.num_phi++;
    }

    // ---- ports: one unknown each, free if current-driven --------------
    //
    // A voltage port keeps an unknown rather than being folded away as a
    // fixed potential, because its row IS the port current: after the solve
    // the residual of that row is I, read from the equation rather than by
    // differentiating the field.
    for (const BoundPort& port : bound.ports) {
        map.port_index.push_back(-1);
        const bool fixed = !is_current_driven(port.type);
        map.port_is_fixed.push_back(fixed);
        map.port_value.push_back(fixed ? port.amplitude : std::complex<double>(0.0, 0.0));
    }

    // ---- global numbering: a, then phi, then V ------------------------
    // Edges were numbered from 0 as they were classified, so the a block is
    // already in place; the phi indices assigned above start from 0 too and
    // are shifted past it here.
    for (int i = 0; i < map.num_p2_nodes; ++i) {
        const std::size_t u = static_cast<std::size_t>(i);
        if (map.phi_index[u] >= 0) map.phi_index[u] += map.num_a;
    }
    map.num_ports = static_cast<int>(bound.ports.size());
    for (int k = 0; k < map.num_ports; ++k) {
        map.port_index[static_cast<std::size_t>(k)] = map.num_a + map.num_phi + k;
    }
    map.num_total = map.num_a + map.num_phi + map.num_ports;

    return map;
}

TetDofs DofMap::local_dofs(int tet, const Mesh& mesh, const BoundProblem& bound) const {
    TetDofs out;
    const std::size_t t = static_cast<std::size_t>(tet);

    // A: the edge sign is not optional. A tet's local vertex order comes
    // from the mesh file and is arbitrary, so roughly 58 % of (tet, local
    // edge) pairs run against the global edge's direction; scattering
    // without this sign gets more than half of them backwards.
    for (int le = 0; le < 6; ++le) {
        const int e = mesh.tet_edges[t][static_cast<std::size_t>(le)];
        out.edge[static_cast<std::size_t>(le)] = {
            edge_index[static_cast<std::size_t>(e)],
            static_cast<double>(mesh.tet_edge_signs[t][static_cast<std::size_t>(le)])};
    }

    // Phi: local nodes 0..3 are the tet's vertices, 4..9 the midpoints of
    // its local edges 0..5 -- the same order, so local Phi node 4 + le sits
    // on local edge le.
    std::array<int, 10> p2{};
    for (int v = 0; v < 4; ++v) p2[static_cast<std::size_t>(v)] = vertex_p2(mesh.tets[t][static_cast<std::size_t>(v)]);
    for (int le = 0; le < 6; ++le) {
        p2[static_cast<std::size_t>(4 + le)] = edge_p2(mesh.tet_edges[t][static_cast<std::size_t>(le)]);
    }

    for (int ln = 0; ln < 10; ++ln) {
        const std::size_t l = static_cast<std::size_t>(ln);
        const std::size_t i = static_cast<std::size_t>(p2[l]);
        switch (phi_state[i]) {
            case PhiDof::Free:
                out.phi[l] = {phi_index[i], 1.0};
                break;
            case PhiDof::Port:
                out.phi[l] = {port_index[static_cast<std::size_t>(phi_port[i])], 1.0};
                break;
            case PhiDof::Cut: {
                // The side decides: the plus side reads the port's unknown,
                // the grounded minus side reads a prescribed zero. This is
                // the one mapping that depends on the tet and not only on
                // the node, which is why the map is built per tet.
                const BoundPort& port = bound.ports[static_cast<std::size_t>(phi_port[i])];
                if (port.side_of_tet(tet) > 0) {
                    out.phi[l] = {port_index[static_cast<std::size_t>(phi_port[i])], 1.0};
                }
                break;
            }
            case PhiDof::Fixed:
                out.phi_fixed[l] = phi_fixed_value[i];
                break;
            case PhiDof::Absent:
                break;
        }
    }
    return out;
}

}  // namespace aphi_solver
