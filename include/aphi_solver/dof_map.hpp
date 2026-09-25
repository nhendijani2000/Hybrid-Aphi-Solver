#pragma once

#include <array>
#include <complex>
#include <vector>

#include "aphi_solver/mesh.hpp"
#include "aphi_solver/problem_binding.hpp"

namespace aphi_solver {

/// The degree-of-freedom map: which unknowns exist, and how a tet's local
/// contributions reach them.
///
/// **Unknown vector layout.** Port unknowns go last because each couples to
/// every DOF of its terminal, so their rows and columns are comparatively
/// dense; keeping them out of the way leaves the A/Phi structure clean.
///
/// ```
/// x = [ a  (free edges) | phi  (free P2 nodes) | V  (one per port) ]
///       0 .. num_a         num_a ..              num_a + num_phi ..
/// ```
///
/// **Phi is P2**, so it has one unknown per mesh *vertex* and one per mesh
/// *edge midpoint* (`docs/FORMULATION.md` Sec. 5.2). Those are indexed
/// together as "P2 nodes": vertex `v` is P2 node `v`, and the midpoint of
/// edge `e` is P2 node `num_nodes + e`. `num_p2_nodes` is the sum.

/// How one A (edge) unknown is treated. Two of the three are eliminated,
/// but for entirely different reasons, and the distinction is worth keeping:
/// the boundary condition is physics and the gauge is an arbitrary choice.
/// Pick a different tree and you get a different A with the *same* B; drop
/// a Dirichlet edge and you are solving a different problem.
enum class EdgeDof {
    Free,       ///< a genuine unknown
    Dirichlet,  ///< a = 0 because n x A = 0 on a surface this edge lies in
    Tree        ///< a = 0 because the tree-cotree gauge says so
};

/// How one Phi (P2 node) unknown is treated.
enum class PhiDof {
    Free,    ///< a genuine unknown
    Absent,  ///< outside Phi's support -- at DC, anywhere sigma == 0
    Fixed,   ///< prescribed: a floating conductor's pinned node
    Port,    ///< on a port terminal: reads that port's unknown
    Cut      ///< on an internal cut: reads 0 or the port's unknown,
             ///< depending which side the tet asking is on
};

/// One entry of the local -> global map: which global unknown a local DOF
/// contributes to, and with what coefficient.
struct DofEntry {
    int index = -1;      ///< -1 means "eliminated": contributes to nothing
    double coeff = 0.0;  ///< the edge orientation sign, or 1 for Phi
};

/// A tet's 16 local DOFs resolved against the global numbering: 6 edges
/// then 10 P2 nodes, in `kTetLocalEdgeVerts` order for both (local Phi node
/// 0..3 are the tet's vertices, 4..9 the midpoints of local edges 0..5).
///
/// `phi_fixed` carries the prescribed value of any local Phi DOF that is
/// eliminated with a non-zero value, which the right-hand side needs. An
/// eliminated DOF whose value is zero leaves it zero, which is the common
/// case -- a Dirichlet edge, a grounded cut side, a pinned node.
struct TetDofs {
    std::array<DofEntry, 6> edge;
    std::array<DofEntry, 10> phi;
    std::array<std::complex<double>, 10> phi_fixed{};
};

/// The map itself. A plain struct of parallel arrays, matching the rest of
/// this project: the data is what the assembly loop reads, and hiding it
/// behind accessors would buy nothing.
struct DofMap {
    int num_nodes = 0;     ///< mesh vertices, so P2 midpoints start here
    int num_p2_nodes = 0;  ///< num_nodes + num_edges

    int num_a = 0;      ///< free edge unknowns
    int num_phi = 0;    ///< free P2-node unknowns
    int num_ports = 0;  ///< one per port, free or fixed
    int num_total = 0;  ///< num_a + num_phi + num_ports

    /// Per mesh edge.
    std::vector<EdgeDof> edge_state;
    std::vector<int> edge_index;  ///< global index, or -1 if eliminated

    /// Per P2 node.
    std::vector<PhiDof> phi_state;
    std::vector<int> phi_index;                        ///< global index, or -1
    std::vector<int> phi_port;                         ///< owning port, or -1
    std::vector<std::complex<double>> phi_fixed_value; ///< for PhiDof::Fixed

    /// Per port, in `BoundProblem::ports` order.
    std::vector<int> port_index;
    std::vector<bool> port_is_fixed;                  ///< true for a voltage port
    std::vector<std::complex<double>> port_value;     ///< prescribed, if fixed

    int vertex_p2(int v) const { return v; }
    int edge_p2(int e) const { return num_nodes + e; }

    /// The local -> global map for one tet, built on demand rather than
    /// stored: 16 entries per tet would be a large array read once, and a
    /// cut node's mapping depends on which side the *tet* is on, so it
    /// cannot be stored per node anyway.
    TetDofs local_dofs(int tet, const Mesh& mesh, const BoundProblem& bound) const;
};

/// Builds the map. Every edge becomes exactly one of Free / Dirichlet /
/// Tree, and every P2 node exactly one PhiDof, so the counts below account
/// for the whole mesh -- which is what the tests check.
DofMap build_dof_map(const BoundProblem& bound, const Mesh& mesh);

}  // namespace aphi_solver
