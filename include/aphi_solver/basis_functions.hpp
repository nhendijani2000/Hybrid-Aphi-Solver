#pragma once

#include <array>

#include "aphi_solver/mesh.hpp"

namespace aphi_solver {

/// Barycentric-coordinate geometry for one tetrahedron: since each L_i is
/// linear (L_i(p) = a[i] + grad_L[i].dot(p)), grad_L is constant over the
/// element and only needs computing once per tet, not once per evaluation
/// point. `volume` is the unsigned tet volume (used later for quadrature
/// weight scaling in Phase 04).
struct TetGeometry {
    std::array<double, 4> a{};
    std::array<Vec3, 4> grad_L{};
    double volume = 0.0;
};

/// Computes TetGeometry for tet t of `mesh`, by solving for each L_i's
/// affine coefficients directly (a 4x4 linear solve, generic and correct
/// for any non-degenerate tet regardless of vertex-ordering orientation --
/// see basis_functions.cpp for the derivation). Standard finite-element
/// geometry, e.g. Jin (2014) Sec. 5.3.2; not specific to any implementation.
TetGeometry compute_tet_geometry(const Mesh& mesh, int t);

/// Evaluates all 4 barycentric coordinates at Cartesian point `p`, using a
/// tet's precomputed TetGeometry.
std::array<double, 4> evaluate_barycentric(const TetGeometry& g, const Vec3& p);

/// First-order (lowest-order) Whitney/Nedelec edge basis function for A --
/// docs/FORMULATION.md Sec 5.1's final basis-function decision. `local_edge`
/// is 0..5, indexing into Mesh::kTetLocalEdgeVerts (i.e. tet_edges[t][local_edge]
/// gives the corresponding global edge/DOF index). `L` is the barycentric
/// coordinates at the evaluation point (e.g. from evaluate_barycentric, or
/// straight from a quadrature rule's own barycentric points).
///
///   N_ij(p) = L_i(p) * grad(L_j) - L_j(p) * grad(L_i)
///
/// verified in this project against Bossavit (1998) Ch. 5 and the user's own
/// prior implementation (BasisFunction.docx / EdModel.cpp).
Vec3 whitney_edge_value(const TetGeometry& g, int local_edge, const std::array<double, 4>& L);

/// Curl of the first-order Whitney edge basis function -- constant over the
/// tet (no evaluation point needed): curl(N_ij) = 2 * grad(L_i) x grad(L_j).
Vec3 whitney_edge_curl(const TetGeometry& g, int local_edge);

/// Second-order (P2, 10-node) nodal basis function for Phi --
/// docs/FORMULATION.md Sec 5.1's final basis-function decision, verified
/// against J.-M. Jin (2014), Ch. 5, Eq. (5.52) / Fig. 5.3. `local_node` is
/// 0..9: 0..3 are the tet's 4 vertices, 4..9 are its 6 edge midpoints in
/// Mesh::kTetLocalEdgeVerts order (so local_node - 4 == the local edge whose
/// midpoint this node is).
///
///   vertex node i   (i = 0..3): N_i = (2*L_i - 1) * L_i
///   edge midpoint e (e = 0..5): N_{4+e} = 4 * L_v0(e) * L_v1(e)
double p2_nodal_value(int local_node, const std::array<double, 4>& L);

/// Gradient of the P2 nodal basis function at the point whose barycentric
/// coordinates are `L` -- linear over the tet (not constant, unlike the
/// first-order edge case above), via the chain rule through grad(L_i).
Vec3 p2_nodal_gradient(const TetGeometry& g, int local_node, const std::array<double, 4>& L);

}  // namespace aphi_solver
