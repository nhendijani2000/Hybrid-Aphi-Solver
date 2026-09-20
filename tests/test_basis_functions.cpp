// Minimal, dependency-free test runner -- see the comment in
// tests/CMakeLists.txt for why no external framework is fetched here.
//
// Sanity checks for the mixed-order basis pair locked in by
// docs/FORMULATION.md Sec 5.1: first-order Whitney/Nedelec edge elements for
// A, second-order (P2, 10-node) Lagrange nodal elements for Phi. These are
// element-level checks (one tet at a time), independent of Phase 02's own
// mesh-level exit criterion (CG = 0), which test_incidence.cpp covers.

#include <cmath>
#include <iostream>
#include <string>

#include "aphi_solver/basis_functions.hpp"
#include "aphi_solver/mesh.hpp"

using aphi_solver::compute_tet_geometry;
using aphi_solver::evaluate_barycentric;
using aphi_solver::Mesh;
using aphi_solver::p2_nodal_gradient;
using aphi_solver::p2_nodal_value;
using aphi_solver::TetGeometry;
using aphi_solver::Vec3;
using aphi_solver::whitney_edge_curl;
using aphi_solver::whitney_edge_value;

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

bool nearly(double a, double b, double tol = 1e-9) {
    return std::abs(a - b) < tol;
}

bool nearly(const Vec3& a, const Vec3& b, double tol = 1e-9) {
    return nearly(a.x, b.x, tol) && nearly(a.y, b.y, tol) && nearly(a.z, b.z, tol);
}

// A single reference tetrahedron, same as test_incidence.cpp's: nodes at
// the origin and the three unit axes.
Mesh make_reference_tet() {
    Mesh m;
    m.nodes = {Vec3(0, 0, 0), Vec3(1, 0, 0), Vec3(0, 1, 0), Vec3(0, 0, 1)};
    m.tets = {{0, 1, 2, 3}};
    m.build_topology();
    return m;
}

// A deliberately skewed, non-axis-aligned tetrahedron -- exercises the
// generic 4x4 solve in compute_tet_geometry on a tet where none of the
// shortcuts available for the reference tet (e.g. grad_L being a unit
// axis vector) apply.
Mesh make_skewed_tet() {
    Mesh m;
    m.nodes = {Vec3(0.3, -0.2, 0.1), Vec3(1.7, 0.4, -0.3), Vec3(0.1, 1.9, 0.6),
               Vec3(-0.4, 0.5, 2.1)};
    m.tets = {{0, 1, 2, 3}};
    m.build_topology();
    return m;
}

// Node position for one of the 10 P2 nodes of tet t in `mesh`: nodes 0..3
// are the tet's own vertices, nodes 4..9 are the midpoints of its 6 edges
// in kTetLocalEdgeVerts order (matching p2_nodal_value's own convention).
Vec3 p2_node_position(const Mesh& mesh, int t, int local_node) {
    const auto& verts = mesh.tets[static_cast<std::size_t>(t)];
    if (local_node < 4) {
        return mesh.nodes[static_cast<std::size_t>(verts[static_cast<std::size_t>(local_node)])];
    }
    const auto& [v0, v1] = aphi_solver::kTetLocalEdgeVerts[static_cast<std::size_t>(local_node - 4)];
    const Vec3& p0 = mesh.nodes[static_cast<std::size_t>(verts[static_cast<std::size_t>(v0)])];
    const Vec3& p1 = mesh.nodes[static_cast<std::size_t>(verts[static_cast<std::size_t>(v1)])];
    return (p0 + p1) * 0.5;
}

// Numerically integrates N_{local_edge}'s tangential circulation along the
// straight segment from mesh vertex m to vertex n of tet t: since N_ij(p(s))
// dotted with the (constant) tangent (p_n - p_m) is an affine function of
// the parameter s in [0,1] (shown in docs/FORMULATION.md's derivation --
// L_i(s), L_j(s) are affine and grad(L) is constant), Simpson's rule with a
// handful of points is exact up to floating point roundoff; more points are
// used anyway so this test does not depend on that being recognized.
double edge_circulation(const TetGeometry& g, int local_edge, const Vec3& pm, const Vec3& pn) {
    const Vec3 tangent = pn - pm;
    constexpr int kIntervals = 8;
    double sum = 0.0;
    for (int k = 0; k <= kIntervals; ++k) {
        double s = static_cast<double>(k) / kIntervals;
        Vec3 p = pm + tangent * s;
        std::array<double, 4> L = evaluate_barycentric(g, p);
        double f = whitney_edge_value(g, local_edge, L).dot(tangent);
        double weight;
        if (k == 0 || k == kIntervals) {
            weight = 1.0;
        } else if (k % 2 == 1) {
            weight = 4.0;
        } else {
            weight = 2.0;
        }
        sum += weight * f;
    }
    return sum * (1.0 / kIntervals) / 3.0;
}

// Finite-difference curl of the Whitney edge value field, at Cartesian
// point p, using central differences with step h. Cross-checks the
// closed-form whitney_edge_curl (constant over the tet) against the actual
// derivative of whitney_edge_value -- independent of whether that closed
// form was transcribed correctly.
Vec3 numeric_curl(const TetGeometry& g, int local_edge, const Vec3& p, double h) {
    auto value_at = [&](const Vec3& q) {
        return whitney_edge_value(g, local_edge, evaluate_barycentric(g, q));
    };
    Vec3 fx_plus = value_at(p + Vec3(h, 0, 0));
    Vec3 fx_minus = value_at(p - Vec3(h, 0, 0));
    Vec3 fy_plus = value_at(p + Vec3(0, h, 0));
    Vec3 fy_minus = value_at(p - Vec3(0, h, 0));
    Vec3 fz_plus = value_at(p + Vec3(0, 0, h));
    Vec3 fz_minus = value_at(p - Vec3(0, 0, h));

    double dFz_dy = (fy_plus.z - fy_minus.z) / (2 * h);
    double dFy_dz = (fz_plus.y - fz_minus.y) / (2 * h);
    double dFx_dz = (fz_plus.x - fz_minus.x) / (2 * h);
    double dFz_dx = (fx_plus.z - fx_minus.z) / (2 * h);
    double dFy_dx = (fx_plus.y - fx_minus.y) / (2 * h);
    double dFx_dy = (fy_plus.x - fy_minus.x) / (2 * h);

    return Vec3(dFz_dy - dFy_dz, dFx_dz - dFz_dx, dFy_dx - dFx_dy);
}

void run_checks_on_tet(const Mesh& mesh, const std::string& label) {
    const TetGeometry g = compute_tet_geometry(mesh, 0);
    const auto& verts = mesh.tets[0];

    // 1. Barycentric coordinates are Kronecker-delta at the tet's own
    //    vertices: L_i(p_j) = delta_ij.
    for (int j = 0; j < 4; ++j) {
        std::array<double, 4> L = evaluate_barycentric(g, mesh.nodes[static_cast<std::size_t>(verts[static_cast<std::size_t>(j)])]);
        for (int i = 0; i < 4; ++i) {
            double expected = (i == j) ? 1.0 : 0.0;
            check(nearly(L[static_cast<std::size_t>(i)], expected, 1e-9),
                  label + ": L_" + std::to_string(i) + " at vertex " + std::to_string(j) +
                      " == " + std::to_string(expected));
        }
    }

    // 2. Partition of unity: sum(L_i) == 1 everywhere, so sum(grad_L_i) ==
    //    0, at an arbitrary sample point (the tet's own centroid).
    Vec3 centroid = (mesh.nodes[static_cast<std::size_t>(verts[0])] + mesh.nodes[static_cast<std::size_t>(verts[1])] +
                      mesh.nodes[static_cast<std::size_t>(verts[2])] + mesh.nodes[static_cast<std::size_t>(verts[3])]) *
                     0.25;
    std::array<double, 4> Lc = evaluate_barycentric(g, centroid);
    double sum_L = Lc[0] + Lc[1] + Lc[2] + Lc[3];
    check(nearly(sum_L, 1.0), label + ": sum(L_i) == 1 at centroid");
    Vec3 sum_grad = g.grad_L[0] + g.grad_L[1] + g.grad_L[2] + g.grad_L[3];
    check(nearly(sum_grad, Vec3(0, 0, 0)), label + ": sum(grad L_i) == 0");
    for (int i = 0; i < 4; ++i) {
        check(nearly(Lc[static_cast<std::size_t>(i)], 0.25), label + ": centroid L_" + std::to_string(i) + " == 0.25");
    }

    // 3. P2 nodal basis is Kronecker-delta at the 10 P2 nodes (4 vertices +
    //    6 edge midpoints).
    for (int k = 0; k < 10; ++k) {
        Vec3 pk = p2_node_position(mesh, 0, k);
        std::array<double, 4> Lk = evaluate_barycentric(g, pk);
        for (int n = 0; n < 10; ++n) {
            double val = p2_nodal_value(n, Lk);
            double expected = (n == k) ? 1.0 : 0.0;
            check(nearly(val, expected, 1e-9), label + ": P2 node " + std::to_string(n) + " at node " +
                                                    std::to_string(k) + " == " + std::to_string(expected));
        }
    }

    // 4. P2 partition of unity at a non-node sample point (the centroid).
    double sum_p2 = 0.0;
    for (int n = 0; n < 10; ++n) {
        sum_p2 += p2_nodal_value(n, Lc);
    }
    check(nearly(sum_p2, 1.0), label + ": sum of all 10 P2 nodal values == 1 at centroid");

    // 5. P2 gradient consistency: finite-difference gradient of
    //    p2_nodal_value must match the closed-form p2_nodal_gradient, at an
    //    interior sample point.
    {
        const double h = 1e-6;
        Vec3 p0 = centroid;
        for (int n = 0; n < 10; ++n) {
            auto value_at = [&](const Vec3& q) { return p2_nodal_value(n, evaluate_barycentric(g, q)); };
            double dfdx = (value_at(p0 + Vec3(h, 0, 0)) - value_at(p0 - Vec3(h, 0, 0))) / (2 * h);
            double dfdy = (value_at(p0 + Vec3(0, h, 0)) - value_at(p0 - Vec3(0, h, 0))) / (2 * h);
            double dfdz = (value_at(p0 + Vec3(0, 0, h)) - value_at(p0 - Vec3(0, 0, h))) / (2 * h);
            Vec3 numeric_grad(dfdx, dfdy, dfdz);
            Vec3 analytic_grad = p2_nodal_gradient(g, n, Lc);
            check(nearly(numeric_grad, analytic_grad, 1e-5),
                  label + ": p2_nodal_gradient(" + std::to_string(n) + ") matches finite-difference gradient");
        }
    }

    // 6. Whitney edge circulation is Kronecker-delta over the tet's own 6
    //    edges: integral of N_j . dl along edge m equals delta_{jm}. This is
    //    the defining DOF-normalization property of Whitney edge elements
    //    (see docs/FORMULATION.md Sec 5.1) and is a strong, independent
    //    check on whitney_edge_value's formula.
    for (int j = 0; j < 6; ++j) {
        for (int m = 0; m < 6; ++m) {
            const auto& [vm, vn] = aphi_solver::kTetLocalEdgeVerts[static_cast<std::size_t>(m)];
            Vec3 pm = mesh.nodes[static_cast<std::size_t>(verts[static_cast<std::size_t>(vm)])];
            Vec3 pn = mesh.nodes[static_cast<std::size_t>(verts[static_cast<std::size_t>(vn)])];
            double circ = edge_circulation(g, j, pm, pn);
            double expected = (j == m) ? 1.0 : 0.0;
            check(nearly(circ, expected, 1e-7),
                  label + ": circulation of N_" + std::to_string(j) + " along edge " + std::to_string(m) +
                      " == " + std::to_string(expected));
        }
    }

    // 7. Whitney curl cross-check: the closed-form (constant) curl must
    //    match a finite-difference curl of the actual value field, at the
    //    centroid.
    for (int j = 0; j < 6; ++j) {
        Vec3 analytic = whitney_edge_curl(g, j);
        Vec3 numeric = numeric_curl(g, j, centroid, 1e-5);
        check(nearly(numeric, analytic, 1e-3),
              label + ": whitney_edge_curl(" + std::to_string(j) + ") matches finite-difference curl");
    }
}

}  // namespace

int main() {
    run_checks_on_tet(make_reference_tet(), "reference tet");
    run_checks_on_tet(make_skewed_tet(), "skewed tet");

    std::cout << g_checks - g_failures << "/" << g_checks << " checks passed.\n";
    return g_failures == 0 ? 0 : 1;
}
