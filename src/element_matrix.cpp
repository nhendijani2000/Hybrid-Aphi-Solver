#include "aphi_solver/element_matrix.hpp"

#include <array>

#include "aphi_solver/quadrature.hpp"

namespace aphi_solver {

ElementCoefficients element_coefficients(const BoundBody& body, double omega) {
    // mu_0 and eps_0 in SI. The mesh is in metres by the time assembly runs
    // (scale_mesh_to_metres), so no length factor belongs here.
    constexpr double mu0 = 4.0e-7 * 3.14159265358979323846;
    constexpr double eps0 = 8.8541878128e-12;

    const double mu = body.mu_r * mu0;
    const double eps = body.eps_r * eps0;

    ElementCoefficients c;
    c.nu = 1.0 / mu;
    // alpha = j*w*sigma - w^2*eps, written so that omega == 0 gives exactly
    // zero rather than a rounding of it.
    c.alpha = std::complex<double>(-omega * omega * eps, omega * body.sigma);
    c.beta = std::complex<double>(body.sigma, omega * eps);
    return c;
}

void kernel_AA(const TetGeometry& g, const ElementCoefficients& c, std::complex<double>* out) {
    // --- curl-curl: constant integrand, so no quadrature loop ------------
    //
    // curl W_i = 2 grad(La) x grad(Lb) is constant over the tet, so this
    // term is exactly nu * V * (curl W_i) . (curl W_j). Computing it by
    // quadrature would give the same answer -- tests/test_element_matrix.cpp
    // checks that it does -- but this is the one block evaluated for every
    // tet in every regime, so it is worth not looping over.
    std::array<Vec3, 6> curl{};
    for (int i = 0; i < 6; ++i) curl[static_cast<std::size_t>(i)] = whitney_edge_curl(g, i);

    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            const double cc = curl[static_cast<std::size_t>(i)].dot(curl[static_cast<std::size_t>(j)]);
            out[static_cast<std::size_t>(i * 6 + j)] = c.nu * g.volume * cc;
        }
    }

    // --- mass: quadratic, exact under the degree-2 rule -------------------
    //
    // Skipped entirely when alpha is zero, which is every DC solve: there
    // the A equation is magnetostatics and has no mass term at all
    // (FORMULATION.md Sec. 2). Not an optimisation so much as declining to
    // add exact zeros.
    if (c.alpha == std::complex<double>(0.0, 0.0)) return;

    const QuadratureRule q = tet_rule_degree_2();
    for (int k = 0; k < q.count; ++k) {
        const std::array<double, 4>& L = q.points[static_cast<std::size_t>(k)];
        const double w = q.weights[static_cast<std::size_t>(k)] * g.volume;

        std::array<Vec3, 6> value{};
        for (int i = 0; i < 6; ++i) value[static_cast<std::size_t>(i)] = whitney_edge_value(g, i, L);

        for (int i = 0; i < 6; ++i) {
            for (int j = 0; j < 6; ++j) {
                const double vv =
                    value[static_cast<std::size_t>(i)].dot(value[static_cast<std::size_t>(j)]);
                out[static_cast<std::size_t>(i * 6 + j)] += c.alpha * (w * vv);
            }
        }
    }
}

void kernel_APhi(const TetGeometry& g, const ElementCoefficients& c, std::complex<double>* out) {
    for (int i = 0; i < 60; ++i) out[static_cast<std::size_t>(i)] = std::complex<double>(0.0, 0.0);
    if (c.beta == std::complex<double>(0.0, 0.0)) return;

    const QuadratureRule q = tet_rule_degree_2();
    for (int k = 0; k < q.count; ++k) {
        const std::array<double, 4>& L = q.points[static_cast<std::size_t>(k)];
        const double w = q.weights[static_cast<std::size_t>(k)] * g.volume;

        std::array<Vec3, 6> edge{};
        for (int i = 0; i < 6; ++i) edge[static_cast<std::size_t>(i)] = whitney_edge_value(g, i, L);
        std::array<Vec3, 10> grad{};
        for (int b = 0; b < 10; ++b) grad[static_cast<std::size_t>(b)] = p2_nodal_gradient(g, b, L);

        for (int i = 0; i < 6; ++i) {
            for (int b = 0; b < 10; ++b) {
                const double wg =
                    edge[static_cast<std::size_t>(i)].dot(grad[static_cast<std::size_t>(b)]);
                out[static_cast<std::size_t>(i * 10 + b)] += c.beta * (w * wg);
            }
        }
    }
}

void kernel_PhiPhi(const TetGeometry& g, const ElementCoefficients& c, std::complex<double>* out) {
    for (int i = 0; i < 100; ++i) out[static_cast<std::size_t>(i)] = std::complex<double>(0.0, 0.0);
    if (c.beta == std::complex<double>(0.0, 0.0)) return;

    const QuadratureRule q = tet_rule_degree_2();
    for (int k = 0; k < q.count; ++k) {
        const std::array<double, 4>& L = q.points[static_cast<std::size_t>(k)];
        const double w = q.weights[static_cast<std::size_t>(k)] * g.volume;

        std::array<Vec3, 10> grad{};
        for (int a = 0; a < 10; ++a) grad[static_cast<std::size_t>(a)] = p2_nodal_gradient(g, a, L);

        for (int a = 0; a < 10; ++a) {
            for (int b = 0; b < 10; ++b) {
                const double gg =
                    grad[static_cast<std::size_t>(a)].dot(grad[static_cast<std::size_t>(b)]);
                out[static_cast<std::size_t>(a * 10 + b)] += c.beta * (w * gg);
            }
        }
    }
}

}  // namespace aphi_solver
