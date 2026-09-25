// Tests for the element coefficients and the A-A element block.
//
// Steps 2 and 3 of `docs/ASSEMBLY_PLAN.md` Sec. 8. The kernels are pure
// functions of a TetGeometry, so everything here is checked against numbers
// worked out by hand rather than against another implementation.

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <iostream>
#include <string>
#include <vector>

#include "aphi_solver/element_matrix.hpp"
#include "aphi_solver/quadrature.hpp"

using namespace aphi_solver;
using Complex = std::complex<double>;

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
bool near(Complex a, Complex b, double tol = 1e-12) { return std::abs(a - b) <= tol; }

/// The reference tet (0,0,0),(1,0,0),(0,1,0),(0,0,1), optionally scaled.
/// Built through a one-tet Mesh so the already-tested geometry code is what
/// produces it -- the kernel itself needs no mesh.
TetGeometry unit_tet(double s = 1.0) {
    Mesh m;
    m.nodes = {{0, 0, 0}, {s, 0, 0}, {0, s, 0}, {0, 0, s}};
    m.tets = {{0, 1, 2, 3}};
    m.build_topology();
    return compute_tet_geometry(m, 0);
}

BoundBody make_body(double sigma, double eps_r = 1.0, double mu_r = 1.0) {
    BoundBody b;
    b.sigma = sigma;
    b.eps_r = eps_r;
    b.mu_r = mu_r;
    return b;
}

// ---------------------------------------------------------------------------

void test_coefficients() {
    const BoundBody copper = make_body(5.8e7);

    // At DC the mass coefficient must be EXACTLY zero, not a rounding of
    // it: that is the decoupling of FORMULATION.md Sec. 2 falling out of the
    // coefficients rather than needing a branch.
    const ElementCoefficients dc = element_coefficients(copper, 0.0);
    check(dc.alpha == Complex(0.0, 0.0), "at omega = 0, alpha is exactly zero");
    check(dc.beta == Complex(5.8e7, 0.0), "at omega = 0, beta is exactly sigma");
    check(dc.nu > 0.0, "nu is positive");

    // nu = 1/mu, and mu_r scales it.
    const double mu0 = 4.0e-7 * M_PI;
    check(near(dc.nu, 1.0 / mu0, 1e-3), "nu = 1/mu_0 for a non-magnetic body");
    const ElementCoefficients magnetic = element_coefficients(make_body(0.0, 1.0, 1000.0), 0.0);
    check(near(magnetic.nu * 1000.0, dc.nu, 1e-3), "mu_r = 1000 divides nu by 1000");

    // The identity the whole conditioning story rests on.
    const double omega = 2.0 * M_PI * 1e6;
    const ElementCoefficients ac = element_coefficients(copper, omega);
    const Complex jw(0.0, omega);
    check(near(ac.alpha, jw * ac.beta, 1e-6 * std::abs(ac.alpha)),
          "alpha == j*omega*beta -- the source of the system's asymmetry");

    // And the literal definitions, so a sign error in either is visible.
    const double eps = 8.8541878128e-12;
    check(near(ac.beta, Complex(5.8e7, omega * eps), 1e-6), "beta = sigma + j*omega*eps");
    check(near(ac.alpha, Complex(-omega * omega * eps, omega * 5.8e7), 1e-3),
          "alpha = j*omega*sigma - omega^2*eps");
}

// The curl-curl block on the unit tet, against values derived by hand.
//
//   curl W_le = 2 grad(La) x grad(Lb) for local edge (a,b), with
//   grad L0 = (-1,-1,-1), grad L1 = (1,0,0), grad L2 = (0,1,0),
//   grad L3 = (0,0,1) and the edge order of kTetLocalEdgeVerts:
//
//     c0 = (0,-2, 2)   c1 = ( 2, 0,-2)   c2 = (-2, 2, 0)
//     c3 = (0, 0, 2)   c4 = ( 2, 0, 0)   c5 = ( 0,-2, 0)
//
//   K[i][j] = nu * V * (ci . cj),  V = 1/6
//
// so every entry is nu/6 times one of 8, -4, 0, 4.
void test_kernel_AA_hand_values() {
    const TetGeometry g = unit_tet();
    check(near(g.volume, 1.0 / 6.0), "the reference tet has volume 1/6");

    ElementCoefficients c;
    c.nu = 1.0;
    c.alpha = Complex(0.0, 0.0);  // curl-curl alone
    c.beta = Complex(0.0, 0.0);

    std::array<Complex, 36> aa{};
    kernel_AA(g, c, aa.data());

    // The dot products c_i . c_j, written out. Not generated from the code
    // under test.
    const int dots[6][6] = {{8, -4, -4, 4, 0, 4},  {-4, 8, -4, -4, 4, 0},
                            {-4, -4, 8, 0, -4, -4}, {4, -4, 0, 4, 0, 0},
                            {0, 4, -4, 0, 4, 0},   {4, 0, -4, 0, 0, 4}};

    bool all_match = true;
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            const Complex want(dots[i][j] / 6.0, 0.0);
            if (!near(aa[static_cast<std::size_t>(i * 6 + j)], want, 1e-14)) all_match = false;
        }
    }
    check(all_match, "the whole 6x6 curl-curl block matches the hand-computed values");
    check(near(aa[0], Complex(4.0 / 3.0, 0.0), 1e-14),
          "K_AA[0][0] = nu*V*|curl W_0|^2 = (1/6)*8 = 4/3");

    // nu multiplies the whole block.
    c.nu = 7.0;
    kernel_AA(g, c, aa.data());
    check(near(aa[0], Complex(7.0 * 4.0 / 3.0, 0.0), 1e-13), "nu scales the block linearly");
}

void test_kernel_AA_symmetry() {
    const TetGeometry g = unit_tet(1.7);
    ElementCoefficients c;
    c.nu = 3.0;
    c.alpha = Complex(2.0, -5.0);

    std::array<Complex, 36> aa{};
    kernel_AA(g, c, aa.data());

    bool symmetric = true;
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) {
            if (!near(aa[static_cast<std::size_t>(i * 6 + j)],
                      aa[static_cast<std::size_t>(j * 6 + i)], 1e-14)) {
                symmetric = false;
            }
        }
    }
    check(symmetric, "the A-A block is symmetric (real basis, symmetric bilinear form)");
}

// curl(grad) = 0, and the discrete gradient of a P1 nodal field lies EXACTLY
// in the Whitney space: the edge coefficient for local edge (a,b) is
// p[b] - p[a]. So the curl-curl block must annihilate every such vector.
//
// This is the sharpest available check on the curl formula and on the edge
// ordering together -- get either wrong and it fails.
void test_kernel_AA_annihilates_gradients() {
    const TetGeometry g = unit_tet(0.4);
    ElementCoefficients c;
    c.nu = 2.5;
    c.alpha = Complex(0.0, 0.0);

    std::array<Complex, 36> aa{};
    kernel_AA(g, c, aa.data());

    const std::array<std::array<double, 4>, 4> nodal = {{{1, 0, 0, 0},
                                                         {0, 1, 0, 0},
                                                         {3, -1, 2, 7},
                                                         {-2.5, 1.25, 0.5, 9}}};
    for (const std::array<double, 4>& p : nodal) {
        std::array<double, 6> a{};
        for (int le = 0; le < 6; ++le) {
            const auto& e = kTetLocalEdgeVerts[static_cast<std::size_t>(le)];
            a[static_cast<std::size_t>(le)] = p[static_cast<std::size_t>(e.second)] -
                                              p[static_cast<std::size_t>(e.first)];
        }
        double worst = 0.0;
        for (int i = 0; i < 6; ++i) {
            Complex row(0.0, 0.0);
            for (int j = 0; j < 6; ++j) {
                row += aa[static_cast<std::size_t>(i * 6 + j)] * a[static_cast<std::size_t>(j)];
            }
            worst = std::max(worst, std::abs(row));
        }
        check(worst < 1e-13,
              "the curl-curl block annihilates the discrete gradient of a P1 field "
              "(worst residual " + std::to_string(worst) + ")");
    }
}

/// Rank by Gaussian elimination with partial pivoting.
int rank_of(const std::array<Complex, 36>& m, double tol) {
    std::array<std::array<double, 6>, 6> a{};
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 6; ++j) a[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] =
            m[static_cast<std::size_t>(i * 6 + j)].real();
    }
    int rank = 0;
    for (int col = 0; col < 6 && rank < 6; ++col) {
        int pivot = -1;
        double best = tol;
        for (int r = rank; r < 6; ++r) {
            if (std::abs(a[static_cast<std::size_t>(r)][static_cast<std::size_t>(col)]) > best) {
                best = std::abs(a[static_cast<std::size_t>(r)][static_cast<std::size_t>(col)]);
                pivot = r;
            }
        }
        if (pivot < 0) continue;
        std::swap(a[static_cast<std::size_t>(rank)], a[static_cast<std::size_t>(pivot)]);
        for (int r = 0; r < 6; ++r) {
            if (r == rank) continue;
            const double f = a[static_cast<std::size_t>(r)][static_cast<std::size_t>(col)] /
                             a[static_cast<std::size_t>(rank)][static_cast<std::size_t>(col)];
            for (int cc = 0; cc < 6; ++cc) {
                a[static_cast<std::size_t>(r)][static_cast<std::size_t>(cc)] -=
                    f * a[static_cast<std::size_t>(rank)][static_cast<std::size_t>(cc)];
            }
        }
        ++rank;
    }
    return rank;
}

void test_kernel_AA_rank() {
    const TetGeometry g = unit_tet();
    ElementCoefficients c;
    c.nu = 1.0;
    c.alpha = Complex(0.0, 0.0);
    std::array<Complex, 36> aa{};
    kernel_AA(g, c, aa.data());

    check(rank_of(aa, 1e-10) == 3,
          "the curl-curl block has rank 3: six edges minus the three-dimensional "
          "gradient subspace");

    // Adding the mass term regularises it -- which is exactly why the
    // low-frequency limit is the hard one.
    c.alpha = Complex(1.0, 0.0);
    kernel_AA(g, c, aa.data());
    check(rank_of(aa, 1e-10) == 6, "adding the mass term makes it full rank");
}

// The mass block against an EXACT integration, done independently here.
//
// This exists because a negative control exposed its absence: replacing the
// degree-2 quadrature with the degree-1 rule left every other test in this
// file green. Symmetry, rank, positive-definiteness and scaling all survive
// a wrong quadrature -- none of them pins a value. Only an exact number does.
//
// W_i . W_j for local edges (a,b) and (c,d) expands to
//
//   L_a L_c (gLb.gLd) - L_a L_d (gLb.gLc) - L_b L_c (gLa.gLd) + L_b L_d (gLa.gLc)
//
// and the barycentric integrals are closed-form: integral L_p L_q dV is
// V/10 when p == q and V/20 otherwise (the simplex formula used in
// tests/test_quadrature.cpp). So the whole 6x6 is computable without any
// quadrature at all.
void test_mass_against_exact_integration() {
    const TetGeometry g = unit_tet(1.4);
    const Complex alpha(2.0, -0.5);

    ElementCoefficients c;
    c.nu = 0.0;  // mass alone
    c.alpha = alpha;
    std::array<Complex, 36> aa{};
    kernel_AA(g, c, aa.data());

    const auto integral_LL = [&](int p, int q) {
        return g.volume * (p == q ? 0.1 : 0.05);
    };

    double worst = 0.0;
    for (int i = 0; i < 6; ++i) {
        const int a = kTetLocalEdgeVerts[static_cast<std::size_t>(i)].first;
        const int b = kTetLocalEdgeVerts[static_cast<std::size_t>(i)].second;
        for (int j = 0; j < 6; ++j) {
            const int cc = kTetLocalEdgeVerts[static_cast<std::size_t>(j)].first;
            const int d = kTetLocalEdgeVerts[static_cast<std::size_t>(j)].second;

            const double exact =
                g.grad_L[static_cast<std::size_t>(b)].dot(g.grad_L[static_cast<std::size_t>(d)]) *
                    integral_LL(a, cc) -
                g.grad_L[static_cast<std::size_t>(b)].dot(g.grad_L[static_cast<std::size_t>(cc)]) *
                    integral_LL(a, d) -
                g.grad_L[static_cast<std::size_t>(a)].dot(g.grad_L[static_cast<std::size_t>(d)]) *
                    integral_LL(b, cc) +
                g.grad_L[static_cast<std::size_t>(a)].dot(g.grad_L[static_cast<std::size_t>(cc)]) *
                    integral_LL(b, d);

            worst = std::max(worst,
                             std::abs(aa[static_cast<std::size_t>(i * 6 + j)] - alpha * exact));
        }
    }
    check(worst < 1e-14,
          "every mass entry matches the closed-form barycentric integral "
          "(worst error " + std::to_string(worst) + ")");
}

void test_mass_is_positive_definite() {
    const TetGeometry g = unit_tet(0.9);
    ElementCoefficients c;
    c.nu = 0.0;  // mass alone
    c.alpha = Complex(1.0, 0.0);

    std::array<Complex, 36> aa{};
    kernel_AA(g, c, aa.data());

    for (int i = 0; i < 6; ++i) {
        check(aa[static_cast<std::size_t>(i * 6 + i)].real() > 0.0,
              "every mass diagonal entry is positive");
    }

    const std::array<std::array<double, 6>, 3> probes = {
        {{1, 1, 1, 1, 1, 1}, {1, -2, 3, -4, 5, -6}, {0.3, 0, -1.7, 2.2, 0, 0.1}}};
    for (const std::array<double, 6>& x : probes) {
        double q = 0.0;
        for (int i = 0; i < 6; ++i) {
            for (int j = 0; j < 6; ++j) {
                q += x[static_cast<std::size_t>(i)] *
                     aa[static_cast<std::size_t>(i * 6 + j)].real() * x[static_cast<std::size_t>(j)];
            }
        }
        check(q > 0.0, "x^T M x > 0 for a non-zero x -- the mass block is positive definite");
    }
}

// Under a uniform scaling of the tet by s, with coefficients held fixed.
// The Whitney function has units of 1/length and its curl of 1/length^2, so
//
//   curl-curl:  (1/s^2)^2 * s^3 = 1/s
//   mass:       (1/s)^2   * s^3 = s
//
// They scale in OPPOSITE directions, which is worth knowing: refining a mesh
// makes the curl-curl term grow and the mass term shrink, and that ratio is
// the low-frequency conditioning problem in miniature.
void test_scaling() {
    ElementCoefficients curl_only;
    curl_only.nu = 1.0;
    curl_only.alpha = Complex(0.0, 0.0);

    std::array<Complex, 36> a1{}, a2{};
    kernel_AA(unit_tet(1.0), curl_only, a1.data());
    kernel_AA(unit_tet(2.0), curl_only, a2.data());
    check(near(a2[0].real(), a1[0].real() / 2.0, 1e-13),
          "the curl-curl block scales as 1/s: doubling the tet halves it");

    ElementCoefficients mass_only;
    mass_only.nu = 0.0;
    mass_only.alpha = Complex(1.0, 0.0);
    kernel_AA(unit_tet(1.0), mass_only, a1.data());
    kernel_AA(unit_tet(2.0), mass_only, a2.data());
    check(near(a2[0].real(), a1[0].real() * 2.0, 1e-13),
          "the mass block scales as s: doubling the tet doubles it");
}

// The curl-curl term takes a shortcut -- nu * V * (curl.curl), no quadrature
// loop -- because its integrand is constant. This checks the shortcut
// against the quadrature path it replaces, which is two routes to one
// answer rather than a restatement.
void test_curl_shortcut_matches_quadrature() {
    const TetGeometry g = unit_tet(1.3);
    ElementCoefficients c;
    c.nu = 2.0;
    c.alpha = Complex(0.0, 0.0);

    std::array<Complex, 36> shortcut{};
    kernel_AA(g, c, shortcut.data());

    const QuadratureRule q = tet_rule_degree_2();
    std::array<double, 36> integrated{};
    for (int k = 0; k < q.count; ++k) {
        const double w = q.weights[static_cast<std::size_t>(k)] * g.volume;
        for (int i = 0; i < 6; ++i) {
            for (int j = 0; j < 6; ++j) {
                integrated[static_cast<std::size_t>(i * 6 + j)] +=
                    w * c.nu * whitney_edge_curl(g, i).dot(whitney_edge_curl(g, j));
            }
        }
    }

    double worst = 0.0;
    for (int i = 0; i < 36; ++i) {
        worst = std::max(worst, std::abs(shortcut[static_cast<std::size_t>(i)].real() -
                                         integrated[static_cast<std::size_t>(i)]));
    }
    check(worst < 1e-14,
          "the constant-integrand shortcut agrees with the quadrature path to round-off");
}

// The two terms must add, not replace one another.
void test_terms_add() {
    const TetGeometry g = unit_tet(0.6);
    std::array<Complex, 36> curl{}, mass{}, both{};

    ElementCoefficients c;
    c.nu = 1.5;
    c.alpha = Complex(0.0, 0.0);
    kernel_AA(g, c, curl.data());

    c.nu = 0.0;
    c.alpha = Complex(4.0, -3.0);
    kernel_AA(g, c, mass.data());

    c.nu = 1.5;
    kernel_AA(g, c, both.data());

    bool adds = true;
    for (int i = 0; i < 36; ++i) {
        if (!near(both[static_cast<std::size_t>(i)],
                  curl[static_cast<std::size_t>(i)] + mass[static_cast<std::size_t>(i)], 1e-13)) {
            adds = false;
        }
    }
    check(adds, "the curl-curl and mass terms add into one block");

    // And at DC the block is the curl-curl term alone.
    const ElementCoefficients dc = element_coefficients(make_body(5.8e7), 0.0);
    std::array<Complex, 36> at_dc{};
    kernel_AA(g, dc, at_dc.data());
    bool real_only = true;
    for (int i = 0; i < 36; ++i) {
        if (at_dc[static_cast<std::size_t>(i)].imag() != 0.0) real_only = false;
    }
    check(real_only, "at DC the A-A block is purely real -- no mass term at all");
}

}  // namespace

int main() {
    test_coefficients();
    test_kernel_AA_hand_values();
    test_kernel_AA_symmetry();
    test_kernel_AA_annihilates_gradients();
    test_kernel_AA_rank();
    test_mass_against_exact_integration();
    test_mass_is_positive_definite();
    test_scaling();
    test_curl_shortcut_matches_quadrature();
    test_terms_add();

    std::cout << (g_checks - g_failures) << "/" << g_checks << " checks passed\n";
    return g_failures == 0 ? 0 : 1;
}
