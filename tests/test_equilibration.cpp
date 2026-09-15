// Minimal, dependency-free test runner -- same approach as test_conditioning.cpp.

#include <cmath>
#include <complex>
#include <iostream>
#include <string>
#include <vector>

#include "aphi_solver/complex_matrix.hpp"
#include "aphi_solver/conditioning.hpp"
#include "aphi_solver/equilibration.hpp"

using aphi_solver::Complex;
using aphi_solver::ComplexMatrix;

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

bool near(Complex a, Complex b, double tol = 1e-6) { return std::abs(a - b) <= tol; }
bool near(double a, double b, double tol = 1e-6) { return std::abs(a - b) <= tol; }

bool matrices_near(const ComplexMatrix& a, const ComplexMatrix& b, double tol = 1e-6) {
    if (a.rows() != b.rows() || a.cols() != b.cols()) return false;
    for (int i = 0; i < a.rows(); ++i)
        for (int j = 0; j < a.cols(); ++j)
            if (!near(a(i, j), b(i, j), tol)) return false;
    return true;
}

// A deliberately badly-scaled 3x3 complex symmetric matrix: diagonal entries span
// nine orders of magnitude, which is the kind of thing low-frequency A-Phi
// assembly is prone to (curl-curl terms vs. tiny mass terms at small omega).
ComplexMatrix make_badly_scaled_symmetric_matrix() {
    ComplexMatrix A(3, 3);
    A(0, 0) = Complex(1.0e6, 0.0);
    A(1, 1) = Complex(1.0e-2, -5.0);
    A(2, 2) = Complex(3.0, 1.0);

    Complex a01(50.0, 20.0);
    Complex a02(0.4, -0.1);
    Complex a12(2.0, 0.5);
    A(0, 1) = a01; A(1, 0) = a01;
    A(0, 2) = a02; A(2, 0) = a02;
    A(1, 2) = a12; A(2, 1) = a12;
    return A;
}

void test_equilibration_preserves_symmetry() {
    ComplexMatrix A = make_badly_scaled_symmetric_matrix();
    std::vector<double> d = aphi_solver::compute_symmetric_equilibration(A, 10);
    ComplexMatrix scaled = aphi_solver::apply_symmetric_equilibration(A, d);

    bool symmetric = true;
    for (int i = 0; i < scaled.rows() && symmetric; ++i)
        for (int j = 0; j < scaled.cols() && symmetric; ++j)
            if (!near(scaled(i, j), scaled(j, i))) symmetric = false;

    check(symmetric, "symmetric equilibration keeps a complex-symmetric matrix symmetric");
}

void test_equilibration_balances_magnitudes() {
    ComplexMatrix A = make_badly_scaled_symmetric_matrix();
    std::vector<double> d = aphi_solver::compute_symmetric_equilibration(A, 15);
    ComplexMatrix scaled = aphi_solver::apply_symmetric_equilibration(A, d);

    // After equilibration, every row's largest-magnitude entry should be close to 1
    // -- that's the whole point of the transform.
    bool balanced = true;
    for (int i = 0; i < scaled.rows(); ++i) {
        double row_max = 0.0;
        for (int j = 0; j < scaled.cols(); ++j) row_max = std::max(row_max, std::abs(scaled(i, j)));
        if (std::abs(row_max - 1.0) > 0.2) balanced = false;  // loose tolerance: a few Ruiz iterations, not full convergence
    }
    check(balanced, "equilibration brings row-max magnitudes close to 1");
}

void test_equilibration_improves_condition_number() {
    ComplexMatrix A = make_badly_scaled_symmetric_matrix();
    double kappa_before = aphi_solver::estimate_condition_number(A, 200);

    std::vector<double> d = aphi_solver::compute_symmetric_equilibration(A, 15);
    ComplexMatrix scaled = aphi_solver::apply_symmetric_equilibration(A, d);
    double kappa_after = aphi_solver::estimate_condition_number(scaled, 200);

    check(kappa_after < kappa_before,
          "equilibration reduces the estimated condition number (before=" +
              std::to_string(kappa_before) + ", after=" + std::to_string(kappa_after) + ")");
}

void test_equilibration_zero_row_is_left_unscaled() {
    ComplexMatrix A(2, 2);
    A(0, 0) = Complex(0.0, 0.0);
    A(0, 1) = Complex(0.0, 0.0);
    A(1, 0) = Complex(0.0, 0.0);
    A(1, 1) = Complex(5.0, 0.0);

    std::vector<double> d;
    bool threw = false;
    try {
        d = aphi_solver::compute_symmetric_equilibration(A, 5);
    } catch (...) {
        threw = true;
    }
    check(!threw, "a zero row does not throw or produce NaN/inf");
    if (!threw) {
        check(std::isfinite(d[0]) && std::isfinite(d[1]), "equilibration factors stay finite with a zero row");
    }
}

void test_round_trip_solve_invariance() {
    ComplexMatrix A = make_badly_scaled_symmetric_matrix();
    ComplexMatrix b(3, 1);
    b(0, 0) = Complex(1.0, 0.5);
    b(1, 0) = Complex(-0.3, 0.2);
    b(2, 0) = Complex(0.7, -0.1);

    ComplexMatrix x_direct = aphi_solver::solve_dense(A, b);

    std::vector<double> d = aphi_solver::compute_symmetric_equilibration(A, 15);
    ComplexMatrix A_scaled = aphi_solver::apply_symmetric_equilibration(A, d);
    ComplexMatrix b_scaled = aphi_solver::scale_rhs(b, d);
    ComplexMatrix y = aphi_solver::solve_dense(A_scaled, b_scaled);
    ComplexMatrix x_recovered = aphi_solver::recover_equilibrated_solution(y, d);

    check(matrices_near(x_direct, x_recovered, 1e-6),
          "equilibrate -> solve -> recover reproduces the direct solve");
}

// Composability check: physics-motivated conditioning (Phase 04, conditioning.hpp)
// and numerical equilibration (this file) are meant to stack. Apply the
// scaled-scalar-potential transform first, then equilibrate the result, solve,
// and undo both layers in reverse order -- the physical (a, Phi) should come back
// unchanged, same as in test_conditioning.cpp's round-trip tests.
void test_composes_with_physics_conditioning() {
    using aphi_solver::APhiBlockSystem;
    const double omega = 2.0 * M_PI * 60.0;
    const Complex j_omega(0.0, omega);

    APhiBlockSystem sys;
    sys.K_AA = ComplexMatrix(2, 2);
    sys.K_AA(0, 0) = Complex(5.0, 0.0);
    sys.K_AA(0, 1) = Complex(0.2, 0.0);
    sys.K_AA(1, 0) = Complex(0.2, 0.0);
    sys.K_AA(1, 1) = Complex(4.0, 0.0);
    sys.K_APhi = ComplexMatrix(2, 1);
    sys.K_APhi(0, 0) = Complex(1.0, 0.0);
    sys.K_APhi(1, 0) = Complex(0.5, 0.0);
    sys.K_PhiA = ComplexMatrix(1, 2);
    sys.K_PhiA(0, 0) = j_omega * Complex(1.0, 0.0);
    sys.K_PhiA(0, 1) = j_omega * Complex(0.5, 0.0);
    sys.K_PhiPhi = ComplexMatrix(1, 1);
    sys.K_PhiPhi(0, 0) = j_omega * Complex(3.0, 0.0);
    sys.rhs_A = ComplexMatrix(2, 1);
    sys.rhs_A(0, 0) = Complex(1.0, 0.5);
    sys.rhs_A(1, 0) = Complex(0.3, -0.2);
    sys.rhs_Phi = ComplexMatrix(1, 1);
    sys.rhs_Phi(0, 0) = Complex(0.7, 0.1);

    auto stack = [](const ComplexMatrix& top, const ComplexMatrix& bottom) {
        ComplexMatrix out(top.rows() + bottom.rows(), top.cols());
        out.set_block(0, 0, top);
        out.set_block(top.rows(), 0, bottom);
        return out;
    };

    // Reference: solve the untouched natural system.
    ComplexMatrix x_natural = aphi_solver::solve_dense(aphi_solver::assemble_dense(sys), stack(sys.rhs_A, sys.rhs_Phi));

    // Layer 1: physics conditioning (scaled scalar potential).
    APhiBlockSystem transformed = aphi_solver::apply_scaled_scalar_potential(sys, omega);
    ComplexMatrix A1 = aphi_solver::assemble_dense(transformed);
    ComplexMatrix b1 = stack(transformed.rhs_A, transformed.rhs_Phi);

    // Layer 2: numerical equilibration on top of the physics-transformed system.
    std::vector<double> d = aphi_solver::compute_symmetric_equilibration(A1, 15);
    ComplexMatrix A2 = aphi_solver::apply_symmetric_equilibration(A1, d);
    ComplexMatrix b2 = aphi_solver::scale_rhs(b1, d);

    ComplexMatrix y = aphi_solver::solve_dense(A2, b2);
    ComplexMatrix x_transformed = aphi_solver::recover_equilibrated_solution(y, d);  // undo equilibration

    ComplexMatrix a_recovered = x_transformed.block(0, 0, sys.num_A(), 1);
    ComplexMatrix phi_prime_recovered = x_transformed.block(sys.num_A(), 0, sys.num_Phi(), 1);
    ComplexMatrix phi_recovered = aphi_solver::recover_scaled_scalar_potential(phi_prime_recovered, omega);  // undo physics layer

    ComplexMatrix a_natural = x_natural.block(0, 0, sys.num_A(), 1);
    ComplexMatrix phi_natural = x_natural.block(sys.num_A(), 0, sys.num_Phi(), 1);

    check(matrices_near(a_recovered, a_natural, 1e-6),
          "stacking equilibration on top of physics conditioning still recovers the correct A");
    check(matrices_near(phi_recovered, phi_natural, 1e-6),
          "stacking equilibration on top of physics conditioning still recovers the correct Phi");
}

}  // namespace

int main() {
    test_equilibration_preserves_symmetry();
    test_equilibration_balances_magnitudes();
    test_equilibration_improves_condition_number();
    test_equilibration_zero_row_is_left_unscaled();
    test_round_trip_solve_invariance();
    test_composes_with_physics_conditioning();

    std::cout << (g_checks - g_failures) << "/" << g_checks << " checks passed\n";
    return g_failures == 0 ? 0 : 1;
}
