// Minimal, dependency-free test runner. No external test framework is fetched
// here on purpose -- see the comment in tests/CMakeLists.txt for why.

#include <cmath>
#include <complex>
#include <iostream>
#include <string>

#include "aphi_solver/complex_matrix.hpp"
#include "aphi_solver/conditioning.hpp"

using aphi_solver::Complex;
using aphi_solver::ComplexMatrix;
using aphi_solver::APhiBlockSystem;
using aphi_solver::ConditioningStrategy;

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

bool near(Complex a, Complex b, double tol = 1e-8) { return std::abs(a - b) <= tol; }
bool near(double a, double b, double tol = 1e-6) { return std::abs(a - b) <= tol; }

bool matrices_near(const ComplexMatrix& a, const ComplexMatrix& b, double tol = 1e-6) {
    if (a.rows() != b.rows() || a.cols() != b.cols()) return false;
    for (int i = 0; i < a.rows(); ++i)
        for (int j = 0; j < a.cols(); ++j)
            if (!near(a(i, j), b(i, j), tol)) return false;
    return true;
}

// A small (n_A = 2, n_Phi = 1) test system. K_PhiA and K_PhiPhi are deliberately
// built as j*omega times a "clean" value, so symmetric row scaling should recover
// exactly K_APhi^T and a real K_PhiPhi -- this is the property the transform is
// meant to produce (see docs/CONDITIONING.md).
APhiBlockSystem make_test_system(double omega) {
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

    return sys;
}

ComplexMatrix stack(const ComplexMatrix& top, const ComplexMatrix& bottom) {
    ComplexMatrix out(top.rows() + bottom.rows(), top.cols());
    out.set_block(0, 0, top);
    out.set_block(top.rows(), 0, bottom);
    return out;
}

void test_symmetric_row_scaling_arithmetic() {
    const double omega = 2.0 * M_PI * 60.0;
    APhiBlockSystem sys = make_test_system(omega);
    APhiBlockSystem scaled = aphi_solver::apply_symmetric_row_scaling(sys, omega);

    // K_PhiA_scaled should equal K_APhi^T exactly, by construction of the test system.
    check(matrices_near(scaled.K_PhiA, sys.K_APhi.hermitian()),
          "symmetric row scaling recovers K_APhi^T in K_PhiA");
    check(near(scaled.K_PhiPhi(0, 0), Complex(3.0, 0.0)),
          "symmetric row scaling recovers a real K_PhiPhi");
    check(matrices_near(scaled.K_AA, sys.K_AA), "symmetric row scaling leaves K_AA untouched");
}

void test_symmetric_row_scaling_zero_omega_throws() {
    APhiBlockSystem sys = make_test_system(1.0);
    bool threw = false;
    try {
        aphi_solver::apply_symmetric_row_scaling(sys, 0.0);
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    check(threw, "symmetric row scaling rejects omega == 0");
}

void test_symmetric_row_scaling_preserves_solution() {
    const double omega = 2.0 * M_PI * 60.0;
    APhiBlockSystem sys = make_test_system(omega);

    ComplexMatrix rhs_natural = stack(sys.rhs_A, sys.rhs_Phi);
    ComplexMatrix x_natural = aphi_solver::solve_dense(aphi_solver::assemble_dense(sys), rhs_natural);

    APhiBlockSystem scaled = aphi_solver::apply_symmetric_row_scaling(sys, omega);
    ComplexMatrix rhs_scaled = stack(scaled.rhs_A, scaled.rhs_Phi);
    ComplexMatrix x_scaled = aphi_solver::solve_dense(aphi_solver::assemble_dense(scaled), rhs_scaled);

    check(matrices_near(x_natural, x_scaled, 1e-6),
          "dividing the Phi row by j*omega does not change the physical solution");
}

void test_scaled_scalar_potential_round_trip() {
    const double omega = 2.0 * M_PI * 60.0;
    APhiBlockSystem sys = make_test_system(omega);

    ComplexMatrix rhs_natural = stack(sys.rhs_A, sys.rhs_Phi);
    ComplexMatrix x_natural = aphi_solver::solve_dense(aphi_solver::assemble_dense(sys), rhs_natural);
    ComplexMatrix phi_natural = x_natural.block(sys.num_A(), 0, sys.num_Phi(), 1);

    APhiBlockSystem transformed = aphi_solver::apply_scaled_scalar_potential(sys, omega);
    ComplexMatrix rhs_transformed = stack(transformed.rhs_A, transformed.rhs_Phi);
    ComplexMatrix x_transformed = aphi_solver::solve_dense(aphi_solver::assemble_dense(transformed), rhs_transformed);
    ComplexMatrix phi_prime = x_transformed.block(sys.num_A(), 0, sys.num_Phi(), 1);

    ComplexMatrix phi_recovered = aphi_solver::recover_scaled_scalar_potential(phi_prime, omega);

    check(matrices_near(phi_recovered, phi_natural, 1e-6),
          "Phi = j*omega*Phi' recovers the same physical Phi as the natural system");

    ComplexMatrix a_natural = x_natural.block(0, 0, sys.num_A(), 1);
    ComplexMatrix a_transformed = x_transformed.block(0, 0, sys.num_A(), 1);
    check(matrices_near(a_natural, a_transformed, 1e-6),
          "the scaled-scalar-potential transform leaves the recovered A unchanged");
}

void test_recommend_strategy() {
    check(aphi_solver::recommend_strategy(10.0, 1000.0) == ConditioningStrategy::SymmetricRowScaled,
          "recommend_strategy picks SymmetricRowScaled below the caller's crossover");
    check(aphi_solver::recommend_strategy(10000.0, 1000.0) == ConditioningStrategy::ScaledScalarPotential,
          "recommend_strategy picks ScaledScalarPotential at/above the caller's crossover");
}

void test_condition_number_known_diagonal() {
    // A diagonal matrix's condition number is just the ratio of its largest to
    // smallest diagonal magnitude -- a direct check of estimate_condition_number
    // against a hand-known answer, independent of the A-Phi system above.
    ComplexMatrix D(3, 3);
    D(0, 0) = Complex(100.0, 0.0);
    D(1, 1) = Complex(10.0, 0.0);
    D(2, 2) = Complex(1.0, 0.0);
    double kappa = aphi_solver::estimate_condition_number(D, 200);
    check(near(kappa, 100.0, 1e-3), "condition number of diag(100,10,1) is ~100");

    ComplexMatrix I = ComplexMatrix::identity(4);
    double kappa_identity = aphi_solver::estimate_condition_number(I, 200);
    check(near(kappa_identity, 1.0, 1e-6), "condition number of the identity is 1");
}

}  // namespace

int main() {
    test_symmetric_row_scaling_arithmetic();
    test_symmetric_row_scaling_zero_omega_throws();
    test_symmetric_row_scaling_preserves_solution();
    test_scaled_scalar_potential_round_trip();
    test_recommend_strategy();
    test_condition_number_known_diagonal();

    std::cout << (g_checks - g_failures) << "/" << g_checks << " checks passed\n";
    return g_failures == 0 ? 0 : 1;
}
