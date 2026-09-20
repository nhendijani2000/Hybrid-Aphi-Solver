// Minimal, dependency-free test runner. No external test framework is fetched
// here on purpose -- see the comment in tests/CMakeLists.txt for why.
//
// Diagonal equilibration (equilibration.hpp): the purely numerical scaling
// applied to the assembled system before factorization, distinct from the
// physics-motivated transforms in conditioning.hpp. Both layers are exercised
// here, including the check that they compose. Operates on the sparse system
// type since Phase 03.5 (docs/ROADMAP.md); the dense ComplexMatrix appears
// only where the small-system ground-truth solve needs it.

#include <cmath>
#include <complex>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>

#include "aphi_solver/complex_matrix.hpp"
#include "aphi_solver/conditioning.hpp"
#include "aphi_solver/equilibration.hpp"

using aphi_solver::Complex;
using aphi_solver::ComplexMatrix;
using aphi_solver::SparseMatrixZ;

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

bool vectors_near(const std::vector<Complex>& a, const std::vector<Complex>& b, double tol = 1e-6) {
    if (a.size() != b.size()) return false;
    for (std::size_t i = 0; i < a.size(); ++i)
        if (!near(a[i], b[i], tol)) return false;
    return true;
}

SparseMatrixZ make_block(int rows, int cols, const std::vector<std::tuple<int, int, Complex>>& entries) {
    SparseMatrixZ m(rows, cols);
    for (const auto& [r, c, v] : entries) m.add(r, c, v);
    m.compress();
    return m;
}

// A deliberately badly-scaled 3x3 complex symmetric matrix: diagonal entries span
// nine orders of magnitude, which is the kind of thing low-frequency A-Phi
// assembly is prone to (curl-curl terms vs. tiny mass terms at small omega).
SparseMatrixZ make_badly_scaled_symmetric_matrix() {
    const Complex a01(50.0, 20.0);
    const Complex a02(0.4, -0.1);
    const Complex a12(2.0, 0.5);
    return make_block(3, 3,
                       {{0, 0, Complex(1.0e6, 0.0)},
                        {1, 1, Complex(1.0e-2, -5.0)},
                        {2, 2, Complex(3.0, 1.0)},
                        {0, 1, a01},
                        {1, 0, a01},
                        {0, 2, a02},
                        {2, 0, a02},
                        {1, 2, a12},
                        {2, 1, a12}});
}

void test_equilibration_preserves_symmetry() {
    const SparseMatrixZ A = make_badly_scaled_symmetric_matrix();
    const std::vector<double> d = aphi_solver::compute_symmetric_equilibration(A, 10);
    const SparseMatrixZ scaled = aphi_solver::apply_symmetric_equilibration(A, d);

    bool symmetric = true;
    for (int i = 0; i < scaled.rows() && symmetric; ++i)
        for (int j = 0; j < scaled.cols() && symmetric; ++j)
            if (!near(scaled.at(i, j), scaled.at(j, i))) symmetric = false;

    check(symmetric, "symmetric equilibration keeps a complex-symmetric matrix symmetric");
    check(scaled.nnz() == A.nnz(), "equilibration leaves the sparsity pattern unchanged");
}

void test_equilibration_balances_magnitudes() {
    const SparseMatrixZ A = make_badly_scaled_symmetric_matrix();
    const std::vector<double> d = aphi_solver::compute_symmetric_equilibration(A, 15);
    const SparseMatrixZ scaled = aphi_solver::apply_symmetric_equilibration(A, d);

    // After equilibration, every row's largest-magnitude entry should be close to 1
    // -- that's the whole point of the transform.
    bool balanced = true;
    for (int i = 0; i < scaled.rows(); ++i) {
        double row_max = 0.0;
        for (int j = 0; j < scaled.cols(); ++j) row_max = std::max(row_max, std::abs(scaled.at(i, j)));
        if (std::abs(row_max - 1.0) > 0.2) balanced = false;  // loose tolerance: a few Ruiz iterations, not full convergence
    }
    check(balanced, "equilibration brings row-max magnitudes close to 1");
}

void test_equilibration_improves_condition_number() {
    const SparseMatrixZ A = make_badly_scaled_symmetric_matrix();
    const double kappa_before = aphi_solver::estimate_condition_number(aphi_solver::to_dense_matrix(A), 200);

    const std::vector<double> d = aphi_solver::compute_symmetric_equilibration(A, 15);
    const SparseMatrixZ scaled = aphi_solver::apply_symmetric_equilibration(A, d);
    const double kappa_after = aphi_solver::estimate_condition_number(aphi_solver::to_dense_matrix(scaled), 200);

    check(kappa_after < kappa_before,
          "equilibration reduces the estimated condition number (before=" +
              std::to_string(kappa_before) + ", after=" + std::to_string(kappa_after) + ")");
}

void test_equilibration_zero_row_is_left_unscaled() {
    // Row 0 is structurally empty -- in CSR that is an empty row segment the
    // row-max loop never enters, which is a different code path from a row of
    // explicitly-stored zeros and the one a real assembled system produces for
    // a DOF that has received no contributions.
    const SparseMatrixZ A = make_block(2, 2, {{1, 1, Complex(5.0, 0.0)}});

    std::vector<double> d;
    bool threw = false;
    try {
        d = aphi_solver::compute_symmetric_equilibration(A, 5);
    } catch (...) {
        threw = true;
    }
    check(!threw, "a structurally empty row does not throw or produce NaN/inf");
    if (!threw) {
        check(std::isfinite(d[0]) && std::isfinite(d[1]), "equilibration factors stay finite with an empty row");
        check(near(d[0], 1.0), "an empty row is left unscaled (d == 1)");
    }

    // The same must hold for a row of explicitly-stored zeros.
    const SparseMatrixZ Z = make_block(2, 2, {{0, 0, Complex(0.0, 0.0)}, {1, 1, Complex(5.0, 0.0)}});
    const std::vector<double> dz = aphi_solver::compute_symmetric_equilibration(Z, 5);
    check(std::isfinite(dz[0]) && near(dz[0], 1.0), "an explicitly-zero row is also left unscaled");
}

void test_round_trip_solve_invariance() {
    const SparseMatrixZ A = make_badly_scaled_symmetric_matrix();
    const std::vector<Complex> b{Complex(1.0, 0.5), Complex(-0.3, 0.2), Complex(0.7, -0.1)};

    const ComplexMatrix x_direct =
        aphi_solver::solve_dense(aphi_solver::to_dense_matrix(A), aphi_solver::to_column(b));

    const std::vector<double> d = aphi_solver::compute_symmetric_equilibration(A, 15);
    const SparseMatrixZ A_scaled = aphi_solver::apply_symmetric_equilibration(A, d);
    const std::vector<Complex> b_scaled = aphi_solver::scale_rhs(b, d);
    const ComplexMatrix y = aphi_solver::solve_dense(aphi_solver::to_dense_matrix(A_scaled),
                                                      aphi_solver::to_column(b_scaled));
    const std::vector<Complex> x_recovered =
        aphi_solver::recover_equilibrated_solution(aphi_solver::from_column(y), d);

    check(vectors_near(aphi_solver::from_column(x_direct), x_recovered, 1e-6),
          "equilibrate -> solve -> recover reproduces the direct solve");
}

// Composability check: physics-motivated conditioning (conditioning.hpp) and
// numerical equilibration (this file) are meant to stack. Apply the
// scaled-scalar-potential transform first, then equilibrate the result, solve,
// and undo both layers in reverse order -- the physical (a, Phi) should come back
// unchanged, same as in test_conditioning.cpp's round-trip tests.
void test_composes_with_physics_conditioning() {
    using aphi_solver::APhiBlockSystem;
    const double omega = 2.0 * M_PI * 60.0;
    const Complex j_omega(0.0, omega);

    APhiBlockSystem sys;
    sys.K_AA = make_block(2, 2,
                           {{0, 0, Complex(5.0, 0.0)},
                            {0, 1, Complex(0.2, 0.0)},
                            {1, 0, Complex(0.2, 0.0)},
                            {1, 1, Complex(4.0, 0.0)}});
    sys.K_APhi = make_block(2, 1, {{0, 0, Complex(1.0, 0.0)}, {1, 0, Complex(0.5, 0.0)}});
    sys.K_PhiA = make_block(1, 2, {{0, 0, j_omega * Complex(1.0, 0.0)}, {0, 1, j_omega * Complex(0.5, 0.0)}});
    sys.K_PhiPhi = make_block(1, 1, {{0, 0, j_omega * Complex(3.0, 0.0)}});
    sys.rhs_A = {Complex(1.0, 0.5), Complex(0.3, -0.2)};
    sys.rhs_Phi = {Complex(0.7, 0.1)};

    // Reference: solve the untouched natural system.
    const ComplexMatrix x_natural =
        aphi_solver::solve_dense(aphi_solver::assemble_dense(sys), aphi_solver::assemble_dense_rhs(sys));

    // Layer 1: physics conditioning (scaled scalar potential).
    const APhiBlockSystem transformed = aphi_solver::apply_scaled_scalar_potential(sys, omega);
    const SparseMatrixZ A1 = aphi_solver::assemble_sparse(transformed);
    const std::vector<Complex> b1 = aphi_solver::from_column(aphi_solver::assemble_dense_rhs(transformed));

    // Layer 2: numerical equilibration on top of the physics-transformed system.
    const std::vector<double> d = aphi_solver::compute_symmetric_equilibration(A1, 15);
    const SparseMatrixZ A2 = aphi_solver::apply_symmetric_equilibration(A1, d);
    const std::vector<Complex> b2 = aphi_solver::scale_rhs(b1, d);

    const ComplexMatrix y =
        aphi_solver::solve_dense(aphi_solver::to_dense_matrix(A2), aphi_solver::to_column(b2));
    const std::vector<Complex> x_transformed =
        aphi_solver::recover_equilibrated_solution(aphi_solver::from_column(y), d);  // undo equilibration

    const std::vector<Complex> a_recovered(x_transformed.begin(), x_transformed.begin() + sys.num_A());
    const std::vector<Complex> phi_prime_recovered(x_transformed.begin() + sys.num_A(), x_transformed.end());
    const std::vector<Complex> phi_recovered =
        aphi_solver::recover_scaled_scalar_potential(phi_prime_recovered, omega);  // undo physics layer

    const std::vector<Complex> x_nat = aphi_solver::from_column(x_natural);
    const std::vector<Complex> a_natural(x_nat.begin(), x_nat.begin() + sys.num_A());
    const std::vector<Complex> phi_natural(x_nat.begin() + sys.num_A(), x_nat.end());

    check(vectors_near(a_recovered, a_natural, 1e-6),
          "stacking equilibration on top of physics conditioning still recovers the correct A");
    check(vectors_near(phi_recovered, phi_natural, 1e-6),
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
