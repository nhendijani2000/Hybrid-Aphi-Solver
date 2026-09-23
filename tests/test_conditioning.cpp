// Minimal, dependency-free test runner. No external test framework is fetched
// here on purpose -- see the comment in tests/CMakeLists.txt for why.

#include <cmath>
#include <complex>
#include <iostream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include "aphi_solver/complex_matrix.hpp"
#include "aphi_solver/conditioning.hpp"
#include "aphi_solver/gauge_variants.hpp"

using aphi_solver::Complex;
using aphi_solver::ComplexMatrix;
using aphi_solver::APhiBlockSystem;
using aphi_solver::ConditioningStrategy;
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

bool near(Complex a, Complex b, double tol = 1e-8) { return std::abs(a - b) <= tol; }
bool near(double a, double b, double tol = 1e-6) { return std::abs(a - b) <= tol; }

bool matrices_near(const ComplexMatrix& a, const ComplexMatrix& b, double tol = 1e-6) {
    if (a.rows() != b.rows() || a.cols() != b.cols()) return false;
    for (int i = 0; i < a.rows(); ++i)
        for (int j = 0; j < a.cols(); ++j)
            if (!near(a(i, j), b(i, j), tol)) return false;
    return true;
}

// Sparse blocks are compared through the dense bridge: at these sizes that
// is the clearest thing to read, and it checks the values rather than the
// storage pattern (two sparse matrices can hold the same matrix with
// different explicit-zero sets).
bool matrices_near(const SparseMatrixZ& a, const SparseMatrixZ& b, double tol = 1e-6) {
    return matrices_near(aphi_solver::to_dense_matrix(a), aphi_solver::to_dense_matrix(b), tol);
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

// A small (n_A = 2, n_Phi = 1) test system. K_PhiA and K_PhiPhi are deliberately
// built as j*omega times a "clean" value, so symmetric row scaling should recover
// exactly K_APhi^T and a real K_PhiPhi -- this is the property the transform is
// meant to produce (see docs/CONDITIONING.md).
APhiBlockSystem make_test_system(double omega) {
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

    return sys;
}

void test_symmetric_row_scaling_arithmetic() {
    const double omega = 2.0 * M_PI * 60.0;
    APhiBlockSystem sys = make_test_system(omega);
    APhiBlockSystem scaled = aphi_solver::apply_symmetric_row_scaling(sys, omega);

    // K_PhiA_scaled should equal K_APhi^T exactly, by construction of the test
    // system. The plain transpose, not the conjugate one: the A-Phi system is
    // complex symmetric, not Hermitian (docs/CONDITIONING.md).
    check(matrices_near(scaled.K_PhiA, sys.K_APhi.transposed()),
          "symmetric row scaling recovers K_APhi^T in K_PhiA");
    check(near(scaled.K_PhiPhi.at(0, 0), Complex(3.0, 0.0)),
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

    ComplexMatrix x_natural =
        aphi_solver::solve_dense(aphi_solver::assemble_dense(sys), aphi_solver::assemble_dense_rhs(sys));

    APhiBlockSystem scaled = aphi_solver::apply_symmetric_row_scaling(sys, omega);
    ComplexMatrix x_scaled =
        aphi_solver::solve_dense(aphi_solver::assemble_dense(scaled), aphi_solver::assemble_dense_rhs(scaled));

    check(matrices_near(x_natural, x_scaled, 1e-6),
          "dividing the Phi row by j*omega does not change the physical solution");
}

void test_scaled_scalar_potential_round_trip() {
    const double omega = 2.0 * M_PI * 60.0;
    APhiBlockSystem sys = make_test_system(omega);

    ComplexMatrix x_natural =
        aphi_solver::solve_dense(aphi_solver::assemble_dense(sys), aphi_solver::assemble_dense_rhs(sys));
    const std::vector<Complex> phi_natural =
        aphi_solver::from_column(x_natural.block(sys.num_A(), 0, sys.num_Phi(), 1));

    APhiBlockSystem transformed = aphi_solver::apply_scaled_scalar_potential(sys, omega);
    ComplexMatrix x_transformed = aphi_solver::solve_dense(aphi_solver::assemble_dense(transformed),
                                                            aphi_solver::assemble_dense_rhs(transformed));
    const std::vector<Complex> phi_prime =
        aphi_solver::from_column(x_transformed.block(sys.num_A(), 0, sys.num_Phi(), 1));

    const std::vector<Complex> phi_recovered = aphi_solver::recover_scaled_scalar_potential(phi_prime, omega);

    check(vectors_near(phi_recovered, phi_natural, 1e-6),
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

// assemble_sparse is the production path (what a sparse direct solver is
// handed); assemble_dense is the small-system ground truth. They must agree
// entry for entry, including getting every block's row/column offset right --
// an off-by-one in one of the four block placements is the obvious way this
// goes wrong, and it would otherwise only show up as a wrong field much
// later.
void test_sparse_and_dense_assembly_agree() {
    const double omega = 2.0 * M_PI * 60.0;
    APhiBlockSystem sys = make_test_system(omega);

    const SparseMatrixZ sparse_full = aphi_solver::assemble_sparse(sys);
    const ComplexMatrix dense_full = aphi_solver::assemble_dense(sys);

    check(sparse_full.rows() == sys.num_A() + sys.num_Phi() && sparse_full.cols() == sparse_full.rows(),
          "assemble_sparse has the combined [a; Phi] shape");
    check(matrices_near(aphi_solver::to_dense_matrix(sparse_full), dense_full),
          "assemble_sparse and assemble_dense agree entry for entry");

    // Spot-check the off-diagonal block offsets directly, so a symmetric
    // mistake in both assembly routines could not pass the comparison above.
    check(near(sparse_full.at(0, sys.num_A()), Complex(1.0, 0.0)),
          "assemble_sparse places K_APhi at (0, n_A)");
    check(near(sparse_full.at(sys.num_A(), 0), Complex(0.0, omega)),
          "assemble_sparse places K_PhiA at (n_A, 0)");
}

// Phase 03.5 step 6: the dense (complex_matrix.hpp) and sparse
// (gauge_variants.hpp) estimate_condition_number overloads are one contract,
// not two estimators that happen to share a name. Before the Sept 2026
// unification the dense one ran a FIXED 100 iterations with no convergence
// test while the sparse one ran up to 500 with one, so the same matrix could
// get two different answers depending only on which type it was stored in --
// and nothing in the suite would have noticed. This pins them together.
void test_condition_number_overloads_agree() {
    // Real, symmetric, well-separated spectrum so both estimators have a
    // single right answer to find: diag(1000, 100, 10, 1), kappa = 1000.
    const std::vector<double> diagonal = {1000.0, 100.0, 10.0, 1.0};
    const int n = static_cast<int>(diagonal.size());

    aphi_solver::SparseMatrix sparse(n, n);
    ComplexMatrix dense(n, n);
    for (int i = 0; i < n; ++i) {
        sparse.add(i, i, diagonal[i]);
        dense(i, i) = Complex(diagonal[i], 0.0);
    }
    sparse.compress();

    const double kappa_sparse = aphi_solver::estimate_condition_number(sparse);
    const double kappa_dense = aphi_solver::estimate_condition_number(dense);

    check(near(kappa_dense, 1000.0, 1e-6), "dense overload finds kappa = 1000 at its default");
    check(near(kappa_sparse, 1000.0, 1e-6), "sparse overload finds kappa = 1000 at its default");
    check(std::abs(kappa_dense - kappa_sparse) < 1e-6 * 1000.0,
          "both overloads agree on the same matrix at their defaults");

    // Pinning the DEFAULT iteration budget needs a matrix whose answer
    // actually depends on it. The spectrum above converges in far under 100
    // iterations, so a check written on it passes whether the default is 100
    // or 500 -- verified by negative control (Sept 2026: the default was
    // temporarily reverted to 100 and the whole suite still passed). The
    // matrix below is built so power iteration converges SLOWLY but steadily:
    // A^H A has eigenvalues (1, 0.95, 0.01), so the Rayleigh quotient's error
    // falls like 0.95^(2k) -- still ~3e-5 at k = 100, reaching the 1e-12
    // tolerance near k = 270. A 500-iteration budget therefore converges and a
    // 100-iteration one does not, and the two give measurably different kappa.
    const double sigma_slow = std::sqrt(0.95);
    ComplexMatrix slow(3, 3);
    slow(0, 0) = Complex(1.0, 0.0);
    slow(1, 1) = Complex(sigma_slow, 0.0);
    slow(2, 2) = Complex(0.1, 0.0);

    const double kappa_default = aphi_solver::estimate_condition_number(slow);
    const double kappa_500 = aphi_solver::estimate_condition_number(slow, 500, 1e-12);
    const double kappa_100 = aphi_solver::estimate_condition_number(slow, 100, 1e-12);

    check(near(kappa_default, 10.0, 1e-6), "slow-converging matrix has kappa = 10");
    check(near(kappa_default, kappa_500, 1e-12), "dense default budget is 500");
    check(!near(kappa_default, kappa_100, 1e-9),
          "the 100-iteration budget is measurably worse (so the check above can fail)");

    // Shared degenerate-input contract: 1.0 for n <= 1, throw if not square.
    check(near(aphi_solver::estimate_condition_number(ComplexMatrix(1, 1)), 1.0, 0.0),
          "dense overload returns 1.0 for n = 1");
    aphi_solver::SparseMatrix one(1, 1);
    one.compress();
    check(near(aphi_solver::estimate_condition_number(one), 1.0, 0.0),
          "sparse overload returns 1.0 for n = 1");

    bool dense_threw = false;
    try {
        aphi_solver::estimate_condition_number(ComplexMatrix(2, 3));
    } catch (const std::invalid_argument&) {
        dense_threw = true;
    }
    check(dense_threw, "dense overload throws on a non-square matrix");

    bool sparse_threw = false;
    try {
        aphi_solver::SparseMatrix oblong(2, 3);
        oblong.compress();
        aphi_solver::estimate_condition_number(oblong);
    } catch (const std::invalid_argument&) {
        sparse_threw = true;
    }
    check(sparse_threw, "sparse overload throws on a non-square matrix");
}

}  // namespace

int main() {
    test_symmetric_row_scaling_arithmetic();
    test_symmetric_row_scaling_zero_omega_throws();
    test_symmetric_row_scaling_preserves_solution();
    test_scaled_scalar_potential_round_trip();
    test_recommend_strategy();
    test_condition_number_known_diagonal();
    test_sparse_and_dense_assembly_agree();
    test_condition_number_overloads_agree();

    std::cout << (g_checks - g_failures) << "/" << g_checks << " checks passed\n";
    return g_failures == 0 ? 0 : 1;
}
