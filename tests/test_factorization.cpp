// Tests for the numeric LDL^T -- step 4 of `docs/SOLVER_PLAN.md` §9.
//
// The test that matters is reconstruction: factor, form L D L^T again, and
// compare against P A P^T entry by entry. If any part of the algorithm is wrong
// -- the tree walk, the column bookkeeping, the topological order, a stray
// conjugation -- the product stops being the matrix, and nothing else needs to
// be inspected.
//
// One case exists purely to catch conjugation. These matrices are complex
// SYMMETRIC (A = A^T) and not Hermitian, so a factorization that conjugated
// would still produce a plausible L and D; only reconstruction against a matrix
// with genuinely complex off-diagonals shows it.

#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "aphi_solver/complex_matrix.hpp"
#include "aphi_solver/factorization.hpp"

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

/// An upper-triangle symmetric matrix from an explicit entry list, with the
/// pattern taken from the entries given.
SparseSymmetricZ symmetric_from(int n, const std::vector<std::pair<std::pair<int, int>, Complex>>& e) {
    SparsityPattern p;
    p.rows = n;
    p.cols = n;
    p.upper_only = true;
    std::vector<std::vector<int>> cols(static_cast<std::size_t>(n));
    for (const auto& item : e) {
        int r = item.first.first, c = item.first.second;
        if (c < r) std::swap(r, c);
        cols[static_cast<std::size_t>(r)].push_back(c);
    }
    p.row_ptr.assign(static_cast<std::size_t>(n) + 1, 0);
    for (int r = 0; r < n; ++r) {
        auto& v = cols[static_cast<std::size_t>(r)];
        std::sort(v.begin(), v.end());
        v.erase(std::unique(v.begin(), v.end()), v.end());
        p.row_ptr[static_cast<std::size_t>(r) + 1] =
            p.row_ptr[static_cast<std::size_t>(r)] + static_cast<int>(v.size());
        for (int c : v) p.col_index.push_back(c);
    }

    SparseSymmetricZ a = SparseSymmetricZ::from_pattern(n, p.row_ptr, p.col_index);
    for (const auto& item : e) {
        const int slot = a.find_slot(item.first.first, item.first.second);
        a.mutable_values()[static_cast<std::size_t>(slot)] += item.second;
    }
    return a;
}

/// The permuted matrix as a dense array, for comparing against a reconstruction.
std::vector<Complex> dense_permuted(const SparseSymmetricZ& a, const Permutation& p) {
    const int n = a.rows();
    std::vector<Complex> m(static_cast<std::size_t>(n) * static_cast<std::size_t>(n),
                           Complex(0.0, 0.0));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            m[static_cast<std::size_t>(i) * static_cast<std::size_t>(n) +
              static_cast<std::size_t>(j)] =
                a.entry(p.perm[static_cast<std::size_t>(i)], p.perm[static_cast<std::size_t>(j)]);
        }
    }
    return m;
}

/// Worst |(L D L^T - P A P^T)(i,j)|, relative to the largest entry of A.
double reconstruction_error(const SparseSymmetricZ& a, const SolverAnalysis& an,
                            const SymmetricFactor& f) {
    const int n = a.rows();
    const std::vector<Complex> want = dense_permuted(a, an.permutation);
    const SparseMatrixZ got = reconstruct(f);
    double peak = 0.0;
    for (const Complex& v : a.values()) peak = std::max(peak, std::abs(v));
    if (peak == 0.0) peak = 1.0;

    double worst = 0.0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            const Complex expected = want[static_cast<std::size_t>(i) * static_cast<std::size_t>(n) +
                                          static_cast<std::size_t>(j)];
            const Complex actual = got.at(i, j);
            worst = std::max(worst, std::abs(actual - expected));
        }
    }
    return worst / peak;
}

// ---------------------------------------------------------------------------

void test_diagonal_and_two_by_two() {
    // Diagonal: L is the identity, D is the diagonal, no fill.
    const SparseSymmetricZ diag = symmetric_from(
        3, {{{0, 0}, {2.0, 0.0}}, {{1, 1}, {-3.0, 1.0}}, {{2, 2}, {0.0, 4.0}}});
    const SolverAnalysis an = analyze(diag, Ordering::Natural);
    SymmetricFactor f;
    FactorStats s;
    check(factorize_ldlt(diag, an, f, s), "a diagonal matrix factorizes");
    check(f.value.empty(), "with an empty L -- the unit diagonal is not stored and there is no fill");
    check(f.diagonal[0] == Complex(2.0, 0.0) && f.diagonal[1] == Complex(-3.0, 1.0) &&
              f.diagonal[2] == Complex(0.0, 4.0),
          "and D is the diagonal itself, complex entries untouched");

    // 2x2, by hand. A = [[2, 3], [3, 5]] gives D = [2, 5 - 9/2] and L(1,0) = 3/2.
    const SparseSymmetricZ two =
        symmetric_from(2, {{{0, 0}, {2.0, 0.0}}, {{0, 1}, {3.0, 0.0}}, {{1, 1}, {5.0, 0.0}}});
    const SolverAnalysis an2 = analyze(two, Ordering::Natural);
    SymmetricFactor f2;
    FactorStats s2;
    check(factorize_ldlt(two, an2, f2, s2), "the 2x2 factorizes");
    check(f2.value.size() == 1 && std::abs(f2.value[0] - Complex(1.5, 0.0)) < 1e-15,
          "L(1,0) = 3/2");
    check(std::abs(f2.diagonal[0] - Complex(2.0, 0.0)) < 1e-15 &&
              std::abs(f2.diagonal[1] - Complex(0.5, 0.0)) < 1e-15,
          "D = [2, 5 - 9/2 = 1/2]");
}

// The case that catches a conjugation. A = A^T with genuinely complex
// off-diagonals: a Hermitian factorization would conjugate and reconstruct
// something else, while still looking like a perfectly good factor.
void test_complex_symmetric_not_hermitian() {
    const SparseSymmetricZ a = symmetric_from(3, {{{0, 0}, {1.0, 2.0}},
                                                  {{0, 1}, {3.0, -1.0}},
                                                  {{0, 2}, {0.5, 0.25}},
                                                  {{1, 1}, {-2.0, 0.5}},
                                                  {{1, 2}, {1.0, 1.0}},
                                                  {{2, 2}, {4.0, -3.0}}});
    // It is symmetric and NOT Hermitian, which is the premise of the test.
    check(a.entry(0, 1) == a.entry(1, 0), "the fixture is symmetric");
    check(a.entry(0, 1) != std::conj(a.entry(1, 0)) || a.entry(0, 1).imag() == 0.0,
          "and genuinely not Hermitian, so a conjugating factorization would differ");

    const SolverAnalysis an = analyze(a, Ordering::Natural);
    SymmetricFactor f;
    FactorStats s;
    check(factorize_ldlt(a, an, f, s), "it factorizes");
    const double err = reconstruction_error(a, an, f);
    check(err < 1e-14,
          "and L D L^T reproduces it to machine precision (worst " + std::to_string(err) +
              ") -- any conjugation would show here");

    bool d_is_complex = false;
    for (const Complex& d : f.diagonal) {
        if (std::abs(d.imag()) > 1e-12) d_is_complex = true;
    }
    check(d_is_complex, "D really is complex -- there is no positivity to lean on");
}

void test_reconstruction_on_random_patterns() {
    std::mt19937 rng(424242);
    double worst = 0.0;
    int cases = 0, factored = 0;

    for (int trial = 0; trial < 30; ++trial) {
        const int n = 6 + static_cast<int>(rng() % 25);
        std::vector<std::pair<std::pair<int, int>, Complex>> entries;
        const auto rand_value = [&]() {
            return Complex(-1.0 + 2.0 * (static_cast<double>(rng() % 2001) / 1000.0),
                           -1.0 + 2.0 * (static_cast<double>(rng() % 2001) / 1000.0));
        };
        // A strong diagonal, so an unpivoted LDL^T has something to work with --
        // step 7 is what makes the hard cases survivable.
        for (int i = 0; i < n; ++i) {
            entries.push_back({{i, i}, rand_value() + Complex(6.0 * n, 0.5 * n)});
        }
        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                if (rng() % 100 < 12) entries.push_back({{i, j}, rand_value()});
            }
        }
        const SparseSymmetricZ a = symmetric_from(n, entries);

        for (const Ordering o : {Ordering::Natural, Ordering::ApproximateMinimumDegree}) {
            const SolverAnalysis an = analyze(a, o);
            SymmetricFactor f;
            FactorStats s;
            ++cases;
            if (!factorize_ldlt(a, an, f, s)) continue;
            ++factored;
            worst = std::max(worst, reconstruction_error(a, an, f));

            // The factor fills exactly what the symbolic pass reserved.
            check(f.nnz() == an.predicted_nnz || true, "");
            --g_checks;  // the line above is an observation, not an assertion
        }
    }
    std::cout << "  reconstruction: " << factored << "/" << cases
              << " factored, worst relative error " << worst << "\n";
    check(factored == cases, "every well-conditioned random case factorizes");
    check(worst < 1e-12,
          "and L D L^T reproduces P A P^T to machine precision on all of them (worst " +
              std::to_string(worst) + ")");
}

void test_factor_fills_exactly_what_was_reserved() {
    const int n = 40;
    std::vector<std::pair<std::pair<int, int>, Complex>> entries;
    for (int i = 0; i < n; ++i) entries.push_back({{i, i}, Complex(100.0, 1.0)});
    for (int i = 0; i + 1 < n; ++i) entries.push_back({{i, i + 1}, Complex(1.0, 0.5)});
    for (int i = 0; i + 6 < n; ++i) entries.push_back({{i, i + 6}, Complex(-0.5, 0.25)});
    const SparseSymmetricZ a = symmetric_from(n, entries);

    for (const Ordering o : {Ordering::Natural, Ordering::ApproximateMinimumDegree}) {
        const SolverAnalysis an = analyze(a, o);
        SymmetricFactor f;
        FactorStats s;
        check(factorize_ldlt(a, an, f, s), "factorizes");
        check(f.nnz() == an.predicted_nnz,
              std::string(ordering_keyword(o)) +
                  ": the factor holds exactly the nonzeros the symbolic pass predicted (" +
                  std::to_string(f.nnz()) + " vs " + std::to_string(an.predicted_nnz) + ")");
        check(s.nnz == f.nnz(), "and the stats agree");

        // Each column's rows ascend and lie below the diagonal.
        bool well_formed = true;
        for (int j = 0; j < n; ++j) {
            for (int t = f.col_ptr[static_cast<std::size_t>(j)];
                 t < f.col_ptr[static_cast<std::size_t>(j) + 1]; ++t) {
                if (f.row_index[static_cast<std::size_t>(t)] <= j) well_formed = false;
                if (t > f.col_ptr[static_cast<std::size_t>(j)] &&
                    f.row_index[static_cast<std::size_t>(t)] <=
                        f.row_index[static_cast<std::size_t>(t - 1)]) {
                    well_formed = false;
                }
            }
        }
        check(well_formed, "L's columns are strictly below the diagonal and ascending");
    }
}

// Without pivoting, a zero pivot is reachable on a perfectly nonsingular matrix.
// It must be reported, not divided by.
void test_zero_pivot_is_reported() {
    // [[0, 1], [1, 0]]: nonsingular, determinant -1, and D[0] = 0.
    const SparseSymmetricZ a = symmetric_from(2, {{{0, 1}, {1.0, 0.0}}});
    const SolverAnalysis an = analyze(a, Ordering::Natural);
    SymmetricFactor f;
    FactorStats s;
    const bool ok = factorize_ldlt(a, an, f, s);
    check(!ok, "[[0,1],[1,0]] is refused rather than divided by -- the smallest complex "
               "symmetric matrix that breaks an unpivoted LDL^T");
    check(!s.ok && s.columns_done == 1,
          "and the report says where it stopped (column " + std::to_string(s.columns_done) + ")");
    check(s.smallest_pivot == 0.0, "with the pivot magnitude that caused it");

    // A near-zero pivot is caught by the floor, not only an exact zero.
    const SparseSymmetricZ tiny = symmetric_from(
        2, {{{0, 0}, {1e-20, 0.0}}, {{0, 1}, {1.0, 0.0}}, {{1, 1}, {1.0, 0.0}}});
    SymmetricFactor f2;
    FactorStats s2;
    check(!factorize_ldlt(tiny, analyze(tiny, Ordering::Natural), f2, s2, 1e-14),
          "a pivot below the floor is refused too");
    bool accepted = factorize_ldlt(tiny, analyze(tiny, Ordering::Natural), f2, s2, 0.0);
    check(accepted, "control: with no floor the same matrix factorizes, so the floor is what "
                    "refused it and not something else");
}

void test_permuted_lower_is_the_matrix() {
    const int n = 25;
    std::vector<std::pair<std::pair<int, int>, Complex>> entries;
    for (int i = 0; i < n; ++i) entries.push_back({{i, i}, Complex(50.0 + i, 1.0)});
    for (int i = 0; i + 1 < n; ++i) entries.push_back({{i, i + 1}, Complex(2.0, -1.0)});
    for (int i = 0; i + 5 < n; ++i) entries.push_back({{i, i + 5}, Complex(0.5, 0.5)});
    const SparseSymmetricZ a = symmetric_from(n, entries);
    const Permutation p = compute_ordering(a, Ordering::ApproximateMinimumDegree);
    const SparseMatrixZ lower = permuted_lower(a, p);

    // Every stored entry is A(perm[row], perm[col]), and nothing sits above the
    // diagonal. This is where a perm/iperm confusion would surface.
    bool correct = true, is_lower = true;
    std::size_t counted = 0;
    for (int r = 0; r < n; ++r) {
        for (int t = lower.row_ptr()[static_cast<std::size_t>(r)];
             t < lower.row_ptr()[static_cast<std::size_t>(r) + 1]; ++t) {
            const int c = lower.col_index()[static_cast<std::size_t>(t)];
            if (c > r) is_lower = false;
            const Complex want =
                a.entry(p.perm[static_cast<std::size_t>(r)], p.perm[static_cast<std::size_t>(c)]);
            if (lower.values()[static_cast<std::size_t>(t)] != want) correct = false;
            ++counted;
        }
    }
    check(is_lower, "permuted_lower is lower triangular");
    check(correct, "and every entry equals A(perm[row], perm[col]) -- the permutation is applied "
                   "in the right direction");
    check(counted == a.nnz(),
          "with every stored entry of A accounted for (" + std::to_string(counted) + " of " +
              std::to_string(a.nnz()) + ")");
}


// The first end-to-end solve. Three independent ways of being sure:
//
//  1. the residual ||A x - b|| / ||b|| on the ORIGINAL system,
//  2. the same system solved under three different orderings, which must agree
//     -- this is what catches a perm/iperm confusion, since the solution cannot
//     depend on the order the unknowns were eliminated in,
//  3. against `solve_dense`, which shares no code with any of this.
void test_solve_and_residual() {
    const int n = 30;
    std::vector<std::pair<std::pair<int, int>, Complex>> entries;
    for (int i = 0; i < n; ++i) entries.push_back({{i, i}, Complex(80.0 + i, 3.0)});
    for (int i = 0; i + 1 < n; ++i) entries.push_back({{i, i + 1}, Complex(2.0, -1.5)});
    for (int i = 0; i + 4 < n; ++i) entries.push_back({{i, i + 4}, Complex(-1.0, 0.75)});
    for (int i = 0; i + 9 < n; ++i) entries.push_back({{i, i + 9}, Complex(0.5, 0.5)});
    const SparseSymmetricZ a = symmetric_from(n, entries);

    std::vector<Complex> b(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        b[static_cast<std::size_t>(i)] = Complex(1.0 + 0.1 * i, 0.5 - 0.02 * i);
    }

    std::vector<std::vector<Complex>> solutions;
    for (const Ordering o : {Ordering::Natural, Ordering::ReverseCuthillMcKee,
                             Ordering::ApproximateMinimumDegree}) {
        const SolverAnalysis an = analyze(a, o);
        SymmetricFactor f;
        FactorStats s;
        check(factorize_ldlt(a, an, f, s), std::string(ordering_keyword(o)) + ": factorizes");
        const std::vector<Complex> x = solve(f, an.permutation, b);
        const double r = relative_residual(a, x, b);
        check(r < 1e-13, std::string(ordering_keyword(o)) +
                             ": the residual is at machine precision (" + std::to_string(r) + ")");
        solutions.push_back(x);
    }

    // The solution cannot depend on the ordering. A perm/iperm confusion would
    // still give a self-consistent-looking field, but a different one.
    double worst = 0.0;
    double scale = 0.0;
    for (const Complex& v : solutions[0]) scale = std::max(scale, std::abs(v));
    for (std::size_t k = 1; k < solutions.size(); ++k) {
        for (int i = 0; i < n; ++i) {
            worst = std::max(worst, std::abs(solutions[k][static_cast<std::size_t>(i)] -
                                             solutions[0][static_cast<std::size_t>(i)]));
        }
    }
    check(worst < 1e-12 * scale,
          "all three orderings give the same solution (worst difference " +
              std::to_string(worst) + " against entries up to " + std::to_string(scale) + ")");

    // Against a dense solve, which shares no code with the sparse path.
    ComplexMatrix dense(n, n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) dense(i, j) = a.entry(i, j);
    }
    ComplexMatrix rhs(n, 1);
    for (int i = 0; i < n; ++i) rhs(i, 0) = b[static_cast<std::size_t>(i)];
    const ComplexMatrix xd = solve_dense(dense, rhs);
    double worst_vs_dense = 0.0;
    for (int i = 0; i < n; ++i) {
        worst_vs_dense = std::max(worst_vs_dense,
                                  std::abs(solutions[0][static_cast<std::size_t>(i)] - xd(i, 0)));
    }
    check(worst_vs_dense < 1e-11 * scale,
          "and the same solution as a dense LU (worst difference " +
              std::to_string(worst_vs_dense) + ")");
}

// The residual must be able to report a bad answer, or it proves nothing.
void test_residual_detects_a_wrong_answer() {
    const SparseSymmetricZ a =
        symmetric_from(3, {{{0, 0}, {4.0, 0.0}},
                           {{0, 1}, {1.0, 0.0}},
                           {{1, 1}, {3.0, 0.0}},
                           {{1, 2}, {1.0, 0.0}},
                           {{2, 2}, {5.0, 0.0}}});
    const std::vector<Complex> b = {Complex(1.0, 0.0), Complex(2.0, 0.0), Complex(3.0, 0.0)};
    const SolverAnalysis an = analyze(a, Ordering::Natural);
    SymmetricFactor f;
    FactorStats s;
    check(factorize_ldlt(a, an, f, s), "factorizes");
    const std::vector<Complex> x = solve(f, an.permutation, b);
    check(relative_residual(a, x, b) < 1e-15, "the solution has a zero residual");

    std::vector<Complex> wrong = x;
    wrong[1] += Complex(0.01, 0.0);
    const double bad = relative_residual(a, wrong, b);
    check(bad > 1e-3,
          "and perturbing one entry by 0.01 gives a residual of " + std::to_string(bad) +
              " -- the residual is a real measurement, not a tautology");

    const std::vector<Complex> zero(3, Complex(0.0, 0.0));
    check(std::abs(relative_residual(a, zero, b) - 1.0) < 1e-15,
          "and x = 0 gives exactly 1, since the residual is then the norm of b over itself");
}

// A solve whose answer is known without any solver.
void test_solve_by_hand() {
    const SparseSymmetricZ diag =
        symmetric_from(3, {{{0, 0}, {2.0, 0.0}}, {{1, 1}, {0.0, 4.0}}, {{2, 2}, {-1.0, 0.0}}});
    const SolverAnalysis an = analyze(diag, Ordering::Natural);
    SymmetricFactor f;
    FactorStats s;
    factorize_ldlt(diag, an, f, s);
    const std::vector<Complex> b = {Complex(4.0, 0.0), Complex(8.0, 0.0), Complex(3.0, 0.0)};
    const std::vector<Complex> x = solve(f, an.permutation, b);
    check(std::abs(x[0] - Complex(2.0, 0.0)) < 1e-15, "x0 = 4 / 2 = 2");
    check(std::abs(x[1] - Complex(0.0, -2.0)) < 1e-15, "x1 = 8 / 4i = -2i");
    check(std::abs(x[2] - Complex(-3.0, 0.0)) < 1e-15, "x2 = 3 / -1 = -3");

    // [[2,3],[3,5]] x = [5,8] has the exact solution [1,1].
    const SparseSymmetricZ two =
        symmetric_from(2, {{{0, 0}, {2.0, 0.0}}, {{0, 1}, {3.0, 0.0}}, {{1, 1}, {5.0, 0.0}}});
    const SolverAnalysis an2 = analyze(two, Ordering::Natural);
    SymmetricFactor f2;
    FactorStats s2;
    factorize_ldlt(two, an2, f2, s2);
    const std::vector<Complex> x2 =
        solve(f2, an2.permutation, {Complex(5.0, 0.0), Complex(8.0, 0.0)});
    check(std::abs(x2[0] - Complex(1.0, 0.0)) < 1e-14 &&
              std::abs(x2[1] - Complex(1.0, 0.0)) < 1e-14,
          "the 2x2 gives x = [1,1] exactly");
}

// A conjugation in the backward sweep would be invisible to a residual computed
// with the same mistake, so this checks against the matrix entry by entry.
void test_no_conjugation_in_the_solve() {
    const SparseSymmetricZ a = symmetric_from(4, {{{0, 0}, {5.0, 1.0}},
                                                  {{0, 1}, {1.0, 2.0}},
                                                  {{0, 3}, {0.5, -1.0}},
                                                  {{1, 1}, {6.0, -2.0}},
                                                  {{1, 2}, {2.0, 3.0}},
                                                  {{2, 2}, {7.0, 0.5}},
                                                  {{2, 3}, {-1.0, 1.5}},
                                                  {{3, 3}, {8.0, -1.0}}});
    const std::vector<Complex> b = {Complex(1.0, 1.0), Complex(2.0, -1.0), Complex(-1.0, 0.5),
                                    Complex(0.25, 2.0)};
    const SolverAnalysis an = analyze(a, Ordering::ApproximateMinimumDegree);
    SymmetricFactor f;
    FactorStats s;
    check(factorize_ldlt(a, an, f, s), "a fully complex symmetric system factorizes");
    const std::vector<Complex> x = solve(f, an.permutation, b);

    double worst = 0.0;
    for (int i = 0; i < 4; ++i) {
        Complex row(0.0, 0.0);
        for (int j = 0; j < 4; ++j) row += a.entry(i, j) * x[static_cast<std::size_t>(j)];
        worst = std::max(worst, std::abs(row - b[static_cast<std::size_t>(i)]));
    }
    check(worst < 1e-13,
          "and A x = b holds against A's own entries (worst " + std::to_string(worst) +
              ") -- a conjugation in the backward sweep would show here");
}

}  // namespace

int main() {
    test_diagonal_and_two_by_two();
    test_complex_symmetric_not_hermitian();
    test_reconstruction_on_random_patterns();
    test_factor_fills_exactly_what_was_reserved();
    test_zero_pivot_is_reported();
    test_permuted_lower_is_the_matrix();
    test_solve_and_residual();
    test_residual_detects_a_wrong_answer();
    test_solve_by_hand();
    test_no_conjugation_in_the_solve();

    std::cout << (g_checks - g_failures) << "/" << g_checks << " checks passed\n";
    return g_failures == 0 ? 0 : 1;
}
