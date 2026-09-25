// Tests for the tetrahedral quadrature rules.
//
// Step 1 of `docs/ASSEMBLY_PLAN.md` Sec. 8, done before any physics uses
// them. A quadrature rule is worth testing on its own because its failure
// mode is silent: a rule that is one degree short still returns a plausible
// number, and the resulting element matrix is wrong by a few percent with
// nothing to distinguish it from a modelling error.
//
// The check is against EXACT integrals, not against a finer rule. For
// monomials in barycentric coordinates there is a closed form:
//
//   integral over the reference tet of L0^p L1^q L2^r L3^s
//     = 6V * p! q! r! s! / (p + q + r + s + 3)!
//
// (the standard simplex formula, with 3! = 6 for a tetrahedron). Testing a
// rule against a finer rule would only show the two agree, which they would
// even if both were wrong.

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "aphi_solver/quadrature.hpp"

using namespace aphi_solver;

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

double factorial(int n) {
    double f = 1.0;
    for (int i = 2; i <= n; ++i) f *= i;
    return f;
}

/// The exact integral of L0^e0 L1^e1 L2^e2 L3^e3 over a tet of volume V.
double exact_monomial(const std::array<int, 4>& e, double volume) {
    const int total = e[0] + e[1] + e[2] + e[3];
    double numerator = 1.0;
    for (int i = 0; i < 4; ++i) numerator *= factorial(e[i]);
    return 6.0 * volume * numerator / factorial(total + 3);
}

/// The same integral, by quadrature.
double quadrature_monomial(const QuadratureRule& r, const std::array<int, 4>& e, double volume) {
    double sum = 0.0;
    for (int k = 0; k < r.count; ++k) {
        double f = 1.0;
        for (int i = 0; i < 4; ++i) {
            f *= std::pow(r.points[static_cast<std::size_t>(k)][static_cast<std::size_t>(i)],
                          e[static_cast<std::size_t>(i)]);
        }
        sum += r.weights[static_cast<std::size_t>(k)] * f;
    }
    return volume * sum;
}

/// Every monomial of exactly total degree `d`, as exponent tuples.
std::vector<std::array<int, 4>> monomials_of_degree(int d) {
    std::vector<std::array<int, 4>> out;
    for (int a = 0; a <= d; ++a) {
        for (int b = 0; b + a <= d; ++b) {
            for (int c = 0; c + b + a <= d; ++c) {
                const int e = d - a - b - c;
                out.push_back({a, b, c, e});
            }
        }
    }
    return out;
}

// ---------------------------------------------------------------------------

void test_rule_shape() {
    for (const QuadratureRule& r : {tet_rule_degree_1(), tet_rule_degree_2()}) {
        double sum = 0.0;
        for (int k = 0; k < r.count; ++k) sum += r.weights[static_cast<std::size_t>(k)];
        check(std::abs(sum - 1.0) < 1e-15,
              "weights sum to 1, so a constant integrand gives exactly the volume");

        for (int k = 0; k < r.count; ++k) {
            double bary = 0.0;
            bool inside = true;
            for (int i = 0; i < 4; ++i) {
                const double L = r.points[static_cast<std::size_t>(k)][static_cast<std::size_t>(i)];
                bary += L;
                if (L < 0.0 || L > 1.0) inside = false;
            }
            check(std::abs(bary - 1.0) < 1e-15, "each point's barycentric coordinates sum to 1");
            check(inside, "each point lies inside the tet");
        }
    }
}

// The claim in the header is "exact for degree <= d". That is what gets
// tested -- every monomial up to that degree, against the closed form.
void test_exactness(const QuadratureRule& r, const std::string& name) {
    const double volume = 1.0 / 6.0;  // the reference tet

    for (int d = 0; d <= r.degree; ++d) {
        double worst = 0.0;
        for (const std::array<int, 4>& e : monomials_of_degree(d)) {
            const double got = quadrature_monomial(r, e, volume);
            const double want = exact_monomial(e, volume);
            worst = std::max(worst, std::abs(got - want));
        }
        check(worst < 1e-15,
              name + " is exact for every monomial of degree " + std::to_string(d) +
                  " (worst error " + std::to_string(worst) + ")");
    }
}

// And the converse: a rule must NOT be exact one degree past its claim, or
// the claim understates it and the degree field is wrong. This is what makes
// the exactness test above meaningful rather than vacuous -- it shows the
// test can distinguish degrees at all.
void test_not_exact_beyond_claim() {
    const double volume = 1.0 / 6.0;
    const QuadratureRule r1 = tet_rule_degree_1();

    double worst = 0.0;
    for (const std::array<int, 4>& e : monomials_of_degree(2)) {
        worst = std::max(worst, std::abs(quadrature_monomial(r1, e, volume) -
                                         exact_monomial(e, volume)));
    }
    check(worst > 1e-3,
          "the 1-point rule is visibly WRONG at degree 2 -- so the degree-2 rule's "
          "exactness is a real property, not an artefact of the test");
}

// The element matrices need degree 2 and no more (ASSEMBLY_PLAN Sec. 2).
// Spot-check the three shapes that actually occur, in barycentric form.
void test_the_integrands_that_matter() {
    const QuadratureRule r = tet_rule_degree_2();
    const double volume = 1.0 / 6.0;

    // W . W and grad(S) . grad(S) reduce to sums of terms Li*Lj; the
    // hardest cases are Li^2 and Li*Lj.
    const double li2 = quadrature_monomial(r, {2, 0, 0, 0}, volume);
    check(std::abs(li2 - exact_monomial({2, 0, 0, 0}, volume)) < 1e-16,
          "L0^2 exactly: V/10");
    check(std::abs(li2 - volume / 10.0) < 1e-16, "and V/10 is what the closed form gives");

    const double lilj = quadrature_monomial(r, {1, 1, 0, 0}, volume);
    check(std::abs(lilj - exact_monomial({1, 1, 0, 0}, volume)) < 1e-16, "L0*L1 exactly");
    check(std::abs(lilj - volume / 20.0) < 1e-16, "which is V/20");

    // A constant integrand -- the curl-curl case -- must give exactly V,
    // since curl W is constant per tet and the whole block is V * (a dot
    // product). This is the shortcut ASSEMBLY_PLAN Sec. 2 takes.
    const double one = quadrature_monomial(r, {0, 0, 0, 0}, volume);
    check(std::abs(one - volume) < 1e-17, "a constant integrand integrates to exactly V");
}

}  // namespace

int main() {
    test_rule_shape();
    test_exactness(tet_rule_degree_1(), "the 1-point rule");
    test_exactness(tet_rule_degree_2(), "the 4-point rule");
    test_not_exact_beyond_claim();
    test_the_integrands_that_matter();

    std::cout << (g_checks - g_failures) << "/" << g_checks << " checks passed\n";
    return g_failures == 0 ? 0 : 1;
}
