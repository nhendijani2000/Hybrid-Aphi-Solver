#pragma once

#include <array>

namespace aphi_solver {

/// Tetrahedral quadrature, in barycentric coordinates.
///
/// A rule is a list of points `L = (L0, L1, L2, L3)` with `sum Li = 1`, and
/// weights that sum to 1. An integral over a tet of volume `V` is then
///
/// ```
///   integral f dV  =  V * sum_k  w[k] * f(L[k])
/// ```
///
/// Barycentric rather than Cartesian because every basis function in this
/// project is written in barycentric coordinates already
/// (`basis_functions.hpp`), and because the rule is then independent of the
/// tet -- the same points serve every element, with only `V` changing.

/// A degree-`d` rule: exact for every polynomial of total degree <= d.
struct QuadratureRule {
    int degree = 0;
    int count = 0;
    std::array<std::array<double, 4>, 4> points{};  // at most 4 points so far
    std::array<double, 4> weights{};
};

/// The 1-point centroid rule. Exact for degree 1.
///
/// Enough for the curl-curl term on its own -- `curl W = 2 grad(Li) x
/// grad(Lj)` is constant per tet -- but for nothing else here.
QuadratureRule tet_rule_degree_1();

/// The classical 4-point rule, exact for degree 2 (and, as it happens, for
/// degree 3 as well, though nothing here needs that).
///
/// Points are the permutations of `(a, b, b, b)` with
/// `a = 0.585410196624968515...`, `b = 0.138196601125010515...`, each
/// weighted 1/4. Those are `(5 + 3*sqrt(5))/20` and `(5 - sqrt(5))/20`.
///
/// **This is exact for every integrand in the element matrices**, not an
/// approximation chosen for accuracy: `W` and `grad(P2)` are linear, so
/// `W.W`, `W.grad(S)` and `grad(S).grad(S)` are all quadratic, and
/// region-wise-constant materials add no degree. See `docs/ASSEMBLY_PLAN.md`
/// Sec. 2.
QuadratureRule tet_rule_degree_2();

}  // namespace aphi_solver
