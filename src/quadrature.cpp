#include "aphi_solver/quadrature.hpp"

#include <cmath>

namespace aphi_solver {

QuadratureRule tet_rule_degree_1() {
    QuadratureRule r;
    r.degree = 1;
    r.count = 1;
    r.points[0] = {0.25, 0.25, 0.25, 0.25};
    r.weights[0] = 1.0;
    return r;
}

QuadratureRule tet_rule_degree_2() {
    QuadratureRule r;
    r.degree = 2;
    r.count = 4;

    // (5 + 3*sqrt(5))/20 and (5 - sqrt(5))/20. Written as expressions rather
    // than decimal literals so they are exact to the last bit and their
    // origin is visible; a mistyped digit in a pasted constant is a silent
    // loss of accuracy that only an exactness test would catch.
    const double s5 = std::sqrt(5.0);
    const double a = (5.0 + 3.0 * s5) / 20.0;
    const double b = (5.0 - s5) / 20.0;

    r.points[0] = {a, b, b, b};
    r.points[1] = {b, a, b, b};
    r.points[2] = {b, b, a, b};
    r.points[3] = {b, b, b, a};
    for (int k = 0; k < 4; ++k) r.weights[static_cast<std::size_t>(k)] = 0.25;
    return r;
}

}  // namespace aphi_solver
