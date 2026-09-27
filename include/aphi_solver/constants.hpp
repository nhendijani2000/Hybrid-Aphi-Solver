#pragma once

namespace aphi_solver {

/// Physical constants, in SI. The mesh is in metres by the time assembly or
/// post-processing runs, whatever `length_unit` the input file used.
///
/// These live here rather than as file-local constants because assembly and
/// post-processing must agree on them exactly: `H = B / mu` computed with a
/// different `mu_0` than the one that built `nu = 1 / mu` is wrong in a way no
/// residual would reveal, and a duplicated literal is how that happens.
constexpr double kMu0 = 4.0e-7 * 3.14159265358979323846;
constexpr double kEps0 = 8.8541878128e-12;

}  // namespace aphi_solver
