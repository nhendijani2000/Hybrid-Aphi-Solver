#pragma once

#include <complex>

#include "aphi_solver/basis_functions.hpp"
#include "aphi_solver/problem_binding.hpp"

namespace aphi_solver {

/// Element matrices: the per-tetrahedron blocks of the A-Phi system.
///
/// See `docs/ASSEMBLY_PLAN.md`. The kernels below are **pure functions of
/// the geometry** -- no `Mesh`, no tet index, no global state -- so they can
/// be checked against hand-computed numbers with a `TetGeometry` built by
/// hand. They use the LOCAL Whitney basis; the global edge orientation is
/// applied later, by the scatter, through `DofMap`'s `coeff` (Sec. 9.1).
///
/// Each writes into a caller-owned buffer. Nothing here allocates, and
/// nothing is returned by value: in assembly these are stack arrays reused
/// for every tet.

/// The material and frequency coefficients for one body, in the form the
/// weak form actually uses (`docs/FORMULATION.md` Sec. 1.4):
///
/// ```
///   nu    = 1 / mu                        the curl-curl coefficient
///   alpha = j*w*sigma - w^2*eps           on A, in both equations
///   beta  = sigma + j*w*eps               on grad(Phi), in both equations
/// ```
///
/// with the identity `alpha = j*w*beta`, which is where the system's
/// asymmetry comes from and what the conditioning formulations exploit.
struct ElementCoefficients {
    double nu = 0.0;
    std::complex<double> alpha;
    std::complex<double> beta;
};

/// Builds them for one body at one angular frequency.
///
/// At `omega == 0` this gives `alpha = 0` and `beta = sigma` exactly, which
/// is the DC decoupling of `FORMULATION.md` Sec. 2 falling out of the
/// coefficients rather than needing a branch: the mass term and the Phi-A
/// coupling vanish, leaving magnetostatics and DC conduction side by side.
ElementCoefficients element_coefficients(const BoundBody& body, double omega);

/// The A-A block, 6x6 row-major, 36 entries:
///
/// ```
///   K[i][j] = nu * integral (curl W_i) . (curl W_j)
///           + alpha * integral W_i . W_j
/// ```
///
/// Both terms are A-A and *add* in the assembled matrix, so they are
/// computed together: one block and one scatter rather than two.
///
/// The curl-curl part needs no quadrature. `curl W = 2 grad(Li) x grad(Lj)`
/// is constant over the tet, so that term is exactly
/// `nu * V * (curl W_i) . (curl W_j)`. The mass term is quadratic and uses
/// the degree-2 rule, which is exact for it.
void kernel_AA(const TetGeometry& g, const ElementCoefficients& c, std::complex<double>* out);

/// The A-Phi coupling block, 6x10 row-major (rows are the tet's 6 edges,
/// columns its 10 P2 nodes), 60 entries:
///
/// ```
///   C[i][b] = beta * integral W_i . grad(S_b)
/// ```
///
/// This is `C_APhi` of `docs/ASSEMBLY_PLAN.md` Sec. 2 -- the block the
/// scatter uses for BOTH couplings. The Phi-A block is not computed
/// separately, because `alpha = j*omega*beta` makes
/// `C_PhiA = j*omega * C_APhi^T` an identity rather than a coincidence.
///
/// Quadratic (W is linear, grad(P2) is linear), so the degree-2 rule is
/// exact for it.
void kernel_APhi(const TetGeometry& g, const ElementCoefficients& c, std::complex<double>* out);

/// The Phi-Phi block, 10x10 row-major, 100 entries:
///
/// ```
///   K[a][b] = beta * integral grad(S_a) . grad(S_b)
/// ```
///
/// Symmetric, positive semi-definite for real positive `beta`, and singular
/// by exactly one dimension: `sum_a S_a = 1`, so `sum_a grad(S_a) = 0` and
/// the constant vector is in its null space. That is the usual stiffness
/// matrix null space, removed globally by a potential reference rather than
/// here.
///
/// Quadratic, so the degree-2 rule is exact.
void kernel_PhiPhi(const TetGeometry& g, const ElementCoefficients& c, std::complex<double>* out);

}  // namespace aphi_solver
