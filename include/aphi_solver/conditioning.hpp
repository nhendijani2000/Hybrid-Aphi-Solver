#pragma once

#include "aphi_solver/complex_matrix.hpp"

namespace aphi_solver {

/// The assembled frequency-domain A-Phi block system, before any conditioning
/// transform is applied:
///
///   [ K_AA    K_APhi  ] [ a   ]   [ rhs_A   ]
///   [ K_PhiA  K_PhiPhi] [ Phi ] = [ rhs_Phi ]
///
/// K_AA carries the (1/mu) curl-curl term plus the (j*w*sigma - w^2*eps) mass
/// term; K_APhi / K_PhiA carry the (sigma + j*w*eps) A-Phi coupling; K_PhiPhi is
/// the scalar block. This is the standard frequency-domain A-Phi system used
/// throughout the public literature (e.g. Dular et al. 2000; Zhao & Fu 2017) --
/// nothing about this block structure is specific to any one implementation.
/// The blocks are sparse as of Phase 03.5 (`docs/ROADMAP.md`): they were
/// dense `ComplexMatrix` while this module only ever saw hand-built test
/// systems, which does not survive contact with a real mesh -- K_AA alone
/// is ~55 MB dense at `meshes/cube_6.msh` and impossible beyond that. The
/// right-hand sides stay dense vectors, because they are dense.
struct APhiBlockSystem {
    SparseMatrixZ K_AA;
    SparseMatrixZ K_APhi;
    SparseMatrixZ K_PhiA;
    SparseMatrixZ K_PhiPhi;
    std::vector<Complex> rhs_A;    // length n_A
    std::vector<Complex> rhs_Phi;  // length n_Phi

    int num_A() const { return K_AA.rows(); }
    int num_Phi() const { return K_PhiPhi.rows(); }
};

/// Conditioning strategy for the coupled A-Phi system:
///
///  - Natural: the system as assembled, no transform.
///  - SymmetricRowScaled: divides the entire scalar-potential row (K_PhiA,
///    K_PhiPhi, rhs_Phi) by j*omega. This makes the two off-diagonal blocks
///    transposes of one another (better structural symmetry) at the cost of a
///    1/omega scaling that blows up as omega -> 0.
///  - ScaledScalarPotential: substitutes Phi = j*omega*Phi' throughout, which
///    instead rescales the Phi *columns* (K_APhi, K_PhiPhi) by j*omega and never
///    divides by omega. The linear solve then returns Phi', which must be mapped
///    back to the physical Phi with `recover_scaled_scalar_potential`.
///
/// Both transforms are direct algebraic consequences of the block system above
/// -- see docs/CONDITIONING.md for the derivation of each. Neither strategy is
/// hard-wired to a particular frequency range: which one behaves better, and at
/// what frequency the crossover falls, is a property of your own mesh and
/// materials and should be measured (Phase 05 of the roadmap), not assumed.
enum class ConditioningStrategy {
    Natural,
    SymmetricRowScaled,
    ScaledScalarPotential,
};

/// Divides the scalar-potential row (K_PhiA, K_PhiPhi, rhs_Phi) by j*omega.
/// Throws std::invalid_argument if omega == 0 (the transform is undefined at DC).
APhiBlockSystem apply_symmetric_row_scaling(const APhiBlockSystem& system, double omega);

/// Substitutes Phi = j*omega*Phi': rescales K_APhi and K_PhiPhi by j*omega and
/// leaves everything else unchanged. Solve the returned system for Phi', then
/// call `recover_scaled_scalar_potential` to get the physical Phi.
APhiBlockSystem apply_scaled_scalar_potential(const APhiBlockSystem& system, double omega);

/// Recovers the physical scalar potential Phi = j*omega*Phi' after solving a
/// system produced by `apply_scaled_scalar_potential`.
std::vector<Complex> recover_scaled_scalar_potential(const std::vector<Complex>& phi_prime, double omega);

/// Dispatches to the requested strategy. `Natural` returns `system` unchanged.
APhiBlockSystem apply_conditioning(const APhiBlockSystem& system, ConditioningStrategy strategy, double omega);

/// Picks SymmetricRowScaled below `crossover_hz`, ScaledScalarPotential at or
/// above it. There is no built-in default for `crossover_hz`: determine it
/// empirically for your own mesh and materials via the frequency sweep in
/// Phase 05 of the project roadmap (plot condition number / iterations vs.
/// frequency for both strategies and read off where they cross). Do not assume
/// a fixed number without measuring it on your own problem.
ConditioningStrategy recommend_strategy(double frequency_hz, double crossover_hz);

/// Assembles the dense (n_A + n_Phi) x (n_A + n_Phi) system matrix, block order
/// [a; Phi]. Useful for condition-number analysis and for the direct-solver
/// baseline in Phase 05; not intended for production-size meshes.
ComplexMatrix assemble_dense(const APhiBlockSystem& system);

/// The matching stacked right-hand side [rhs_A; rhs_Phi] as an
/// (n_A + n_Phi) x 1 column, so `solve_dense(assemble_dense(sys),
/// assemble_dense_rhs(sys))` is the complete small-system baseline solve.
ComplexMatrix assemble_dense_rhs(const APhiBlockSystem& system);

/// Assembles the whole system as ONE sparse matrix in the same [a; Phi]
/// block order -- the form a sparse direct solver is handed (Phase 05).
/// Unlike `assemble_dense` this is a production path, so it is built
/// block-by-block through the triplet phase rather than through any dense
/// intermediate.
SparseMatrixZ assemble_sparse(const APhiBlockSystem& system);

}  // namespace aphi_solver
