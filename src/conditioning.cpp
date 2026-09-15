#include "aphi_solver/conditioning.hpp"

#include <stdexcept>

namespace aphi_solver {

APhiBlockSystem apply_symmetric_row_scaling(const APhiBlockSystem& system, double omega) {
    if (omega == 0.0) {
        throw std::invalid_argument(
            "apply_symmetric_row_scaling: omega must be nonzero (the 1/omega scaling is undefined at DC)");
    }
    const Complex j_omega(0.0, omega);
    APhiBlockSystem out = system;
    out.K_PhiA = system.K_PhiA / j_omega;
    out.K_PhiPhi = system.K_PhiPhi / j_omega;
    out.rhs_Phi = system.rhs_Phi / j_omega;
    return out;
}

APhiBlockSystem apply_scaled_scalar_potential(const APhiBlockSystem& system, double omega) {
    const Complex j_omega(0.0, omega);
    APhiBlockSystem out = system;
    out.K_APhi = system.K_APhi * j_omega;
    out.K_PhiPhi = system.K_PhiPhi * j_omega;
    // K_AA, K_PhiA, rhs_A, rhs_Phi are unchanged: only the Phi' columns are rescaled.
    return out;
}

ComplexMatrix recover_scaled_scalar_potential(const ComplexMatrix& phi_prime, double omega) {
    const Complex j_omega(0.0, omega);
    return phi_prime * j_omega;
}

APhiBlockSystem apply_conditioning(const APhiBlockSystem& system, ConditioningStrategy strategy, double omega) {
    switch (strategy) {
        case ConditioningStrategy::Natural:
            return system;
        case ConditioningStrategy::SymmetricRowScaled:
            return apply_symmetric_row_scaling(system, omega);
        case ConditioningStrategy::ScaledScalarPotential:
            return apply_scaled_scalar_potential(system, omega);
    }
    throw std::invalid_argument("apply_conditioning: unknown strategy");
}

ConditioningStrategy recommend_strategy(double frequency_hz, double crossover_hz) {
    return frequency_hz < crossover_hz ? ConditioningStrategy::SymmetricRowScaled
                                        : ConditioningStrategy::ScaledScalarPotential;
}

ComplexMatrix assemble_dense(const APhiBlockSystem& system) {
    const int nA = system.num_A();
    const int nPhi = system.num_Phi();
    ComplexMatrix full(nA + nPhi, nA + nPhi);
    full.set_block(0, 0, system.K_AA);
    full.set_block(0, nA, system.K_APhi);
    full.set_block(nA, 0, system.K_PhiA);
    full.set_block(nA, nA, system.K_PhiPhi);
    return full;
}

}  // namespace aphi_solver
