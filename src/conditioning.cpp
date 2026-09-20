#include "aphi_solver/conditioning.hpp"

#include <stdexcept>

namespace aphi_solver {

namespace {

std::vector<Complex> scale_vector(const std::vector<Complex>& v, const Complex& s) {
    std::vector<Complex> out(v.size());
    for (std::size_t i = 0; i < v.size(); ++i) out[i] = v[i] * s;
    return out;
}

/// Copies one block's entries into a combined matrix at (row_offset,
/// col_offset). The destination is still in its triplet phase, so this is a
/// plain append -- no structural search, no dense intermediate.
void scatter_block(SparseMatrixZ& dst, const SparseMatrixZ& src, int row_offset, int col_offset) {
    const auto& row_ptr = src.row_ptr();
    const auto& col_index = src.col_index();
    const auto& values = src.values();
    for (int r = 0; r < src.rows(); ++r) {
        for (int k = row_ptr[static_cast<std::size_t>(r)]; k < row_ptr[static_cast<std::size_t>(r) + 1]; ++k) {
            dst.add(row_offset + r, col_offset + col_index[static_cast<std::size_t>(k)],
                    values[static_cast<std::size_t>(k)]);
        }
    }
}

}  // namespace

APhiBlockSystem apply_symmetric_row_scaling(const APhiBlockSystem& system, double omega) {
    if (omega == 0.0) {
        throw std::invalid_argument(
            "apply_symmetric_row_scaling: omega must be nonzero (the 1/omega scaling is undefined at DC)");
    }
    const Complex inv_j_omega = Complex(1.0, 0.0) / Complex(0.0, omega);
    APhiBlockSystem out = system;
    out.K_PhiA = system.K_PhiA.scaled(inv_j_omega);
    out.K_PhiPhi = system.K_PhiPhi.scaled(inv_j_omega);
    out.rhs_Phi = scale_vector(system.rhs_Phi, inv_j_omega);
    return out;
}

APhiBlockSystem apply_scaled_scalar_potential(const APhiBlockSystem& system, double omega) {
    const Complex j_omega(0.0, omega);
    APhiBlockSystem out = system;
    out.K_APhi = system.K_APhi.scaled(j_omega);
    out.K_PhiPhi = system.K_PhiPhi.scaled(j_omega);
    // K_AA, K_PhiA, rhs_A, rhs_Phi are unchanged: only the Phi' columns are rescaled.
    return out;
}

std::vector<Complex> recover_scaled_scalar_potential(const std::vector<Complex>& phi_prime, double omega) {
    return scale_vector(phi_prime, Complex(0.0, omega));
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

SparseMatrixZ assemble_sparse(const APhiBlockSystem& system) {
    const int nA = system.num_A();
    const int nPhi = system.num_Phi();
    SparseMatrixZ full(nA + nPhi, nA + nPhi);
    full.reserve(system.K_AA.nnz() + system.K_APhi.nnz() + system.K_PhiA.nnz() + system.K_PhiPhi.nnz());
    scatter_block(full, system.K_AA, 0, 0);
    scatter_block(full, system.K_APhi, 0, nA);
    scatter_block(full, system.K_PhiA, nA, 0);
    scatter_block(full, system.K_PhiPhi, nA, nA);
    full.compress();
    return full;
}

ComplexMatrix assemble_dense(const APhiBlockSystem& system) {
    // Deliberately NOT to_dense_matrix(assemble_sparse(system)): this is the
    // ground-truth path, and deriving it from the production path would make
    // tests/test_conditioning.cpp's "the two assemblies agree" check
    // tautological -- a block-offset bug in assemble_sparse would corrupt
    // both sides equally and still compare equal. Placing the four blocks
    // independently here is what gives that cross-check something to catch.
    const int nA = system.num_A();
    const int nPhi = system.num_Phi();
    ComplexMatrix full(nA + nPhi, nA + nPhi);
    full.set_block(0, 0, to_dense_matrix(system.K_AA));
    full.set_block(0, nA, to_dense_matrix(system.K_APhi));
    full.set_block(nA, 0, to_dense_matrix(system.K_PhiA));
    full.set_block(nA, nA, to_dense_matrix(system.K_PhiPhi));
    return full;
}

ComplexMatrix assemble_dense_rhs(const APhiBlockSystem& system) {
    std::vector<Complex> stacked;
    stacked.reserve(system.rhs_A.size() + system.rhs_Phi.size());
    stacked.insert(stacked.end(), system.rhs_A.begin(), system.rhs_A.end());
    stacked.insert(stacked.end(), system.rhs_Phi.begin(), system.rhs_Phi.end());
    return to_column(stacked);
}

}  // namespace aphi_solver
