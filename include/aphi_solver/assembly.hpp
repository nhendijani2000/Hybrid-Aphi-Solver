#pragma once

#include <complex>
#include <vector>

#include "aphi_solver/dof_map.hpp"
#include "aphi_solver/element_matrix.hpp"
#include "aphi_solver/sparse_matrix.hpp"
#include "aphi_solver/sparsity.hpp"

namespace aphi_solver {

/// Global assembly: the numeric pass. See `docs/ASSEMBLY_PLAN.md` §4.

/// Which of the three conditioning formulations to assemble
/// (`docs/CONDITIONING.md`). They differ only in a row scale `r` applied to
/// every Phi-row contribution and a column scale `c` from the substitution
/// `Phi = c * Phi'`:
///
/// | | r | c | (A,Phi) | (Phi,A) | symmetric |
/// |---|---|---|---|---|---|
/// | Natural | 1 | 1 | beta*C | alpha*C^T | no |
/// | RowScaled | 1/(j*w) | 1 | beta*C | beta*C^T | yes |
/// | ScaledPhi | 1 | j*w | alpha*C | alpha*C^T | yes |
///
/// Both symmetric cases work because `alpha = j*w*beta`: RowScaled divides
/// the `j*w` out of the Phi row, ScaledPhi supplies it to the Phi column
/// instead. The condition is exactly `c == r * j*w`.
enum class Formulation3 {
    Natural,    ///< as assembled; needs a general (LU) solver
    RowScaled,  ///< CONDITIONING.md Formulation 1
    ScaledPhi   ///< CONDITIONING.md Formulation 2; the unknown becomes Phi'
};

/// The scales for one formulation at one frequency.
struct FormulationScales {
    std::complex<double> row{1.0, 0.0};     ///< r
    std::complex<double> column{1.0, 0.0};  ///< c

    /// True when the assembled matrix should come out symmetric, which is
    /// exactly `c == r * j*omega`.
    bool symmetric = false;
};

/// Throws std::invalid_argument for RowScaled or ScaledPhi at omega == 0:
/// both need a division by `j*omega` somewhere (the row here, a prescribed
/// value in ScaledPhi's case), and at DC the system decouples anyway so
/// there is no coupling block to symmetrise.
FormulationScales formulation_scales(Formulation3 f, double omega);

/// What one assembly produced.
struct AssembledSystem {
    SparseMatrixZ matrix;
    std::vector<std::complex<double>> rhs;
};

/// Assembles the global system: one pass over the bodies, then over each
/// body's tets, evaluating the three kernels into stack buffers and
/// scattering them straight into `pattern`'s slots. Nothing per-element is
/// kept.
///
/// Material coefficients are computed **once per body**, which is why the
/// loop is bodies-then-tets rather than over tets directly.
///
/// A prescribed DOF has no column, so its contribution moves to the
/// right-hand side instead -- carrying `Phi_given / c`, the value of the
/// SCALED unknown, which cancels the `c` in that column's coefficient.
///
/// The two port types are duals. A current port's terminal potential is an
/// unknown and its row is the terminal's current balance, so the driving
/// current goes into that row's right-hand side. A voltage port's terminal
/// potential is prescribed, so it has no column and nothing is scattered
/// into its row; the row carries the constraint `V' = V_given / c` with a
/// unit diagonal, which keeps the matrix square and symmetric. The price is
/// that the port current is no longer a residual of that row -- extraction
/// must re-form the terminal's current balance.
///
/// Throws std::logic_error if a scatter finds no slot: that means the
/// symbolic and numeric passes disagree about the matrix's shape, and
/// dropping the term instead would give a quietly wrong matrix.
AssembledSystem assemble(const BoundProblem& bound, const Mesh& mesh, const DofMap& dofs,
                         const SparsityPattern& pattern, double omega, Formulation3 formulation);

}  // namespace aphi_solver
