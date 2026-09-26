#pragma once

#include <complex>
#include <cstddef>
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

/// Every tet's scatter positions, computed once and reused at every
/// frequency -- candidate A of `docs/ASSEMBLY_PLAN.md` Sec. 11.
///
/// `refill` otherwise works out where each of a tet's ~225 entries belongs
/// as it goes, which measured 43 % of its time. That bookkeeping depends on
/// the mesh and the DOF map only, never on omega, so across a frequency
/// sweep it is the same work repeated. This holds the answer instead.
///
/// **Optional, and deliberately so.** The cost is `live^2` ints per tet, up
/// to about 1 KB, and it scales with the mesh rather than with the sweep.
/// `bytes()` reports the real figure. Pass it to `refill` for a sweep; leave
/// it out for a single solve, where building it would cost more than it
/// saves.
class ScatterMap {
public:
    /// Builds the map for every tet of the mesh. The bodies partition the
    /// tets, so indexing by tet id covers exactly what `refill` visits.
    ///
    /// Throws std::logic_error on a pair the pattern has no slot for -- the
    /// same condition, and for the same reason, as `refill`.
    static ScatterMap build(const SparsityPattern& pattern, const DofMap& dofs, const Mesh& mesh,
                            const BoundProblem& bound);

    /// Heap bytes held, for reporting before deciding to build one.
    std::size_t bytes() const;

    int num_tets() const { return static_cast<int>(live_.size()); }

    /// Distinct live global DOFs of one tet, at most 16.
    int live_count(int tet) const { return live_[static_cast<std::size_t>(tet)]; }

    /// Local edge (0..5) and local P2 node (0..9) to its position among this
    /// tet's live DOFs, or -1 when that DOF is eliminated.
    const int* edge_at(int tet) const { return &pos_[static_cast<std::size_t>(tet) * 16]; }
    const int* phi_at(int tet) const { return &pos_[static_cast<std::size_t>(tet) * 16 + 6]; }

    /// This tet's `live_count` x `live_count` block of value-array indices,
    /// row-major.
    const int* slots(int tet) const {
        return slot_.data() + offset_[static_cast<std::size_t>(tet)];
    }

private:
    std::vector<int> live_;          ///< per tet: distinct live DOFs
    std::vector<int> pos_;           ///< per tet: 6 edge positions then 10 Phi
    std::vector<std::size_t> offset_;  ///< per tet + 1: its block's start in slot_
    std::vector<int> slot_;          ///< sum over tets of live^2
};

/// What one assembly produced.
struct AssembledSystem {
    SparseMatrixZ matrix;
    std::vector<std::complex<double>> rhs;
};

/// Allocates a system for `pattern` with every value zero. The structure is
/// validated here, so once per mesh rather than once per frequency.
///
/// Separate from `refill` because allocating it is not cheap: 23.8 MB of
/// fresh pages on the cylinder and the first write to each, which measured
/// 11.4 ms against 0.6 ms to re-zero an array already held
/// (`docs/ASSEMBLY_PLAN.md` Sec. 12). A sweep should do it once.
AssembledSystem make_system(const SparsityPattern& pattern, int num_unknowns);

/// Assembles into an existing system, **clearing it first**. Safe to call
/// repeatedly on one `AssembledSystem` at different frequencies, which is
/// the point: the matrix's structure is frequency-independent even though
/// essentially every value is not.
///
/// `scatter`, when given, must have been built from this same `pattern`,
/// `dofs`, `mesh` and `bound`; it replaces working out each entry's position.
///
/// One pass over the bodies, then over each body's tets, evaluating the
/// three kernels into stack buffers and scattering them straight into
/// `pattern`'s slots. Nothing per-element is kept.
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
void refill(AssembledSystem& system, const BoundProblem& bound, const Mesh& mesh,
            const DofMap& dofs, const SparsityPattern& pattern, double omega,
            Formulation3 formulation, const ScatterMap* scatter = nullptr);

/// `make_system` then `refill`, for a single solve. A sweep should call the
/// two separately and keep the system between frequencies.
AssembledSystem assemble(const BoundProblem& bound, const Mesh& mesh, const DofMap& dofs,
                         const SparsityPattern& pattern, double omega, Formulation3 formulation);

}  // namespace aphi_solver
