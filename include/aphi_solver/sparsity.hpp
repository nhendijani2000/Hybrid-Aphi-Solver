#pragma once

#include <vector>

#include "aphi_solver/dof_map.hpp"
#include "aphi_solver/mesh.hpp"
#include "aphi_solver/problem_binding.hpp"

namespace aphi_solver {

/// Where the assembled matrix's nonzeros are, worked out before any value
/// exists. See `docs/ASSEMBLY_PLAN.md` §3.
///
/// **Why this is a separate pass.** The alternative is the triplet path
/// `Sparse<T>` already offers: `add()` every contribution, then `compress()`
/// sorts and sums duplicates. That works, and stays the right tool for the
/// small hand-built matrices the tests use, but for assembly it holds every
/// contribution in memory at once and sorts all of them.
///
/// It also cannot be reused. The pattern depends only on the mesh's
/// connectivity and which DOFs exist; the *values* depend on materials and
/// frequency. A 41-point sweep therefore builds this once and refills values
/// 41 times, where the triplet path would rebuild and re-sort every time.
struct SparsityPattern {
    int rows = 0;
    int cols = 0;

    /// CSR structure. `col_index[row_ptr[r] .. row_ptr[r+1])` are row `r`'s
    /// columns, **ascending**, each appearing once -- which is what lets
    /// `find_slot` binary-search them.
    std::vector<int> row_ptr;
    std::vector<int> col_index;

    std::size_t nnz() const { return col_index.size(); }

    /// The position in `col_index` (and so in a values array laid out
    /// alongside it) of entry (row, col), or -1 if the pattern has no such
    /// entry.
    ///
    /// A -1 during assembly means the symbolic and numeric passes disagree
    /// about the matrix's shape, which must abort rather than quietly drop
    /// a term -- see `ASSEMBLY_PLAN.md` §7.4.
    int find_slot(int row, int col) const;
};

/// Builds the pattern from the DOF map: for every tet, every ordered pair of
/// its live global DOFs is a possible nonzero.
///
/// **Ordered** pairs, both (i,j) and (j,i). The assembled matrix is not
/// symmetric in general -- `C_PhiA = jω·C_APhiᵀ`, not `C_APhiᵀ` -- so a
/// pattern storing one triangle would leave every Φ–A contribution with
/// nowhere to go. That is a silent half-loss of the coupling, and it is the
/// easiest mistake to make here.
///
/// Built by counting into `row_ptr`, prefix-summing, filling, then sorting
/// and uniquing each row in place: no `std::set`, and no per-row
/// `std::vector`, either of which would allocate once per row.
SparsityPattern build_sparsity(const DofMap& dofs, const BoundProblem& bound, const Mesh& mesh);

}  // namespace aphi_solver
