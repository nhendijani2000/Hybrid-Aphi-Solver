#include "aphi_solver/sparsity.hpp"

#include <algorithm>
#include <array>
#include <stdexcept>
#include <string>

namespace aphi_solver {

int SparsityPattern::find_slot(int row, int col) const {
    if (row < 0 || row >= rows) return -1;
    const auto begin = col_index.begin() + row_ptr[static_cast<std::size_t>(row)];
    const auto end = col_index.begin() + row_ptr[static_cast<std::size_t>(row) + 1];
    const auto it = std::lower_bound(begin, end, col);
    if (it == end || *it != col) return -1;
    return static_cast<int>(it - col_index.begin());
}

SparsityPattern build_sparsity(const DofMap& dofs, const BoundProblem& bound, const Mesh& mesh,
                               SparsityStorage storage) {
    SparsityPattern p;
    p.rows = dofs.num_total;
    p.cols = dofs.num_total;
    p.upper_only = storage == SparsityStorage::UpperTriangle;
    p.row_ptr.assign(static_cast<std::size_t>(p.rows) + 1, 0);

    // A tet's live DOFs: at most 16 (6 edges + 10 P2 nodes), fewer once
    // Dirichlet edges, tree edges and absent Phi nodes drop out. Several
    // local DOFs can map to the SAME global index -- every node of a port
    // terminal reads that port's single unknown -- so the list is uniqued
    // before pairs are formed, or the duplicate columns would be counted
    // several times.
    std::array<int, 16> live{};

    const auto gather = [&](int tet, int& count) {
        const TetDofs d = dofs.local_dofs(tet, mesh, bound);
        count = 0;
        for (const DofEntry& e : d.edge) {
            if (e.index >= 0) live[static_cast<std::size_t>(count++)] = e.index;
        }
        for (const DofEntry& e : d.phi) {
            if (e.index >= 0) live[static_cast<std::size_t>(count++)] = e.index;
        }
        std::sort(live.begin(), live.begin() + count);
        count = static_cast<int>(std::unique(live.begin(), live.begin() + count) -
                                 live.begin());
    };

    // Pass 1: count entries per row. Upper bounds, with duplicates across
    // tets still counted -- they are removed after the fill.
    for (int t = 0; t < mesh.num_tets(); ++t) {
        int n = 0;
        gather(t, n);
        for (int i = 0; i < n; ++i) {
            // Upper triangle: row `live[i]` keeps only the columns at or
            // after it, and the list is sorted, so that is the tail from i.
            const int kept = p.upper_only ? n - i : n;
            p.row_ptr[static_cast<std::size_t>(live[static_cast<std::size_t>(i)]) + 1] += kept;
        }
    }
    // A voltage port's row gets nothing from the tets: its terminal is a
    // prescribed Phi, so the port has no column and no current-balance
    // equation was scattered into its row. Assembly writes the constraint
    // `V = V_given` there instead, which needs a diagonal slot.
    for (std::size_t k = 0; k < dofs.port_is_fixed.size(); ++k) {
        if (!dofs.port_is_fixed[k]) continue;
        p.row_ptr[static_cast<std::size_t>(dofs.port_index[k]) + 1] += 1;
    }
    for (int r = 0; r < p.rows; ++r) {
        p.row_ptr[static_cast<std::size_t>(r) + 1] += p.row_ptr[static_cast<std::size_t>(r)];
    }

    // Pass 2: fill. `cursor` walks each row's slice as entries land in it.
    std::vector<int> cursor(p.row_ptr.begin(), p.row_ptr.end() - 1);
    p.col_index.assign(static_cast<std::size_t>(p.row_ptr.back()), 0);
    for (int t = 0; t < mesh.num_tets(); ++t) {
        int n = 0;
        gather(t, n);
        for (int i = 0; i < n; ++i) {
            const int row = live[static_cast<std::size_t>(i)];
            for (int j = p.upper_only ? i : 0; j < n; ++j) {
                // Both (i,j) and (j,i): the matrix is not symmetric in
                // general, so one triangle would lose the Phi-A coupling.
                // Pass 1 counted this row's upper bound; writing past it is a
                // heap overrun, which is how a miscount in pass 1 showed up
                // when it was tried as a negative control. One comparison per
                // entry turns that into a diagnosable exception, and
                // build_sparsity runs once per mesh so it costs nothing that
                // matters.
                if (cursor[static_cast<std::size_t>(row)] >=
                    p.row_ptr[static_cast<std::size_t>(row) + 1]) {
                    throw std::logic_error(
                        "build_sparsity: row " + std::to_string(row) +
                        " received more entries than the counting pass reserved for it. The "
                        "count and the fill disagree.");
                }
                p.col_index[static_cast<std::size_t>(cursor[static_cast<std::size_t>(row)]++)] =
                    live[static_cast<std::size_t>(j)];
            }
        }
    }

    for (std::size_t k = 0; k < dofs.port_is_fixed.size(); ++k) {
        if (!dofs.port_is_fixed[k]) continue;
        const int row = dofs.port_index[k];
        if (cursor[static_cast<std::size_t>(row)] >=
            p.row_ptr[static_cast<std::size_t>(row) + 1]) {
            throw std::logic_error("build_sparsity: no room reserved for port row " +
                                   std::to_string(row) + "'s constraint diagonal.");
        }
        p.col_index[static_cast<std::size_t>(cursor[static_cast<std::size_t>(row)]++)] = row;
    }

    // Pass 3: sort and unique each row in place, compacting as we go. The
    // rows are short (tens of entries), so this is a handful of comparisons
    // each rather than a sort of the whole array.
    std::vector<int> compacted_ptr(static_cast<std::size_t>(p.rows) + 1, 0);
    int write = 0;
    for (int r = 0; r < p.rows; ++r) {
        const int begin = p.row_ptr[static_cast<std::size_t>(r)];
        const int end = cursor[static_cast<std::size_t>(r)];
        std::sort(p.col_index.begin() + begin, p.col_index.begin() + end);
        const auto last = std::unique(p.col_index.begin() + begin, p.col_index.begin() + end);

        compacted_ptr[static_cast<std::size_t>(r)] = write;
        for (auto it = p.col_index.begin() + begin; it != last; ++it) {
            p.col_index[static_cast<std::size_t>(write++)] = *it;
        }
    }
    compacted_ptr[static_cast<std::size_t>(p.rows)] = write;

    p.col_index.resize(static_cast<std::size_t>(write));
    p.col_index.shrink_to_fit();
    p.row_ptr = std::move(compacted_ptr);
    return p;
}

}  // namespace aphi_solver
