#pragma once

#include <algorithm>
#include <complex>
#include <stdexcept>
#include <string>
#include <vector>

#include "aphi_solver/sparse_matrix.hpp"

namespace aphi_solver {

/// A symmetric sparse matrix, storing only the upper triangle in CSR.
///
/// `A = A^T` exactly, not `A = A^H`: the A-Phi system's symmetric forms are
/// **complex symmetric**. That permits a complex `LDL^T` but not Cholesky,
/// and COCG rather than CG. Nothing here assumes a real or positive matrix.
///
/// **Why a separate class rather than a flag on `Sparse<T>`.** An
/// upper-triangle CSR array is indistinguishable from a general CSR array by
/// inspection, so handing one to a general `matvec` gives a wrong answer with
/// no complaint -- exactly the kind of silent half-loss this project keeps
/// finding. A distinct type makes that a compile error instead.
///
/// On `cylinder_box.msh` at AC this stores 797445 of the full matrix's
/// 1557522 entries: 12.2 MB against 23.8 MB for complex values.
template <typename T>
class SparseSymmetric {
public:
    SparseSymmetric() = default;

    /// Adopts an upper-triangle CSR structure with every value zero.
    ///
    /// Validates squareness, that each row's columns ascend strictly and lie
    /// in range, and -- the one that matters here -- that no column precedes
    /// its row. A lower-triangle entry would be stored but never reached by
    /// `matvec`, which mirrors instead.
    static SparseSymmetric<T> from_pattern(int rows, std::vector<int> row_ptr,
                                           std::vector<int> col_index) {
        if (rows < 0) throw std::invalid_argument("SparseSymmetric::from_pattern: negative rows");
        if (static_cast<int>(row_ptr.size()) != rows + 1) {
            throw std::invalid_argument(
                "SparseSymmetric::from_pattern: row_ptr must have rows + 1 entries");
        }
        if (row_ptr.front() != 0 || row_ptr.back() != static_cast<int>(col_index.size())) {
            throw std::invalid_argument(
                "SparseSymmetric::from_pattern: row_ptr must run from 0 to col_index.size()");
        }
        for (int r = 0; r < rows; ++r) {
            const int begin = row_ptr[static_cast<std::size_t>(r)];
            const int end = row_ptr[static_cast<std::size_t>(r) + 1];
            if (end < begin) {
                throw std::invalid_argument(
                    "SparseSymmetric::from_pattern: row_ptr is not non-decreasing");
            }
            for (int k = begin; k < end; ++k) {
                const int c = col_index[static_cast<std::size_t>(k)];
                if (c < 0 || c >= rows) {
                    throw std::invalid_argument(
                        "SparseSymmetric::from_pattern: a column index is out of range");
                }
                if (c < r) {
                    throw std::invalid_argument(
                        "SparseSymmetric::from_pattern: entry (" + std::to_string(r) + ", " +
                        std::to_string(c) +
                        ") is below the diagonal. Only the upper triangle is stored; a "
                        "lower-triangle entry would never be read, since matvec mirrors.");
                }
                if (k > begin && c <= col_index[static_cast<std::size_t>(k - 1)]) {
                    throw std::invalid_argument(
                        "SparseSymmetric::from_pattern: each row's columns must be strictly "
                        "ascending");
                }
            }
        }

        SparseSymmetric<T> m;
        m.rows_ = rows;
        m.row_ptr_ = std::move(row_ptr);
        m.col_index_ = std::move(col_index);
        m.values_.assign(m.col_index_.size(), T{});
        return m;
    }

    int rows() const { return rows_; }
    int cols() const { return rows_; }

    /// Entries **stored**, which is the upper triangle. The full matrix has
    /// `2 * nnz() - (diagonal entries)`.
    std::size_t nnz() const { return values_.size(); }

    const std::vector<int>& row_ptr() const { return row_ptr_; }
    const std::vector<int>& col_index() const { return col_index_; }
    const std::vector<T>& values() const { return values_; }
    std::vector<T>& mutable_values() { return values_; }

    /// Where (row, col) is stored, or -1. Either order works: the pair is
    /// normalised, since `A(i,j)` and `A(j,i)` are the same number and only
    /// one is held.
    int find_slot(int row, int col) const {
        int r = row, c = col;
        if (c < r) std::swap(r, c);
        if (r < 0 || r >= rows_) return -1;
        const auto begin = col_index_.begin() + row_ptr_[static_cast<std::size_t>(r)];
        const auto end = col_index_.begin() + row_ptr_[static_cast<std::size_t>(r) + 1];
        const auto it = std::lower_bound(begin, end, c);
        if (it == end || *it != c) return -1;
        return static_cast<int>(it - col_index_.begin());
    }

    T entry(int row, int col) const {
        const int slot = find_slot(row, col);
        return slot < 0 ? T{} : values_[static_cast<std::size_t>(slot)];
    }

    /// `y = A x` for the full symmetric matrix. Each stored off-diagonal
    /// entry contributes twice, once to each of the rows it belongs to; the
    /// diagonal contributes once. Getting that wrong is the whole reason this
    /// is a type and not a flag.
    std::vector<T> matvec(const std::vector<T>& x) const {
        if (static_cast<int>(x.size()) != rows_) {
            throw std::invalid_argument("SparseSymmetric::matvec: x has the wrong length");
        }
        std::vector<T> y(static_cast<std::size_t>(rows_), T{});
        for (int r = 0; r < rows_; ++r) {
            const std::size_t ur = static_cast<std::size_t>(r);
            for (int k = row_ptr_[ur]; k < row_ptr_[ur + 1]; ++k) {
                const std::size_t uk = static_cast<std::size_t>(k);
                const int c = col_index_[uk];
                const T a = values_[uk];
                y[ur] += a * x[static_cast<std::size_t>(c)];
                if (c != r) y[static_cast<std::size_t>(c)] += a * x[ur];
            }
        }
        return y;
    }

    /// The same matrix with both triangles stored, for anything that needs a
    /// general matrix -- and for tests, where comparing against an
    /// independently assembled full matrix is the strongest check available.
    Sparse<T> to_full() const {
        Sparse<T> out(rows_, rows_);
        for (int r = 0; r < rows_; ++r) {
            const std::size_t ur = static_cast<std::size_t>(r);
            for (int k = row_ptr_[ur]; k < row_ptr_[ur + 1]; ++k) {
                const std::size_t uk = static_cast<std::size_t>(k);
                const int c = col_index_[uk];
                out.add(r, c, values_[uk]);
                if (c != r) out.add(c, r, values_[uk]);
            }
        }
        out.compress();
        return out;
    }

private:
    int rows_ = 0;
    std::vector<int> row_ptr_;
    std::vector<int> col_index_;
    std::vector<T> values_;
};

using SparseSymmetricD = SparseSymmetric<double>;
using SparseSymmetricZ = SparseSymmetric<std::complex<double>>;

}  // namespace aphi_solver
