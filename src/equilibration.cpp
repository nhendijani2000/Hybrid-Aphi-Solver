#include "aphi_solver/equilibration.hpp"

#include <cmath>
#include <stdexcept>

namespace aphi_solver {

std::vector<double> compute_symmetric_equilibration(const SparseMatrixZ& A, int iterations) {
    const int n = A.rows();
    if (A.cols() != n) {
        throw std::invalid_argument("compute_symmetric_equilibration: A must be square");
    }

    std::vector<double> d(static_cast<std::size_t>(n), 1.0);

    const auto& row_ptr = A.row_ptr();
    const auto& col_index = A.col_index();
    // Only the values change between iterations -- the sparsity pattern is
    // fixed -- so the working copy is the value array alone, not the matrix.
    std::vector<Complex> working = A.values();

    for (int iter = 0; iter < iterations; ++iter) {
        std::vector<double> row_max(static_cast<std::size_t>(n), 0.0);
        for (int i = 0; i < n; ++i) {
            for (int k = row_ptr[static_cast<std::size_t>(i)]; k < row_ptr[static_cast<std::size_t>(i) + 1]; ++k) {
                const double mag = std::abs(working[static_cast<std::size_t>(k)]);
                if (mag > row_max[static_cast<std::size_t>(i)]) row_max[static_cast<std::size_t>(i)] = mag;
            }
        }

        std::vector<double> factor(static_cast<std::size_t>(n), 1.0);
        for (int i = 0; i < n; ++i) {
            // A zero row (or a row that has already converged to ~0, which
            // shouldn't happen for a well-posed system) is left unscaled rather
            // than dividing by a vanishing number. A structurally empty row --
            // which CSR represents as an empty segment, so the loop above never
            // runs for it -- lands here too, with row_max still 0.
            factor[static_cast<std::size_t>(i)] =
                (row_max[static_cast<std::size_t>(i)] > 0.0) ? 1.0 / std::sqrt(row_max[static_cast<std::size_t>(i)]) : 1.0;
            d[static_cast<std::size_t>(i)] *= factor[static_cast<std::size_t>(i)];
        }

        for (int i = 0; i < n; ++i) {
            for (int k = row_ptr[static_cast<std::size_t>(i)]; k < row_ptr[static_cast<std::size_t>(i) + 1]; ++k) {
                working[static_cast<std::size_t>(k)] *=
                    factor[static_cast<std::size_t>(i)] * factor[static_cast<std::size_t>(col_index[static_cast<std::size_t>(k)])];
            }
        }
    }

    return d;
}

SparseMatrixZ apply_symmetric_equilibration(const SparseMatrixZ& A, const std::vector<double>& d) {
    const int n = A.rows();
    if (A.cols() != n || static_cast<int>(d.size()) != n) {
        throw std::invalid_argument("apply_symmetric_equilibration: dimension mismatch");
    }
    SparseMatrixZ out = A;
    const auto& row_ptr = out.row_ptr();
    const auto& col_index = out.col_index();
    std::vector<Complex>& values = out.mutable_values();
    for (int i = 0; i < n; ++i) {
        for (int k = row_ptr[static_cast<std::size_t>(i)]; k < row_ptr[static_cast<std::size_t>(i) + 1]; ++k) {
            values[static_cast<std::size_t>(k)] *=
                Complex(d[static_cast<std::size_t>(i)] * d[static_cast<std::size_t>(col_index[static_cast<std::size_t>(k)])], 0.0);
        }
    }
    return out;
}

std::vector<Complex> scale_rhs(const std::vector<Complex>& b, const std::vector<double>& d) {
    if (d.size() != b.size()) {
        throw std::invalid_argument("scale_rhs: dimension mismatch");
    }
    std::vector<Complex> out(b.size());
    for (std::size_t i = 0; i < b.size(); ++i) {
        out[i] = b[i] * Complex(d[i], 0.0);
    }
    return out;
}

std::vector<Complex> recover_equilibrated_solution(const std::vector<Complex>& y, const std::vector<double>& d) {
    // Same operation as scale_rhs (elementwise diag(d) multiply) -- named
    // separately because it means something different at the call site
    // (undoing the substitution x = diag(d) y, not scaling a right-hand side).
    return scale_rhs(y, d);
}

}  // namespace aphi_solver
