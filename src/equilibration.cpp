#include "aphi_solver/equilibration.hpp"

#include <cmath>
#include <stdexcept>

namespace aphi_solver {

std::vector<double> compute_symmetric_equilibration(const ComplexMatrix& A, int iterations) {
    const int n = A.rows();
    if (A.cols() != n) {
        throw std::invalid_argument("compute_symmetric_equilibration: A must be square");
    }

    std::vector<double> d(n, 1.0);
    ComplexMatrix working = A;  // progressively rescaled working copy

    for (int iter = 0; iter < iterations; ++iter) {
        std::vector<double> row_max(n, 0.0);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                double mag = std::abs(working(i, j));
                if (mag > row_max[i]) row_max[i] = mag;
            }
        }

        std::vector<double> factor(n, 1.0);
        for (int i = 0; i < n; ++i) {
            // A zero row (or a row that has already converged to ~0, which
            // shouldn't happen for a well-posed system) is left unscaled rather
            // than dividing by a vanishing number.
            factor[i] = (row_max[i] > 0.0) ? 1.0 / std::sqrt(row_max[i]) : 1.0;
            d[i] *= factor[i];
        }

        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                working(i, j) *= factor[i] * factor[j];
            }
        }
    }

    return d;
}

ComplexMatrix apply_symmetric_equilibration(const ComplexMatrix& A, const std::vector<double>& d) {
    const int n = A.rows();
    if (A.cols() != n || static_cast<int>(d.size()) != n) {
        throw std::invalid_argument("apply_symmetric_equilibration: dimension mismatch");
    }
    ComplexMatrix out(n, n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            out(i, j) = A(i, j) * Complex(d[i] * d[j], 0.0);
        }
    }
    return out;
}

ComplexMatrix scale_rhs(const ComplexMatrix& b, const std::vector<double>& d) {
    if (static_cast<int>(d.size()) != b.rows()) {
        throw std::invalid_argument("scale_rhs: dimension mismatch");
    }
    ComplexMatrix out(b.rows(), b.cols());
    for (int i = 0; i < b.rows(); ++i) {
        for (int j = 0; j < b.cols(); ++j) {
            out(i, j) = b(i, j) * Complex(d[i], 0.0);
        }
    }
    return out;
}

ComplexMatrix recover_equilibrated_solution(const ComplexMatrix& y, const std::vector<double>& d) {
    // Same operation as scale_rhs (elementwise diag(d) multiply) -- named
    // separately because it means something different at the call site
    // (undoing the substitution x = diag(d) y, not scaling a right-hand side).
    return scale_rhs(y, d);
}

}  // namespace aphi_solver
