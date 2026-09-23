#include "aphi_solver/complex_matrix.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace aphi_solver {

ComplexMatrix ComplexMatrix::operator+(const ComplexMatrix& other) const {
    if (rows_ != other.rows_ || cols_ != other.cols_) {
        throw std::invalid_argument("ComplexMatrix::operator+: shape mismatch");
    }
    ComplexMatrix out(rows_, cols_);
    for (size_t i = 0; i < data_.size(); ++i) out.data_[i] = data_[i] + other.data_[i];
    return out;
}

ComplexMatrix ComplexMatrix::operator-(const ComplexMatrix& other) const {
    if (rows_ != other.rows_ || cols_ != other.cols_) {
        throw std::invalid_argument("ComplexMatrix::operator-: shape mismatch");
    }
    ComplexMatrix out(rows_, cols_);
    for (size_t i = 0; i < data_.size(); ++i) out.data_[i] = data_[i] - other.data_[i];
    return out;
}

ComplexMatrix ComplexMatrix::operator*(Complex scalar) const {
    ComplexMatrix out(rows_, cols_);
    for (size_t i = 0; i < data_.size(); ++i) out.data_[i] = data_[i] * scalar;
    return out;
}

ComplexMatrix ComplexMatrix::operator/(Complex scalar) const {
    if (scalar == Complex(0.0, 0.0)) {
        throw std::invalid_argument("ComplexMatrix::operator/: division by zero");
    }
    ComplexMatrix out(rows_, cols_);
    for (size_t i = 0; i < data_.size(); ++i) out.data_[i] = data_[i] / scalar;
    return out;
}

ComplexMatrix ComplexMatrix::operator*(const ComplexMatrix& other) const {
    if (cols_ != other.rows_) {
        throw std::invalid_argument("ComplexMatrix::operator*: inner dimensions mismatch");
    }
    ComplexMatrix out(rows_, other.cols_);
    for (int i = 0; i < rows_; ++i) {
        for (int k = 0; k < cols_; ++k) {
            const Complex a_ik = (*this)(i, k);
            if (a_ik == Complex(0.0, 0.0)) continue;
            for (int j = 0; j < other.cols_; ++j) {
                out(i, j) += a_ik * other(k, j);
            }
        }
    }
    return out;
}

ComplexMatrix ComplexMatrix::hermitian() const {
    ComplexMatrix out(cols_, rows_);
    for (int i = 0; i < rows_; ++i) {
        for (int j = 0; j < cols_; ++j) {
            out(j, i) = std::conj((*this)(i, j));
        }
    }
    return out;
}

void ComplexMatrix::set_block(int row_offset, int col_offset, const ComplexMatrix& src) {
    if (row_offset + src.rows() > rows_ || col_offset + src.cols() > cols_) {
        throw std::out_of_range("ComplexMatrix::set_block: block does not fit");
    }
    for (int i = 0; i < src.rows(); ++i) {
        for (int j = 0; j < src.cols(); ++j) {
            (*this)(row_offset + i, col_offset + j) = src(i, j);
        }
    }
}

ComplexMatrix ComplexMatrix::block(int row_offset, int col_offset, int rows, int cols) const {
    if (row_offset + rows > rows_ || col_offset + cols > cols_) {
        throw std::out_of_range("ComplexMatrix::block: out of range");
    }
    ComplexMatrix out(rows, cols);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            out(i, j) = (*this)(row_offset + i, col_offset + j);
        }
    }
    return out;
}

double ComplexMatrix::frobenius_norm() const {
    double sum = 0.0;
    for (const auto& v : data_) sum += std::norm(v);
    return std::sqrt(sum);
}

ComplexMatrix to_dense_matrix(const SparseMatrixZ& m) {
    ComplexMatrix out(m.rows(), m.cols());
    const auto& row_ptr = m.row_ptr();
    const auto& col_index = m.col_index();
    const auto& values = m.values();
    for (int r = 0; r < m.rows(); ++r) {
        for (int k = row_ptr[static_cast<std::size_t>(r)]; k < row_ptr[static_cast<std::size_t>(r) + 1]; ++k) {
            out(r, col_index[static_cast<std::size_t>(k)]) = values[static_cast<std::size_t>(k)];
        }
    }
    return out;
}

ComplexMatrix to_column(const std::vector<Complex>& v) {
    ComplexMatrix out(static_cast<int>(v.size()), 1);
    for (std::size_t i = 0; i < v.size(); ++i) {
        out(static_cast<int>(i), 0) = v[i];
    }
    return out;
}

std::vector<Complex> from_column(const ComplexMatrix& m) {
    if (m.cols() != 1) {
        throw std::invalid_argument("from_column: matrix must have exactly one column");
    }
    std::vector<Complex> out(static_cast<std::size_t>(m.rows()));
    for (int i = 0; i < m.rows(); ++i) {
        out[static_cast<std::size_t>(i)] = m(i, 0);
    }
    return out;
}

ComplexMatrix solve_dense(const ComplexMatrix& A, const ComplexMatrix& B) {
    const int n = A.rows();
    if (A.cols() != n) throw std::invalid_argument("solve_dense: A must be square");
    if (B.rows() != n) throw std::invalid_argument("solve_dense: B row count must match A");

    // Augmented working copy for Gaussian elimination with partial pivoting.
    ComplexMatrix M = A;
    ComplexMatrix X = B;
    const int k = B.cols();

    for (int col = 0; col < n; ++col) {
        // Partial pivot: largest magnitude entry in this column, at or below the diagonal.
        int pivot_row = col;
        double best = std::abs(M(col, col));
        for (int r = col + 1; r < n; ++r) {
            double mag = std::abs(M(r, col));
            if (mag > best) {
                best = mag;
                pivot_row = r;
            }
        }
        if (best == 0.0) {
            throw std::runtime_error("solve_dense: matrix is singular (zero pivot)");
        }
        if (pivot_row != col) {
            for (int j = 0; j < n; ++j) std::swap(M(col, j), M(pivot_row, j));
            for (int j = 0; j < k; ++j) std::swap(X(col, j), X(pivot_row, j));
        }

        const Complex pivot = M(col, col);
        for (int r = col + 1; r < n; ++r) {
            const Complex factor = M(r, col) / pivot;
            if (factor == Complex(0.0, 0.0)) continue;
            for (int j = col; j < n; ++j) M(r, j) -= factor * M(col, j);
            for (int j = 0; j < k; ++j) X(r, j) -= factor * X(col, j);
        }
    }

    // Back-substitution.
    ComplexMatrix result(n, k);
    for (int row = n - 1; row >= 0; --row) {
        for (int j = 0; j < k; ++j) {
            Complex sum = X(row, j);
            for (int c = row + 1; c < n; ++c) sum -= M(row, c) * result(c, j);
            result(row, j) = sum / M(row, row);
        }
    }
    return result;
}

namespace {

// Rayleigh quotient v^H M v / v^H v for a candidate vector v (n x 1).
double rayleigh_quotient(const ComplexMatrix& M, const ComplexMatrix& v) {
    ComplexMatrix Mv = M * v;
    Complex numerator = (v.hermitian() * Mv)(0, 0);
    Complex denominator = (v.hermitian() * v)(0, 0);
    return numerator.real() / denominator.real();
}

ComplexMatrix normalize(const ComplexMatrix& v) {
    double norm = v.frobenius_norm();
    if (norm == 0.0) throw std::runtime_error("normalize: zero vector");
    return v / Complex(norm, 0.0);
}

}  // namespace

// True once the Rayleigh quotient has settled: the same relative test the
// sparse overload in gauge_variants.cpp uses, kept identical on purpose so
// both estimators stop at the same point for the same matrix. See the
// shared-contract note in complex_matrix.hpp.
static bool converged(double next, double previous, double tol) {
    return std::abs(next - previous) < tol * std::max(1.0, std::abs(next));
}

double estimate_condition_number(const ComplexMatrix& A, int max_iterations, double tol) {
    const int n = A.rows();
    if (A.cols() != n) {
        throw std::invalid_argument("estimate_condition_number: A must be square");
    }
    if (n <= 1) return 1.0;
    const ComplexMatrix Ah = A.hermitian();
    const ComplexMatrix M = Ah * A;  // Hermitian positive semi-definite; eigenvalues = singular values^2 of A

    // Power iteration for the largest eigenvalue of M.
    ComplexMatrix v_max(n, 1);
    for (int i = 0; i < n; ++i) v_max(i, 0) = Complex(1.0, 0.0);
    v_max = normalize(v_max);
    double lambda_max = 0.0;
    for (int it = 0; it < max_iterations; ++it) {
        v_max = normalize(M * v_max);
        const double next = rayleigh_quotient(M, v_max);
        const bool done = converged(next, lambda_max, tol);
        lambda_max = next;
        if (done) break;
    }

    // Inverse power iteration for the smallest eigenvalue of M.
    ComplexMatrix v_min(n, 1);
    for (int i = 0; i < n; ++i) v_min(i, 0) = Complex(1.0, 0.0);
    v_min = normalize(v_min);
    double lambda_min = 0.0;
    for (int it = 0; it < max_iterations; ++it) {
        v_min = normalize(solve_dense(M, v_min));
        const double next = rayleigh_quotient(M, v_min);
        const bool done = converged(next, lambda_min, tol);
        lambda_min = next;
        if (done) break;
    }

    if (lambda_min <= 0.0) {
        return std::numeric_limits<double>::infinity();
    }
    return std::sqrt(lambda_max / lambda_min);
}

}  // namespace aphi_solver
