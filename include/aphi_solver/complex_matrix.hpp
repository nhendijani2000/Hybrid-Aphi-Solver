#pragma once

#include <complex>
#include <cstddef>
#include <stdexcept>
#include <vector>

#include "aphi_solver/sparse_matrix.hpp"

namespace aphi_solver {

using Complex = std::complex<double>;

/// A small, dependency-free dense complex matrix (row-major). Deliberately minimal:
/// this project's real assembly will eventually need a sparse solver for
/// production-size meshes (see Phase 05 of the roadmap), but the conditioning
/// module and its tests only need to operate on small block systems, so a
/// hand-written dense type keeps the build free of any third-party dependency
/// while that's true. A vector is just an n x 1 matrix here.
class ComplexMatrix {
public:
    ComplexMatrix() : rows_(0), cols_(0) {}
    ComplexMatrix(int rows, int cols) : rows_(rows), cols_(cols), data_(static_cast<size_t>(rows) * cols, Complex(0.0, 0.0)) {}

    static ComplexMatrix zero(int rows, int cols) { return ComplexMatrix(rows, cols); }
    static ComplexMatrix identity(int n) {
        ComplexMatrix m(n, n);
        for (int i = 0; i < n; ++i) m(i, i) = Complex(1.0, 0.0);
        return m;
    }

    int rows() const { return rows_; }
    int cols() const { return cols_; }

    Complex& operator()(int r, int c) {
        check_bounds(r, c);
        return data_[static_cast<size_t>(r) * cols_ + c];
    }
    const Complex& operator()(int r, int c) const {
        check_bounds(r, c);
        return data_[static_cast<size_t>(r) * cols_ + c];
    }

    ComplexMatrix operator+(const ComplexMatrix& other) const;
    ComplexMatrix operator-(const ComplexMatrix& other) const;
    ComplexMatrix operator*(Complex scalar) const;
    ComplexMatrix operator/(Complex scalar) const;
    ComplexMatrix operator*(const ComplexMatrix& other) const;  // matrix product

    /// Conjugate transpose (Hermitian adjoint).
    ComplexMatrix hermitian() const;

    /// Copies `src` into this matrix starting at (row_offset, col_offset).
    void set_block(int row_offset, int col_offset, const ComplexMatrix& src);

    /// Returns the sub-block [row_offset, row_offset+rows) x [col_offset, col_offset+cols).
    ComplexMatrix block(int row_offset, int col_offset, int rows, int cols) const;

    double frobenius_norm() const;

private:
    void check_bounds(int r, int c) const {
        if (r < 0 || r >= rows_ || c < 0 || c >= cols_) {
            throw std::out_of_range("ComplexMatrix index out of range");
        }
    }

    int rows_;
    int cols_;
    std::vector<Complex> data_;
};

/// Bridges from the sparse production type to this dense one. Since Phase
/// 03.5 the assembled system is sparse (`sparse_matrix.hpp`), while this
/// dense type remains the small-system ground truth that `solve_dense`
/// below operates on -- these two conversions are the seam between them,
/// and like `solve_dense` itself they are for verification-scale problems
/// only, never a production path.
ComplexMatrix to_dense_matrix(const SparseMatrixZ& m);

/// Wraps a dense vector as an n x 1 ComplexMatrix, for handing a
/// right-hand side to `solve_dense`.
ComplexMatrix to_column(const std::vector<Complex>& v);

/// Reads an n x 1 ComplexMatrix back out as a plain vector.
std::vector<Complex> from_column(const ComplexMatrix& m);

/// Solves A X = B for X via dense Gaussian elimination with partial pivoting.
/// A must be square and non-singular (throws std::runtime_error if a pivot is
/// numerically zero). Intended for the small systems this module and its tests
/// use -- swap for a real sparse solver (Phase 05) once problems grow past a few
/// hundred unknowns.
ComplexMatrix solve_dense(const ComplexMatrix& A, const ComplexMatrix& B);

/// Estimate of the 2-norm condition number (largest / smallest singular value) of
/// `A`, obtained by power iteration and inverse power iteration on A^H A. This is
/// a diagnostic for small test/benchmark systems, not a production-scale routine.
double estimate_condition_number(const ComplexMatrix& A, int iterations = 100);

}  // namespace aphi_solver
