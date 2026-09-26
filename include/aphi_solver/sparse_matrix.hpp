#pragma once

#include <algorithm>
#include <complex>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

namespace aphi_solver {

/// A sparse matrix with a two-phase lifecycle: triplet accumulation during
/// assembly, then a single `compress()` into CSR for everything afterwards.
///
/// This replaces the COO/`std::map`-based `SparseMatrix` that
/// `incidence.hpp` carried through Phases 02-03, per
/// `docs/ENGINEERING_STANDARDS.md` item 3 ("compressed sparse formats
/// (CSR/CSC) ... rather than the COO/`std::map`-based `SparseMatrix`") and
/// `docs/ROADMAP.md` Phase 03.5, step 2. Two things forced the change at
/// this point in the project:
///
///  - **Complex.** The frequency-domain A-Phi blocks are complex
///    (`docs/FORMULATION.md` Sec 5.2), but the old type was real-only while
///    the only complex type available (`ComplexMatrix`) was *dense*. A
///    dense K_AA is already ~55 MB at `meshes/cube_6.msh` and impossible at
///    real EDA mesh sizes. Templating on the scalar covers both from one
///    implementation rather than a second complex-typed copy of everything.
///  - **Speed.** The old `coalesce()` summed duplicates through a
///    `std::map<pair<int,int>, double>` -- a node-based container with a
///    cache-hostile access pattern, in the one code path that runs once per
///    assembled entry. `compress()` below uses a counting sort into CSR
///    instead (O(nnz + rows), plus a per-row sort over segments that are
///    ~15-50 entries wide for a tetrahedral FEM matrix).
///
/// Scalar type `T` is `double` (the real incidence operators C and G) or
/// `std::complex<double>` (the assembled system blocks). Header-only because
/// it is a template; the project's usual .hpp/.cpp split does not apply
/// without explicit instantiation, which would defeat the point of having
/// one implementation serve both scalars.
template <typename T>
class Sparse {
public:
    struct Triplet {
        int row;
        int col;
        T value;
    };

    Sparse() = default;
    Sparse(int rows, int cols) : rows_(rows), cols_(cols) {}

    int rows() const { return rows_; }
    int cols() const { return cols_; }

    /// True once `compress()` has run. The CSR accessors and every operation
    /// below require this; `add()` requires the opposite. The two phases are
    /// deliberately not interchangeable -- silently re-expanding a
    /// compressed matrix so one more entry can be appended is how an
    /// assembly loop ends up quadratic without anyone noticing.
    bool compressed() const { return compressed_; }

    std::size_t nnz() const { return compressed_ ? values_.size() : triplets_.size(); }

    // ---------------------------------------------------------------- build

    /// Accumulates one contribution. Duplicates are expected and summed by
    /// `compress()` -- that is exactly what element-matrix scatter produces,
    /// since every interior mesh entity is shared by several tets.
    void add(int row, int col, const T& value) {
        if (compressed_) {
            throw std::logic_error("Sparse::add called after compress()");
        }
        check_index(row, col);
        triplets_.push_back({row, col, value});
    }

    void reserve(std::size_t n) { triplets_.reserve(n); }

    /// Sorts into CSR and sums duplicate (row, col) pairs. Idempotent.
    void compress() {
        if (compressed_) return;

        // Counting sort by row: count, prefix-sum into row_ptr_, then
        // scatter. Avoids sorting the whole triplet array by (row, col),
        // which would be O(nnz log nnz) over the full width instead of
        // O(nnz) plus a sort confined to each row's own short segment.
        row_ptr_.assign(static_cast<std::size_t>(rows_) + 1, 0);
        for (const Triplet& t : triplets_) {
            ++row_ptr_[static_cast<std::size_t>(t.row) + 1];
        }
        for (int r = 0; r < rows_; ++r) {
            row_ptr_[static_cast<std::size_t>(r) + 1] += row_ptr_[static_cast<std::size_t>(r)];
        }

        std::vector<int> scattered_col(triplets_.size());
        std::vector<T> scattered_val(triplets_.size());
        std::vector<int> cursor(row_ptr_.begin(), row_ptr_.end() - 1);
        for (const Triplet& t : triplets_) {
            const int dst = cursor[static_cast<std::size_t>(t.row)]++;
            scattered_col[static_cast<std::size_t>(dst)] = t.col;
            scattered_val[static_cast<std::size_t>(dst)] = t.value;
        }

        // Sort each row's segment by column and collapse duplicates,
        // compacting in place as we go.
        col_index_.clear();
        values_.clear();
        col_index_.reserve(scattered_col.size());
        values_.reserve(scattered_val.size());

        std::vector<int> order;
        std::vector<int> new_row_ptr(static_cast<std::size_t>(rows_) + 1, 0);
        for (int r = 0; r < rows_; ++r) {
            const int begin = row_ptr_[static_cast<std::size_t>(r)];
            const int end = row_ptr_[static_cast<std::size_t>(r) + 1];

            order.clear();
            for (int k = begin; k < end; ++k) order.push_back(k);
            std::sort(order.begin(), order.end(), [&](int a, int b) {
                return scattered_col[static_cast<std::size_t>(a)] < scattered_col[static_cast<std::size_t>(b)];
            });

            int prev_col = -1;
            for (int k : order) {
                const int c = scattered_col[static_cast<std::size_t>(k)];
                const T& v = scattered_val[static_cast<std::size_t>(k)];
                if (c == prev_col) {
                    values_.back() += v;
                } else {
                    col_index_.push_back(c);
                    values_.push_back(v);
                    prev_col = c;
                }
            }
            new_row_ptr[static_cast<std::size_t>(r) + 1] = static_cast<int>(col_index_.size());
        }

        row_ptr_.swap(new_row_ptr);
        triplets_.clear();
        triplets_.shrink_to_fit();
        compressed_ = true;
    }

    /// Builds an already-compressed matrix with a given structure and all
    /// values zero, ready to be filled in place.
    ///
    /// This is the assembly path, and it is the opposite way round from
    /// `add()` + `compress()`: the structure is known first (from
    /// `build_sparsity`, which reads the DOF map) and only the numbers are
    /// missing. That buys two things the triplet path cannot give. It never
    /// holds every contribution in memory at once -- on the cylinder that
    /// would be ~2.2 M triplets against 941 665 nonzeros -- and it is
    /// **reusable**: the structure depends on connectivity while the values
    /// depend on materials and frequency, so a 41-point sweep builds this
    /// once and refills `mutable_values()` 41 times.
    ///
    /// `row_ptr` must have `rows + 1` entries and be non-decreasing, and
    /// each row's slice of `col_index` must be strictly ascending and in
    /// range -- which is what `find_slot` and the CSR operations here
    /// assume. Throws std::invalid_argument otherwise, because a malformed
    /// structure would otherwise surface as a wrong answer rather than a
    /// failure.
    static Sparse<T> from_pattern(int rows, int cols, std::vector<int> row_ptr,
                                  std::vector<int> col_index) {
        if (static_cast<int>(row_ptr.size()) != rows + 1) {
            throw std::invalid_argument("Sparse::from_pattern: row_ptr must have rows + 1 entries");
        }
        if (row_ptr.front() != 0 || row_ptr.back() != static_cast<int>(col_index.size())) {
            throw std::invalid_argument(
                "Sparse::from_pattern: row_ptr must run from 0 to col_index.size()");
        }
        for (int r = 0; r < rows; ++r) {
            const int begin = row_ptr[static_cast<std::size_t>(r)];
            const int end = row_ptr[static_cast<std::size_t>(r) + 1];
            if (end < begin) {
                throw std::invalid_argument("Sparse::from_pattern: row_ptr is not non-decreasing");
            }
            for (int k = begin; k < end; ++k) {
                const int c = col_index[static_cast<std::size_t>(k)];
                if (c < 0 || c >= cols) {
                    throw std::invalid_argument("Sparse::from_pattern: a column index is out of range");
                }
                if (k > begin && c <= col_index[static_cast<std::size_t>(k - 1)]) {
                    throw std::invalid_argument(
                        "Sparse::from_pattern: each row's columns must be strictly ascending");
                }
            }
        }

        Sparse<T> m(rows, cols);
        m.row_ptr_ = std::move(row_ptr);
        m.col_index_ = std::move(col_index);
        m.values_.assign(m.col_index_.size(), T{});
        m.compressed_ = true;
        return m;
    }

    // ------------------------------------------------------------ CSR access

    const std::vector<int>& row_ptr() const {
        require_compressed("row_ptr");
        return row_ptr_;
    }
    const std::vector<int>& col_index() const {
        require_compressed("col_index");
        return col_index_;
    }
    const std::vector<T>& values() const {
        require_compressed("values");
        return values_;
    }

    /// Mutable value access for in-place scaling (equilibration, the
    /// frequency-scaling transforms in `conditioning.hpp`). The sparsity
    /// pattern stays fixed, so only the values are exposed.
    std::vector<T>& mutable_values() {
        require_compressed("mutable_values");
        return values_;
    }

    // ------------------------------------------------------------ operations

    /// y = A * x
    std::vector<T> matvec(const std::vector<T>& x) const {
        require_compressed("matvec");
        if (static_cast<int>(x.size()) != cols_) {
            throw std::invalid_argument("Sparse::matvec: vector length does not match column count");
        }
        std::vector<T> y(static_cast<std::size_t>(rows_), T{});
        for (int r = 0; r < rows_; ++r) {
            T sum{};
            for (int k = row_ptr_[static_cast<std::size_t>(r)]; k < row_ptr_[static_cast<std::size_t>(r) + 1]; ++k) {
                sum += values_[static_cast<std::size_t>(k)] * x[static_cast<std::size_t>(col_index_[static_cast<std::size_t>(k)])];
            }
            y[static_cast<std::size_t>(r)] = sum;
        }
        return y;
    }

    /// y = A^T * x (the plain transpose, not the conjugate transpose -- the
    /// A-Phi system is complex *symmetric*, not Hermitian, per
    /// `docs/CONDITIONING.md`, so this is the transpose that matters here).
    /// Computed without forming A^T: each stored entry contributes to one
    /// output slot, so this is a scatter rather than a gather.
    std::vector<T> matvec_transpose(const std::vector<T>& x) const {
        require_compressed("matvec_transpose");
        if (static_cast<int>(x.size()) != rows_) {
            throw std::invalid_argument("Sparse::matvec_transpose: vector length does not match row count");
        }
        std::vector<T> y(static_cast<std::size_t>(cols_), T{});
        for (int r = 0; r < rows_; ++r) {
            const T& xr = x[static_cast<std::size_t>(r)];
            for (int k = row_ptr_[static_cast<std::size_t>(r)]; k < row_ptr_[static_cast<std::size_t>(r) + 1]; ++k) {
                y[static_cast<std::size_t>(col_index_[static_cast<std::size_t>(k)])] += values_[static_cast<std::size_t>(k)] * xr;
            }
        }
        return y;
    }

    Sparse<T> transposed() const {
        require_compressed("transposed");
        Sparse<T> t(cols_, rows_);
        t.reserve(values_.size());
        for (int r = 0; r < rows_; ++r) {
            for (int k = row_ptr_[static_cast<std::size_t>(r)]; k < row_ptr_[static_cast<std::size_t>(r) + 1]; ++k) {
                t.add(col_index_[static_cast<std::size_t>(k)], r, values_[static_cast<std::size_t>(k)]);
            }
        }
        t.compress();
        return t;
    }

    /// C = A * B, row-by-row through B's rows (Gustavson's method): for each
    /// nonzero A(i,k), accumulate A(i,k) * B(k,:) into a dense scratch row,
    /// then emit that row's nonzeros. O(flops) rather than the old
    /// implementation's per-entry `std::map` lookups.
    Sparse<T> multiply(const Sparse<T>& b) const {
        require_compressed("multiply");
        b.require_compressed("multiply (right operand)");
        if (cols_ != b.rows_) {
            throw std::invalid_argument("Sparse::multiply: inner dimensions do not match");
        }
        Sparse<T> c(rows_, b.cols_);
        std::vector<T> accum(static_cast<std::size_t>(b.cols_), T{});
        // Row-stamped marker rather than testing `accum[col] == 0` to decide
        // whether a column is new to this row: a running sum that passes
        // back through exactly zero would otherwise be pushed onto `touched`
        // a second time and redo work it has already done.
        std::vector<int> touched_in_row(static_cast<std::size_t>(b.cols_), -1);
        std::vector<int> touched;
        for (int r = 0; r < rows_; ++r) {
            touched.clear();
            for (int k = row_ptr_[static_cast<std::size_t>(r)]; k < row_ptr_[static_cast<std::size_t>(r) + 1]; ++k) {
                const int mid = col_index_[static_cast<std::size_t>(k)];
                const T& a_ik = values_[static_cast<std::size_t>(k)];
                for (int m = b.row_ptr_[static_cast<std::size_t>(mid)];
                     m < b.row_ptr_[static_cast<std::size_t>(mid) + 1]; ++m) {
                    const int col = b.col_index_[static_cast<std::size_t>(m)];
                    if (touched_in_row[static_cast<std::size_t>(col)] != r) {
                        touched_in_row[static_cast<std::size_t>(col)] = r;
                        touched.push_back(col);
                    }
                    accum[static_cast<std::size_t>(col)] += a_ik * b.values_[static_cast<std::size_t>(m)];
                }
            }
            std::sort(touched.begin(), touched.end());
            for (int col : touched) {
                // Entries that cancelled to exactly zero are dropped rather
                // than stored. They are not part of the product's structure,
                // and keeping them inflates nnz -- which matters here beyond
                // memory, because the tree-cotree fill-in comparison in
                // tools/compare_gauges.cpp reports nnz_D / nnz_A as a
                // headline number. Only exact zeros are pruned; a
                // near-cancellation is left alone rather than silently
                // swallowed by a tolerance.
                const T& v = accum[static_cast<std::size_t>(col)];
                if (!(v == T{})) c.add(r, col, v);
                accum[static_cast<std::size_t>(col)] = T{};
            }
        }
        c.compress();
        return c;
    }

    /// The principal submatrix on the kept indices: keeps row/column i iff
    /// `full_to_reduced[i] != -1`, placing it at `full_to_reduced[i]`.
    ///
    /// This is the whole of the Albanese-Rubinacci gauge reduction --
    /// "eliminating those rows and columns which correspond to the tree
    /// edges" (Munteanu, Sec. IV) is exactly a principal submatrix, which is
    /// why the method preserves the sparsity pattern with no fill-in. It
    /// lives here, on the matrix type, rather than in gauge_variants.cpp
    /// because the index map is the only thing that differs between the two
    /// callers: a bare edge-indexed curl-curl matrix (where the eliminated
    /// set is the tree edges) and a coupled [a; Phi] system (where it is the
    /// tree-edge A-DOFs, with every Phi row and column passing through
    /// untouched). Templated on the scalar, so the real and complex cases
    /// share one implementation -- see docs/ROADMAP.md Phase 03.5, step 4.
    Sparse<T> principal_submatrix(const std::vector<int>& full_to_reduced, int reduced_size) const {
        require_compressed("principal_submatrix");
        if (rows_ != cols_) {
            throw std::invalid_argument("Sparse::principal_submatrix: matrix must be square");
        }
        if (static_cast<int>(full_to_reduced.size()) != rows_) {
            throw std::invalid_argument("Sparse::principal_submatrix: index map length must match matrix size");
        }
        Sparse<T> out(reduced_size, reduced_size);
        for (int r = 0; r < rows_; ++r) {
            const int rr = full_to_reduced[static_cast<std::size_t>(r)];
            if (rr < 0) continue;  // eliminated row: skipped before its entries are read at all
            for (int k = row_ptr_[static_cast<std::size_t>(r)]; k < row_ptr_[static_cast<std::size_t>(r) + 1]; ++k) {
                const int cc = full_to_reduced[static_cast<std::size_t>(col_index_[static_cast<std::size_t>(k)])];
                if (cc < 0) continue;  // eliminated column
                out.add(rr, cc, values_[static_cast<std::size_t>(k)]);
            }
        }
        out.compress();
        return out;
    }

    /// A copy with every stored value multiplied by `s`. The sparsity
    /// pattern is unchanged, so this is O(nnz) with no structural work --
    /// which is what the frequency-scaling transforms in `conditioning.hpp`
    /// need (they rescale whole blocks by j*omega or 1/(j*omega)).
    Sparse<T> scaled(const T& s) const {
        require_compressed("scaled");
        Sparse<T> out(*this);
        for (T& v : out.values_) v *= s;
        return out;
    }

    /// True iff every stored entry has magnitude <= tol. `std::abs` covers
    /// both scalar types (absolute value for double, modulus for complex).
    bool is_zero(double tol = 1e-9) const {
        require_compressed("is_zero");
        for (const T& v : values_) {
            if (std::abs(v) > tol) return false;
        }
        return true;
    }

    // -------------------------------------------------------------- interop

    /// Exports the triplet (coordinate) arrays a sparse direct solver wants.
    /// MUMPS's assembled centralized format takes IRN/JCN/A with **1-based**
    /// Fortran indexing, which is why `index_base` defaults to 1 rather than
    /// 0 -- see `docs/LINEAR_SOLVER.md`. Keeping this alongside CSR is not
    /// redundancy: CSR is what matvec and the gauge reduction need, triplets
    /// are what the solver's own interface needs.
    void to_triplets(std::vector<int>& irn, std::vector<int>& jcn, std::vector<T>& a,
                      int index_base = 1) const {
        require_compressed("to_triplets");
        irn.clear();
        jcn.clear();
        a.clear();
        irn.reserve(values_.size());
        jcn.reserve(values_.size());
        a.reserve(values_.size());
        for (int r = 0; r < rows_; ++r) {
            for (int k = row_ptr_[static_cast<std::size_t>(r)]; k < row_ptr_[static_cast<std::size_t>(r) + 1]; ++k) {
                irn.push_back(r + index_base);
                jcn.push_back(col_index_[static_cast<std::size_t>(k)] + index_base);
                a.push_back(values_[static_cast<std::size_t>(k)]);
            }
        }
    }

    /// Dense conversion, for small test meshes and debugging only.
    std::vector<std::vector<T>> to_dense() const {
        require_compressed("to_dense");
        std::vector<std::vector<T>> d(static_cast<std::size_t>(rows_),
                                       std::vector<T>(static_cast<std::size_t>(cols_), T{}));
        for (int r = 0; r < rows_; ++r) {
            for (int k = row_ptr_[static_cast<std::size_t>(r)]; k < row_ptr_[static_cast<std::size_t>(r) + 1]; ++k) {
                d[static_cast<std::size_t>(r)][static_cast<std::size_t>(col_index_[static_cast<std::size_t>(k)])] =
                    values_[static_cast<std::size_t>(k)];
            }
        }
        return d;
    }

    /// Single-entry lookup, O(log nnz_in_row) via binary search on the row's
    /// sorted column segment. For tests and assertions -- an assembly or
    /// solver loop should walk the CSR arrays directly.
    T at(int row, int col) const {
        require_compressed("at");
        check_index(row, col);
        const int begin = row_ptr_[static_cast<std::size_t>(row)];
        const int end = row_ptr_[static_cast<std::size_t>(row) + 1];
        const auto first = col_index_.begin() + begin;
        const auto last = col_index_.begin() + end;
        const auto it = std::lower_bound(first, last, col);
        if (it == last || *it != col) return T{};
        return values_[static_cast<std::size_t>(it - col_index_.begin())];
    }

private:
    void check_index(int row, int col) const {
        if (row < 0 || row >= rows_ || col < 0 || col >= cols_) {
            throw std::out_of_range("Sparse: index (" + std::to_string(row) + "," + std::to_string(col) +
                                     ") outside " + std::to_string(rows_) + "x" + std::to_string(cols_));
        }
    }

    void require_compressed(const char* what) const {
        if (!compressed_) {
            throw std::logic_error(std::string("Sparse::") + what + " requires compress() first");
        }
    }

    int rows_ = 0;
    int cols_ = 0;
    bool compressed_ = false;

    std::vector<Triplet> triplets_;  // build phase only

    std::vector<int> row_ptr_;    // size rows_ + 1
    std::vector<int> col_index_;  // size nnz
    std::vector<T> values_;       // size nnz
};

using SparseMatrixD = Sparse<double>;
using SparseMatrixZ = Sparse<std::complex<double>>;

}  // namespace aphi_solver
