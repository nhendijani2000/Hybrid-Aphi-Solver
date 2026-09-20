// Minimal, dependency-free test runner -- see the comment in
// tests/CMakeLists.txt for why no external framework is fetched here.
//
// Checks for the templated CSR sparse type (docs/ROADMAP.md Phase 03.5,
// step 2). Exercised at both scalar types it exists to serve: double (the
// real incidence operators C and G) and std::complex<double> (the assembled
// frequency-domain A-Phi blocks).

#include <complex>
#include <iostream>
#include <string>
#include <vector>

#include "aphi_solver/sparse_matrix.hpp"

using aphi_solver::Sparse;
using aphi_solver::SparseMatrixD;
using aphi_solver::SparseMatrixZ;

using Complex = std::complex<double>;

namespace {

int g_checks = 0;
int g_failures = 0;

void check(bool condition, const std::string& what) {
    ++g_checks;
    if (!condition) {
        ++g_failures;
        std::cerr << "FAIL: " << what << "\n";
    }
}

bool nearly(double a, double b, double tol = 1e-12) { return std::abs(a - b) < tol; }
bool nearly(const Complex& a, const Complex& b, double tol = 1e-12) { return std::abs(a - b) < tol; }

// A = [[1, 0, 2],
//      [0, 3, 0],
//      [4, 0, 5]]
// Built with entries deliberately out of column order and with (0,0) split
// across two contributions, mimicking what element-matrix scatter produces:
// every interior mesh entity is shared by several tets, so duplicate
// (row, col) pairs are the normal case, not an edge case.
SparseMatrixD make_test_matrix() {
    SparseMatrixD a(3, 3);
    a.add(0, 2, 2.0);
    a.add(2, 2, 5.0);
    a.add(0, 0, 0.5);
    a.add(1, 1, 3.0);
    a.add(2, 0, 4.0);
    a.add(0, 0, 0.5);  // duplicate, must sum with the one above
    a.compress();
    return a;
}

void test_build_and_compress() {
    const SparseMatrixD a = make_test_matrix();

    check(a.rows() == 3 && a.cols() == 3, "compress: dimensions preserved");
    check(a.compressed(), "compress: reports compressed");
    check(a.nnz() == 5, "compress: duplicate (0,0) summed into a single stored entry (nnz == 5)");

    check(nearly(a.at(0, 0), 1.0), "compress: duplicates at (0,0) summed to 1.0");
    check(nearly(a.at(0, 2), 2.0), "compress: (0,2) == 2");
    check(nearly(a.at(1, 1), 3.0), "compress: (1,1) == 3");
    check(nearly(a.at(2, 0), 4.0), "compress: (2,0) == 4");
    check(nearly(a.at(2, 2), 5.0), "compress: (2,2) == 5");
    check(nearly(a.at(0, 1), 0.0), "compress: structural zero (0,1) reads as 0");
    check(nearly(a.at(1, 0), 0.0), "compress: structural zero (1,0) reads as 0");

    // CSR invariants: row_ptr is non-decreasing, spans [0, nnz], and each
    // row's column segment is sorted ascending (which `at`'s binary search
    // and the row-merge in `multiply` both depend on).
    const auto& rp = a.row_ptr();
    check(rp.size() == 4, "CSR: row_ptr has rows + 1 entries");
    check(rp.front() == 0 && rp.back() == static_cast<int>(a.nnz()), "CSR: row_ptr spans [0, nnz]");
    bool monotone = true, sorted_rows = true;
    for (int r = 0; r < a.rows(); ++r) {
        if (rp[static_cast<std::size_t>(r)] > rp[static_cast<std::size_t>(r) + 1]) monotone = false;
        for (int k = rp[static_cast<std::size_t>(r)] + 1; k < rp[static_cast<std::size_t>(r) + 1]; ++k) {
            if (a.col_index()[static_cast<std::size_t>(k - 1)] >= a.col_index()[static_cast<std::size_t>(k)]) {
                sorted_rows = false;
            }
        }
    }
    check(monotone, "CSR: row_ptr is non-decreasing");
    check(sorted_rows, "CSR: each row's columns are strictly ascending");

    // compress() is idempotent.
    SparseMatrixD b = make_test_matrix();
    b.compress();
    check(b.nnz() == 5, "compress: idempotent (second call is a no-op)");
}

void test_matvec() {
    const SparseMatrixD a = make_test_matrix();

    const std::vector<double> x{1.0, 1.0, 1.0};
    const std::vector<double> y = a.matvec(x);
    check(y.size() == 3, "matvec: result length == rows");
    check(nearly(y[0], 3.0) && nearly(y[1], 3.0) && nearly(y[2], 9.0),
          "matvec: A * [1,1,1] == [3,3,9] (row sums)");

    const std::vector<double> yt = a.matvec_transpose(x);
    check(nearly(yt[0], 5.0) && nearly(yt[1], 3.0) && nearly(yt[2], 7.0),
          "matvec_transpose: A^T * [1,1,1] == [5,3,7] (column sums)");

    // A non-uniform vector, so a transposed index mix-up cannot pass by
    // symmetry of the input.
    const std::vector<double> v{1.0, 2.0, 3.0};
    const std::vector<double> yv = a.matvec(v);
    check(nearly(yv[0], 7.0) && nearly(yv[1], 6.0) && nearly(yv[2], 19.0),
          "matvec: A * [1,2,3] == [7,6,19]");
}

void test_transpose_and_multiply() {
    const SparseMatrixD a = make_test_matrix();

    const SparseMatrixD at = a.transposed();
    check(at.rows() == 3 && at.cols() == 3, "transposed: dimensions swapped");
    check(nearly(at.at(2, 0), 2.0), "transposed: (0,2) -> (2,0)");
    check(nearly(at.at(0, 2), 4.0), "transposed: (2,0) -> (0,2)");
    check(at.nnz() == a.nnz(), "transposed: nnz preserved");

    // A * A, computed by hand:
    //   [[9, 0, 12], [0, 9, 0], [24, 0, 33]]
    const SparseMatrixD aa = a.multiply(a);
    check(nearly(aa.at(0, 0), 9.0), "multiply: (A*A)(0,0) == 9");
    check(nearly(aa.at(0, 2), 12.0), "multiply: (A*A)(0,2) == 12");
    check(nearly(aa.at(1, 1), 9.0), "multiply: (A*A)(1,1) == 9");
    check(nearly(aa.at(2, 0), 24.0), "multiply: (A*A)(2,0) == 24");
    check(nearly(aa.at(2, 2), 33.0), "multiply: (A*A)(2,2) == 33");
    check(nearly(aa.at(0, 1), 0.0) && nearly(aa.at(1, 0), 0.0) && nearly(aa.at(2, 1), 0.0),
          "multiply: (A*A) zeros where expected");

    // Cross-check the product against matvec applied twice: (A*A)x == A(Ax)
    // for an arbitrary x -- independent of the hand-computed entries above.
    const std::vector<double> x{2.0, -1.0, 0.5};
    const std::vector<double> direct = aa.matvec(x);
    const std::vector<double> twice = a.matvec(a.matvec(x));
    check(nearly(direct[0], twice[0]) && nearly(direct[1], twice[1]) && nearly(direct[2], twice[2]),
          "multiply: (A*A)x agrees with A(Ax)");

    // An exactly-cancelling product must leave no entry above tolerance --
    // the same property the Phase 02 CG = 0 identity check relies on.
    SparseMatrixD p(2, 2);
    p.add(0, 0, 1.0);
    p.add(0, 1, 1.0);
    p.compress();
    SparseMatrixD q(2, 2);
    q.add(0, 0, 1.0);
    q.add(1, 0, -1.0);
    q.compress();
    const SparseMatrixD pq = p.multiply(q);
    check(pq.is_zero(), "multiply: exactly cancelling product is_zero()");
    check(pq.nnz() == 0,
          "multiply: exact cancellation is pruned, not stored as an explicit zero");
    check(!a.is_zero(), "is_zero: a nonzero matrix is not reported zero");

    // A column whose running sum passes back THROUGH zero before ending
    // nonzero must still come out with its final value, and exactly once.
    // Row 0 of L is [1,1,1]; column 0 of R is [1,-1,5], so the accumulator
    // goes 1 -> 0 -> 5. Testing "is this column new?" by checking the
    // accumulator against zero (rather than a row stamp) mistakes the
    // midpoint for a fresh column.
    SparseMatrixD L(1, 3);
    L.add(0, 0, 1.0);
    L.add(0, 1, 1.0);
    L.add(0, 2, 1.0);
    L.compress();
    SparseMatrixD R(3, 1);
    R.add(0, 0, 1.0);
    R.add(1, 0, -1.0);
    R.add(2, 0, 5.0);
    R.compress();
    const SparseMatrixD LR = L.multiply(R);
    check(LR.nnz() == 1, "multiply: a sum passing through zero yields exactly one stored entry");
    check(nearly(LR.at(0, 0), 5.0), "multiply: a sum passing through zero keeps its final value");
}

void test_triplet_export() {
    const SparseMatrixD a = make_test_matrix();

    std::vector<int> irn, jcn;
    std::vector<double> vals;

    // Default is 1-based: MUMPS's assembled centralized format uses Fortran
    // indexing (docs/LINEAR_SOLVER.md).
    a.to_triplets(irn, jcn, vals);
    check(irn.size() == a.nnz() && jcn.size() == a.nnz() && vals.size() == a.nnz(),
          "to_triplets: array lengths == nnz");
    bool one_based_ok = true;
    for (std::size_t k = 0; k < irn.size(); ++k) {
        if (irn[k] < 1 || irn[k] > 3 || jcn[k] < 1 || jcn[k] > 3) one_based_ok = false;
    }
    check(one_based_ok, "to_triplets: default export is 1-based (MUMPS convention)");

    a.to_triplets(irn, jcn, vals, 0);
    bool round_trips = true;
    for (std::size_t k = 0; k < irn.size(); ++k) {
        if (!nearly(a.at(irn[k], jcn[k]), vals[k])) round_trips = false;
    }
    check(round_trips, "to_triplets: 0-based export round-trips against at()");

    const auto dense = a.to_dense();
    check(dense.size() == 3 && dense[0].size() == 3, "to_dense: shape");
    check(nearly(dense[2][2], 5.0) && nearly(dense[0][1], 0.0), "to_dense: values and structural zeros");
}

void test_complex_instantiation() {
    // B = [[1+i, 2  ],
    //      [0,   3-i]]
    SparseMatrixZ b(2, 2);
    b.add(0, 0, Complex(1.0, 1.0));
    b.add(0, 1, Complex(2.0, 0.0));
    b.add(1, 1, Complex(3.0, -1.0));
    b.compress();

    check(b.nnz() == 3, "complex: nnz == 3");
    check(nearly(b.at(0, 0), Complex(1.0, 1.0)), "complex: at(0,0) == 1+i");

    // B * [1, i] == [1+3i, 1+3i]  (by hand: (1+i)+2i = 1+3i; (3-i)i = 1+3i)
    const std::vector<Complex> x{Complex(1.0, 0.0), Complex(0.0, 1.0)};
    const std::vector<Complex> y = b.matvec(x);
    check(nearly(y[0], Complex(1.0, 3.0)), "complex: matvec component 0 == 1+3i");
    check(nearly(y[1], Complex(1.0, 3.0)), "complex: matvec component 1 == 1+3i");

    // Plain (non-conjugate) transpose: the A-Phi system is complex
    // *symmetric*, not Hermitian (docs/CONDITIONING.md), so the entries must
    // move without conjugation.
    const SparseMatrixZ bt = b.transposed();
    check(nearly(bt.at(1, 0), Complex(2.0, 0.0)), "complex: transposed moves (0,1) -> (1,0)");
    check(nearly(bt.at(0, 0), Complex(1.0, 1.0)), "complex: transposed does NOT conjugate (1+i stays 1+i)");

    // Complex duplicate summation.
    SparseMatrixZ d(1, 1);
    d.add(0, 0, Complex(1.0, 2.0));
    d.add(0, 0, Complex(3.0, -5.0));
    d.compress();
    check(nearly(d.at(0, 0), Complex(4.0, -3.0)), "complex: duplicates sum to 4-3i");

    // is_zero uses std::abs, i.e. the modulus, for complex scalars.
    SparseMatrixZ z(1, 1);
    z.add(0, 0, Complex(1e-15, -1e-15));
    z.compress();
    check(z.is_zero(), "complex: is_zero uses modulus");
}

void test_phase_guards() {
    SparseMatrixD a(2, 2);
    a.add(0, 0, 1.0);

    bool threw = false;
    try {
        (void)a.matvec(std::vector<double>{1.0, 1.0});
    } catch (const std::logic_error&) {
        threw = true;
    }
    check(threw, "guard: matvec before compress() throws");

    a.compress();
    threw = false;
    try {
        a.add(1, 1, 1.0);
    } catch (const std::logic_error&) {
        threw = true;
    }
    check(threw, "guard: add after compress() throws");

    threw = false;
    try {
        SparseMatrixD b(2, 2);
        b.add(2, 0, 1.0);  // row out of range
    } catch (const std::out_of_range&) {
        threw = true;
    }
    check(threw, "guard: out-of-range index throws");

    threw = false;
    try {
        SparseMatrixD p(2, 3);
        p.compress();
        SparseMatrixD q(2, 2);
        q.compress();
        (void)p.multiply(q);  // 3 != 2
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    check(threw, "guard: multiply with mismatched inner dimensions throws");
}

// An empty row in the middle is a real case (a DOF with no contributions
// yet, e.g. before boundary terms are added) and a classic off-by-one in
// hand-written CSR construction.
void test_empty_rows() {
    SparseMatrixD a(4, 4);
    a.add(0, 0, 1.0);
    a.add(3, 3, 2.0);
    a.compress();

    check(a.nnz() == 2, "empty rows: nnz == 2");
    check(a.row_ptr()[1] == a.row_ptr()[2] && a.row_ptr()[2] == a.row_ptr()[3],
          "empty rows: interior empty rows have equal row_ptr bounds");
    const std::vector<double> y = a.matvec(std::vector<double>{1.0, 1.0, 1.0, 1.0});
    check(nearly(y[1], 0.0) && nearly(y[2], 0.0), "empty rows: matvec yields 0 there");
    check(nearly(y[0], 1.0) && nearly(y[3], 2.0), "empty rows: populated rows still correct");

    SparseMatrixD e(3, 3);
    e.compress();
    check(e.nnz() == 0 && e.is_zero(), "empty matrix: compresses to nnz 0 and is_zero");
}

}  // namespace

int main() {
    test_build_and_compress();
    test_matvec();
    test_transpose_and_multiply();
    test_triplet_export();
    test_complex_instantiation();
    test_phase_guards();
    test_empty_rows();

    std::cout << g_checks - g_failures << "/" << g_checks << " checks passed.\n";
    return g_failures == 0 ? 0 : 1;
}
