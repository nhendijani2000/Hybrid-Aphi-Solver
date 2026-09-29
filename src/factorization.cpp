#include "aphi_solver/factorization.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <stdexcept>
#include <string>

namespace aphi_solver {

using Complex = std::complex<double>;

SparseMatrixZ permuted_lower(const SparseSymmetricZ& a, const Permutation& p) {
    const int n = a.rows();
    if (p.size() != n) {
        throw std::invalid_argument("permuted_lower: the permutation is sized " +
                                    std::to_string(p.size()) + " for a matrix of " +
                                    std::to_string(n) + " rows");
    }

    // Every stored entry (r, c) of A's upper triangle becomes (iperm[r],
    // iperm[c]) and is normalised below the diagonal. A permutation moves
    // entries across it freely, so both orders occur.
    std::vector<int> row_ptr(static_cast<std::size_t>(n) + 1, 0);
    const auto placed = [&](int r, int c, int& row, int& col) {
        row = p.iperm[static_cast<std::size_t>(r)];
        col = p.iperm[static_cast<std::size_t>(c)];
        if (col > row) std::swap(row, col);
    };

    for (int r = 0; r < n; ++r) {
        for (int k = a.row_ptr()[static_cast<std::size_t>(r)];
             k < a.row_ptr()[static_cast<std::size_t>(r) + 1]; ++k) {
            int row = 0, col = 0;
            placed(r, a.col_index()[static_cast<std::size_t>(k)], row, col);
            ++row_ptr[static_cast<std::size_t>(row) + 1];
        }
    }
    for (int r = 0; r < n; ++r) {
        row_ptr[static_cast<std::size_t>(r) + 1] += row_ptr[static_cast<std::size_t>(r)];
    }

    std::vector<int> cursor(row_ptr.begin(), row_ptr.end() - 1);
    std::vector<int> col_index(static_cast<std::size_t>(row_ptr.back()), 0);
    std::vector<Complex> value(static_cast<std::size_t>(row_ptr.back()), Complex(0.0, 0.0));
    for (int r = 0; r < n; ++r) {
        for (int k = a.row_ptr()[static_cast<std::size_t>(r)];
             k < a.row_ptr()[static_cast<std::size_t>(r) + 1]; ++k) {
            int row = 0, col = 0;
            placed(r, a.col_index()[static_cast<std::size_t>(k)], row, col);
            const std::size_t at = static_cast<std::size_t>(cursor[static_cast<std::size_t>(row)]++);
            col_index[at] = col;
            value[at] = a.values()[static_cast<std::size_t>(k)];
        }
    }

    // Sort each row by column, carrying the values along.
    std::vector<int> order;
    SparseMatrixZ out(n, n);
    for (int r = 0; r < n; ++r) {
        const int begin = row_ptr[static_cast<std::size_t>(r)];
        const int end = row_ptr[static_cast<std::size_t>(r) + 1];
        order.resize(static_cast<std::size_t>(end - begin));
        for (int t = 0; t < end - begin; ++t) order[static_cast<std::size_t>(t)] = begin + t;
        std::sort(order.begin(), order.end(), [&](int x, int y) {
            return col_index[static_cast<std::size_t>(x)] < col_index[static_cast<std::size_t>(y)];
        });
        for (int t : order) {
            out.add(r, col_index[static_cast<std::size_t>(t)], value[static_cast<std::size_t>(t)]);
        }
    }
    out.compress();
    return out;
}

namespace {

/// The nonzero pattern of row `k` of `L`, in topological order, on `stack`
/// starting at the returned index.
///
/// Walk up the elimination tree from each `A(k, j<k)` until a node already
/// marked for this row is met, then reverse that path onto the output stack. The
/// reversal is what makes the result topological -- a descendant always precedes
/// its ancestor -- and the numeric loop below depends on that, because it uses
/// column `j` of `L` only after column `j` is complete.
int ereach(const SparseMatrixZ& lower, int k, const std::vector<int>& parent,
           std::vector<int>& flag, std::vector<int>& stack) {
    const int n = static_cast<int>(parent.size());
    int top = n;
    flag[static_cast<std::size_t>(k)] = k;
    for (int t = lower.row_ptr()[static_cast<std::size_t>(k)];
         t < lower.row_ptr()[static_cast<std::size_t>(k) + 1]; ++t) {
        int i = lower.col_index()[static_cast<std::size_t>(t)];
        if (i >= k) continue;  // the diagonal, handled separately
        int len = 0;
        while (flag[static_cast<std::size_t>(i)] != k) {
            stack[static_cast<std::size_t>(len++)] = i;
            flag[static_cast<std::size_t>(i)] = k;
            i = parent[static_cast<std::size_t>(i)];
            if (i < 0) break;
        }
        while (len > 0) stack[static_cast<std::size_t>(--top)] = stack[static_cast<std::size_t>(--len)];
    }
    return top;
}

}  // namespace

bool factorize_ldlt(const SparseSymmetricZ& a, const SolverAnalysis& analysis,
                    SymmetricFactor& out, FactorStats& stats, double pivot_floor) {
    const auto started = std::chrono::steady_clock::now();
    const int n = a.rows();
    if (analysis.permutation.size() != n || static_cast<int>(analysis.col_count.size()) != n) {
        throw std::invalid_argument(
            "factorize_ldlt: the analysis describes " +
            std::to_string(analysis.permutation.size()) + " unknowns, the matrix has " +
            std::to_string(n) + ". An analysis of a different matrix would lay out L's columns "
            "wrongly and the factorization would write into the wrong places.");
    }

    const SparseMatrixZ lower = permuted_lower(a, analysis.permutation);

    // L's storage, laid out from the column counts. The unit diagonal is not
    // stored, hence the -1 per column.
    out.rows = n;
    out.col_ptr.assign(static_cast<std::size_t>(n) + 1, 0);
    for (int j = 0; j < n; ++j) {
        out.col_ptr[static_cast<std::size_t>(j) + 1] =
            out.col_ptr[static_cast<std::size_t>(j)] +
            std::max(0, analysis.col_count[static_cast<std::size_t>(j)] - 1);
    }
    out.row_index.assign(static_cast<std::size_t>(out.col_ptr.back()), 0);
    out.value.assign(static_cast<std::size_t>(out.col_ptr.back()), Complex(0.0, 0.0));
    out.diagonal.assign(static_cast<std::size_t>(n), Complex(0.0, 0.0));

    std::vector<Complex> y(static_cast<std::size_t>(n), Complex(0.0, 0.0));
    std::vector<int> flag(static_cast<std::size_t>(n), -1);
    std::vector<int> stack(static_cast<std::size_t>(n), 0);
    std::vector<int> filled(static_cast<std::size_t>(n), 0);  // entries so far in each column

    stats = FactorStats{};
    stats.smallest_pivot = 0.0;
    bool first_pivot = true;

    for (int k = 0; k < n; ++k) {
        // Scatter row k of the permuted lower triangle, and take the diagonal
        // out of `y` so the loop below sees only the off-diagonal part.
        const int top = ereach(lower, k, analysis.parent, flag, stack);
        Complex d(0.0, 0.0);
        for (int t = lower.row_ptr()[static_cast<std::size_t>(k)];
             t < lower.row_ptr()[static_cast<std::size_t>(k) + 1]; ++t) {
            const int j = lower.col_index()[static_cast<std::size_t>(t)];
            const Complex v = lower.values()[static_cast<std::size_t>(t)];
            if (j == k) {
                d += v;
            } else if (j < k) {
                y[static_cast<std::size_t>(j)] += v;
            }
        }

        for (int s = top; s < n; ++s) {
            const int j = stack[static_cast<std::size_t>(s)];
            const Complex yj = y[static_cast<std::size_t>(j)];
            y[static_cast<std::size_t>(j)] = Complex(0.0, 0.0);

            // Apply column j of L, which is complete because the pattern is in
            // topological order and rows are processed in increasing k.
            const int begin = out.col_ptr[static_cast<std::size_t>(j)];
            const int end = begin + filled[static_cast<std::size_t>(j)];
            for (int t = begin; t < end; ++t) {
                y[static_cast<std::size_t>(out.row_index[static_cast<std::size_t>(t)])] -=
                    out.value[static_cast<std::size_t>(t)] * yj;
            }

            // No conjugation anywhere: A = A^T, not A^H.
            const Complex lkj = yj / out.diagonal[static_cast<std::size_t>(j)];
            d -= lkj * yj;

            if (end >= out.col_ptr[static_cast<std::size_t>(j) + 1]) {
                throw std::logic_error(
                    "factorize_ldlt: column " + std::to_string(j) +
                    " received more entries than the symbolic pass reserved for it. The "
                    "elimination tree and the column counts disagree.");
            }
            out.row_index[static_cast<std::size_t>(end)] = k;
            out.value[static_cast<std::size_t>(end)] = lkj;
            stats.largest_multiplier = std::max(stats.largest_multiplier, std::abs(lkj));
            ++filled[static_cast<std::size_t>(j)];
        }

        out.diagonal[static_cast<std::size_t>(k)] = d;
        const double magnitude = std::abs(d);
        stats.largest_pivot = std::max(stats.largest_pivot, magnitude);
        if (first_pivot || magnitude < stats.smallest_pivot) {
            stats.smallest_pivot = magnitude;
            first_pivot = false;
        }
        stats.columns_done = k + 1;

        // Without pivoting there is nothing to be done with a pivot this small,
        // and dividing by it would make every later entry meaningless. Stop and
        // say so. Step 7 perturbs instead and recovers by refinement.
        if (!(magnitude > pivot_floor)) {
            stats.ok = false;
            stats.nnz = out.nnz();
            stats.milliseconds =
                std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - started)
                    .count();
            return false;
        }
    }

    stats.ok = true;
    stats.nnz = out.nnz();
    stats.milliseconds =
        std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - started).count();
    return true;
}

SparseMatrixZ reconstruct(const SymmetricFactor& f) {
    const int n = f.rows;
    SparseMatrixZ out(n, n);
    // (L D L^T)(r, c) = sum_j L(r,j) D(j) L(c,j), with L's unit diagonal
    // implied. Build L by rows first so the sum is easy to form.
    std::vector<std::vector<std::pair<int, Complex>>> rows(static_cast<std::size_t>(n));
    for (int r = 0; r < n; ++r) rows[static_cast<std::size_t>(r)].push_back({r, Complex(1.0, 0.0)});
    for (int j = 0; j < n; ++j) {
        for (int t = f.col_ptr[static_cast<std::size_t>(j)];
             t < f.col_ptr[static_cast<std::size_t>(j) + 1]; ++t) {
            rows[static_cast<std::size_t>(f.row_index[static_cast<std::size_t>(t)])].push_back(
                {j, f.value[static_cast<std::size_t>(t)]});
        }
    }
    for (int r = 0; r < n; ++r) {
        for (int c = 0; c <= r; ++c) {
            Complex sum(0.0, 0.0);
            for (const auto& ra : rows[static_cast<std::size_t>(r)]) {
                for (const auto& ca : rows[static_cast<std::size_t>(c)]) {
                    if (ra.first != ca.first) continue;
                    sum += ra.second * f.diagonal[static_cast<std::size_t>(ra.first)] * ca.second;
                }
            }
            if (sum != Complex(0.0, 0.0)) {
                out.add(r, c, sum);
                if (r != c) out.add(c, r, sum);
            }
        }
    }
    out.compress();
    return out;
}


void solve_in_place(const SymmetricFactor& f, std::vector<Complex>& y) {
    const int n = f.rows;
    if (static_cast<int>(y.size()) != n) {
        throw std::invalid_argument("solve_in_place: the right-hand side has " +
                                    std::to_string(y.size()) + " entries for a factor of " +
                                    std::to_string(n) + " rows");
    }

    // L z = y. L is UNIT lower triangular, so there is nothing to divide by;
    // column-oriented because that is how L is stored.
    for (int j = 0; j < n; ++j) {
        const Complex yj = y[static_cast<std::size_t>(j)];
        for (int t = f.col_ptr[static_cast<std::size_t>(j)];
             t < f.col_ptr[static_cast<std::size_t>(j) + 1]; ++t) {
            y[static_cast<std::size_t>(f.row_index[static_cast<std::size_t>(t)])] -=
                f.value[static_cast<std::size_t>(t)] * yj;
        }
    }

    // D w = z.
    for (int j = 0; j < n; ++j) {
        const Complex d = f.diagonal[static_cast<std::size_t>(j)];
        if (d == Complex(0.0, 0.0)) {
            throw std::logic_error("solve_in_place: D[" + std::to_string(j) +
                                   "] is zero, so this factor cannot be used. factorize_ldlt "
                                   "returns false rather than producing one.");
        }
        y[static_cast<std::size_t>(j)] /= d;
    }

    // L^T x = w. The TRANSPOSE, not the conjugate transpose: A = A^T here.
    for (int j = n - 1; j >= 0; --j) {
        Complex acc = y[static_cast<std::size_t>(j)];
        for (int t = f.col_ptr[static_cast<std::size_t>(j)];
             t < f.col_ptr[static_cast<std::size_t>(j) + 1]; ++t) {
            acc -= f.value[static_cast<std::size_t>(t)] *
                   y[static_cast<std::size_t>(f.row_index[static_cast<std::size_t>(t)])];
        }
        y[static_cast<std::size_t>(j)] = acc;
    }
}

std::vector<Complex> solve(const SymmetricFactor& f, const Permutation& p,
                           const std::vector<Complex>& b) {
    const int n = f.rows;
    if (p.size() != n || static_cast<int>(b.size()) != n) {
        throw std::invalid_argument("solve: the factor, the permutation and the right-hand side "
                                    "must all describe the same number of unknowns");
    }
    std::vector<Complex> y(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        y[static_cast<std::size_t>(i)] = b[static_cast<std::size_t>(p.perm[static_cast<std::size_t>(i)])];
    }
    solve_in_place(f, y);
    std::vector<Complex> x(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        x[static_cast<std::size_t>(p.perm[static_cast<std::size_t>(i)])] = y[static_cast<std::size_t>(i)];
    }
    return x;
}

double relative_residual(const SparseSymmetricZ& a, const std::vector<Complex>& x,
                         const std::vector<Complex>& b) {
    if (static_cast<int>(x.size()) != a.rows() || x.size() != b.size()) {
        throw std::invalid_argument("relative_residual: lengths must match the matrix");
    }
    const std::vector<Complex> ax = a.matvec(x);
    double numerator = 0.0, denominator = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) {
        numerator += std::norm(ax[i] - b[i]);
        denominator += std::norm(b[i]);
    }
    numerator = std::sqrt(numerator);
    denominator = std::sqrt(denominator);
    return denominator > 0.0 ? numerator / denominator : numerator;
}


bool solve_symmetric(const SparseSymmetricZ& a, const std::vector<Complex>& b,
                     std::vector<Complex>& x, SolveReport& report, const SolveOptions& options,
                     const SolverAnalysis* precomputed) {
    const int n = a.rows();
    if (static_cast<int>(b.size()) != n) {
        throw std::invalid_argument("solve_symmetric: the right-hand side has " +
                                    std::to_string(b.size()) + " entries for a matrix of " +
                                    std::to_string(n) + " rows");
    }
    report = SolveReport{};

    // --- scale -----------------------------------------------------------
    // Symmetric, so A = A^T survives it -- which is what keeps LDL^T and the
    // upper-triangle storage applicable at all.
    std::vector<double> d(static_cast<std::size_t>(n), 1.0);
    if (options.equilibration_iterations > 0) {
        d = compute_symmetric_equilibration(a, options.equilibration_iterations);
    }
    report.equilibration_min = d.empty() ? 1.0 : *std::min_element(d.begin(), d.end());
    report.equilibration_max = d.empty() ? 1.0 : *std::max_element(d.begin(), d.end());

    const SparseSymmetricZ scaled =
        options.equilibration_iterations > 0 ? apply_symmetric_equilibration(a, d) : a;
    const std::vector<Complex> scaled_rhs =
        options.equilibration_iterations > 0 ? scale_rhs(b, d) : b;

    // --- order and factorize --------------------------------------------
    // The pattern is unchanged by scaling, and `analyze` reads only the
    // pattern -- so a caller that has already analysed this system can hand
    // the result in and skip it. `solve_mesh` does, because it analyses once
    // to report nnz(L) and would otherwise pay for the same AMD ordering
    // twice on every run: 9 s each on case 02.
    auto t0 = std::chrono::steady_clock::now();
    SolverAnalysis owned;
    if (precomputed == nullptr) owned = analyze(scaled, options.ordering);
    const SolverAnalysis& analysis = precomputed != nullptr ? *precomputed : owned;
    auto t1 = std::chrono::steady_clock::now();
    report.analyze_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

    SymmetricFactor factor;
    FactorStats stats;
    const bool factored = factorize_ldlt(scaled, analysis, factor, stats, options.pivot_floor);
    report.factorize_ms = stats.milliseconds;
    report.smallest_pivot = stats.smallest_pivot;
    report.largest_pivot = stats.largest_pivot;
    report.largest_multiplier = stats.largest_multiplier;
    report.factor_nnz = stats.nnz;
    report.columns_done = stats.columns_done;
    if (!factored) {
        report.ok = false;
        return false;
    }

    // --- solve, then UNDO the scaling ------------------------------------
    auto s0 = std::chrono::steady_clock::now();
    const std::vector<Complex> y = solve(factor, analysis.permutation, scaled_rhs);
    // x = diag(d) y. A multiplication, not a division, and not the identity.
    x = options.equilibration_iterations > 0 ? recover_equilibrated_solution(y, d) : y;
    auto s1 = std::chrono::steady_clock::now();
    report.solve_ms = std::chrono::duration<double, std::milli>(s1 - s0).count();

    // --- measure, against the ORIGINAL system ----------------------------
    report.residual = relative_residual(a, x, b);
    report.backward_error = backward_error(a, x, b);
    report.ok = true;
    return true;
}


double backward_error(const SparseSymmetricZ& a, const std::vector<Complex>& x,
                      const std::vector<Complex>& b) {
    if (static_cast<int>(x.size()) != a.rows() || x.size() != b.size()) {
        throw std::invalid_argument("backward_error: lengths must match the matrix");
    }
    const std::vector<Complex> ax = a.matvec(x);
    double numerator = 0.0, bnorm = 0.0, xnorm = 0.0, amax = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) {
        numerator += std::norm(ax[i] - b[i]);
        bnorm += std::norm(b[i]);
        xnorm += std::norm(x[i]);
    }
    for (const Complex& v : a.values()) amax = std::max(amax, std::abs(v));
    const double denominator = amax * std::sqrt(xnorm) + std::sqrt(bnorm);
    return denominator > 0.0 ? std::sqrt(numerator) / denominator : std::sqrt(numerator);
}

}  // namespace aphi_solver
