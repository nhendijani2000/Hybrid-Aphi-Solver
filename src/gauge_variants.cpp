#include "aphi_solver/gauge_variants.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace aphi_solver {

std::vector<double> dense_solve(std::vector<std::vector<double>> A, std::vector<double> b) {
    const int n = static_cast<int>(A.size());
    for (int col = 0; col < n; ++col) {
        int pivot_row = col;
        double pivot_val = std::abs(A[static_cast<std::size_t>(col)][static_cast<std::size_t>(col)]);
        for (int row = col + 1; row < n; ++row) {
            const double v = std::abs(A[static_cast<std::size_t>(row)][static_cast<std::size_t>(col)]);
            if (v > pivot_val) {
                pivot_row = row;
                pivot_val = v;
            }
        }
        if (pivot_val < 1e-12) {
            throw std::runtime_error("dense_solve: matrix is numerically singular");
        }
        if (pivot_row != col) {
            std::swap(A[static_cast<std::size_t>(pivot_row)], A[static_cast<std::size_t>(col)]);
            std::swap(b[static_cast<std::size_t>(pivot_row)], b[static_cast<std::size_t>(col)]);
        }
        const double pivot = A[static_cast<std::size_t>(col)][static_cast<std::size_t>(col)];
        for (int k = col; k < n; ++k) {
            A[static_cast<std::size_t>(col)][static_cast<std::size_t>(k)] /= pivot;
        }
        b[static_cast<std::size_t>(col)] /= pivot;
        for (int row = 0; row < n; ++row) {
            if (row == col) continue;
            const double factor = A[static_cast<std::size_t>(row)][static_cast<std::size_t>(col)];
            if (factor == 0.0) continue;
            for (int k = col; k < n; ++k) {
                A[static_cast<std::size_t>(row)][static_cast<std::size_t>(k)] -= factor * A[static_cast<std::size_t>(col)][static_cast<std::size_t>(k)];
            }
            b[static_cast<std::size_t>(row)] -= factor * b[static_cast<std::size_t>(col)];
        }
    }
    return b;
}

namespace {

/// Solves G_t x = b for the tree's own (implicit) node-incidence matrix, via
/// a single top-down walk of the spanning tree/forest in group-discovery
/// order: TreeCotreeResult::parent_group is guaranteed to only ever point to
/// a strictly smaller group id (a parent is always discovered before its
/// child), so processing groups 0..num_groups-1 in order always has a
/// group's parent value already computed. `b_by_edge` is indexed by GLOBAL
/// edge index; only its entries at tree edges are read. O(V) rather than an
/// O(V^3) dense inversion of G_t -- see docs/ENGINEERING_STANDARDS.md.
std::vector<double> solve_tree_system(const Mesh& mesh, const TreeCotreeResult& tc,
                                       const std::vector<int>& free_group_index, int num_free_groups,
                                       const std::vector<double>& b_by_edge) {
    std::vector<double> x(static_cast<std::size_t>(num_free_groups), 0.0);
    for (int g = 0; g < tc.num_groups; ++g) {
        if (tc.parent_group[static_cast<std::size_t>(g)] == -1) continue;  // reference group: x implicitly 0
        const int edge = tc.discovering_edge[static_cast<std::size_t>(g)];
        const auto& e = mesh.edges[static_cast<std::size_t>(edge)];
        const int v = (tc.node_group[static_cast<std::size_t>(e.first)] == g) ? e.first : e.second;
        const int p = tc.parent_group[static_cast<std::size_t>(g)];
        const int p_free = free_group_index[static_cast<std::size_t>(p)];
        const double xp = (p_free == -1) ? 0.0 : x[static_cast<std::size_t>(p_free)];
        const double bval = b_by_edge[static_cast<std::size_t>(edge)];
        // Canonical edge orientation is (e.first, e.second) with
        // e.first < e.second, giving G(edge, e.first) = -1, G(edge,
        // e.second) = +1 (build_gradient_matrix's own convention, reused
        // here group-wise): v == e.second means v sits at the "+1" end.
        const double xg = (v == e.second) ? (xp + bval) : (xp - bval);
        x[static_cast<std::size_t>(free_group_index[static_cast<std::size_t>(g)])] = xg;
    }
    return x;
}

}  // namespace

EssentialIncidenceMatrix compute_essential_incidence_matrix(const Mesh& mesh, const TreeCotreeResult& tc) {
    EssentialIncidenceMatrix result;
    const int num_groups = tc.num_groups;
    result.free_group_index.assign(static_cast<std::size_t>(num_groups), -1);
    int num_free = 0;
    for (int g = 0; g < num_groups; ++g) {
        if (tc.parent_group[static_cast<std::size_t>(g)] != -1) {
            result.free_group_index[static_cast<std::size_t>(g)] = num_free++;
        }
    }
    result.num_free_groups = num_free;

    const int num_edges = mesh.num_edges();
    result.cotree_local_index.assign(static_cast<std::size_t>(num_edges), -1);
    int num_cotree = 0;
    for (int e = 0; e < num_edges; ++e) {
        if (!tc.is_tree_edge[static_cast<std::size_t>(e)]) {
            result.cotree_local_index[static_cast<std::size_t>(e)] = num_cotree++;
        }
    }
    result.num_cotree_edges = num_cotree;
    result.values.assign(static_cast<std::size_t>(num_cotree),
                          std::vector<double>(static_cast<std::size_t>(num_free), 0.0));

    // One O(V) tree walk per free group gives one COLUMN of G_t^{-1}
    // (g_t_inv_columns[k][m] == G_t^{-1}[m, k]).
    std::vector<double> unit_rhs(static_cast<std::size_t>(num_edges), 0.0);
    std::vector<std::vector<double>> g_t_inv_columns(static_cast<std::size_t>(num_free));
    for (int g = 0; g < num_groups; ++g) {
        if (tc.parent_group[static_cast<std::size_t>(g)] == -1) continue;
        const int col = result.free_group_index[static_cast<std::size_t>(g)];
        std::fill(unit_rhs.begin(), unit_rhs.end(), 0.0);
        unit_rhs[static_cast<std::size_t>(tc.discovering_edge[static_cast<std::size_t>(g)])] = 1.0;
        g_t_inv_columns[static_cast<std::size_t>(col)] =
            solve_tree_system(mesh, tc, result.free_group_index, num_free, unit_rhs);
    }

    // F = G_c * G_t^{-1}: F[row, col] = -G_t^{-1}[fi, col] + G_t^{-1}[fj, col]
    // = -g_t_inv_columns[col][fi] + g_t_inv_columns[col][fj], for cotree
    // edge `row` = (i, j) with free-group indices fi, fj (a reference-group
    // endpoint contributes 0, since its potential is fixed, not a free
    // unknown).
    for (int e = 0; e < num_edges; ++e) {
        if (tc.is_tree_edge[static_cast<std::size_t>(e)]) continue;
        const int row = result.cotree_local_index[static_cast<std::size_t>(e)];
        const auto& edge = mesh.edges[static_cast<std::size_t>(e)];
        const int gi = tc.node_group[static_cast<std::size_t>(edge.first)];
        const int gj = tc.node_group[static_cast<std::size_t>(edge.second)];
        const int fi = result.free_group_index[static_cast<std::size_t>(gi)];
        const int fj = result.free_group_index[static_cast<std::size_t>(gj)];
        for (int col = 0; col < num_free; ++col) {
            double val = 0.0;
            if (fi != -1) val -= g_t_inv_columns[static_cast<std::size_t>(col)][static_cast<std::size_t>(fi)];
            if (fj != -1) val += g_t_inv_columns[static_cast<std::size_t>(col)][static_cast<std::size_t>(fj)];
            result.values[static_cast<std::size_t>(row)][static_cast<std::size_t>(col)] = val;
        }
    }
    return result;
}

std::vector<double> select_cotree_entries(const std::vector<double>& full, const std::vector<int>& cotree_local_index) {
    int num_cotree = 0;
    for (int idx : cotree_local_index) {
        if (idx != -1) ++num_cotree;
    }
    std::vector<double> result(static_cast<std::size_t>(num_cotree));
    for (std::size_t e = 0; e < cotree_local_index.size(); ++e) {
        if (cotree_local_index[e] != -1) {
            result[static_cast<std::size_t>(cotree_local_index[e])] = full[e];
        }
    }
    return result;
}

GaugeVariant build_albanese_rubinacci_gauge(const SparseMatrix& M, const TreeCotreeResult& tc) {
    GaugeVariant variant;
    variant.name = "Albanese-Rubinacci";
    const int num_edges = M.rows();
    variant.cotree_local_index.assign(static_cast<std::size_t>(num_edges), -1);
    int num_cotree = 0;
    for (int e = 0; e < num_edges; ++e) {
        if (!tc.is_tree_edge[static_cast<std::size_t>(e)]) {
            variant.cotree_local_index[static_cast<std::size_t>(e)] = num_cotree++;
        }
    }

    // Walk M's CSR rows directly, skipping whole tree rows before touching
    // their entries at all -- a tree row contributes nothing under a_t = 0.
    variant.reduced_matrix = SparseMatrix(num_cotree, num_cotree);
    const auto& row_ptr = M.row_ptr();
    const auto& col_index = M.col_index();
    const auto& values = M.values();
    for (int e = 0; e < num_edges; ++e) {
        const int r = variant.cotree_local_index[static_cast<std::size_t>(e)];
        if (r == -1) continue;  // tree row: dropped, per a_t = 0
        for (int k = row_ptr[static_cast<std::size_t>(e)]; k < row_ptr[static_cast<std::size_t>(e) + 1]; ++k) {
            const int c = variant.cotree_local_index[static_cast<std::size_t>(col_index[static_cast<std::size_t>(k)])];
            if (c == -1) continue;  // tree column: dropped for the same reason
            variant.reduced_matrix.add(r, c, values[static_cast<std::size_t>(k)]);
        }
    }
    variant.reduced_matrix.compress();
    return variant;
}

std::vector<double> recover_albanese_rubinacci_solution(const std::vector<double>& a_c, const GaugeVariant& variant,
                                                          int num_edges) {
    std::vector<double> a(static_cast<std::size_t>(num_edges), 0.0);
    for (int e = 0; e < num_edges; ++e) {
        const int c = variant.cotree_local_index[static_cast<std::size_t>(e)];
        if (c != -1) a[static_cast<std::size_t>(e)] = a_c[static_cast<std::size_t>(c)];
    }
    return a;
}

GaugeVariant build_munteanu_unsymmetric_gauge(const SparseMatrix& M, const TreeCotreeResult& tc,
                                               const EssentialIncidenceMatrix& F) {
    GaugeVariant variant;
    variant.name = "Munteanu unsymmetric";
    const int num_edges = M.rows();
    const int num_cotree = F.num_cotree_edges;
    variant.cotree_local_index = F.cotree_local_index;

    // Dense L^T (E x num_cotree): identity on cotree rows, -F^T on tree rows
    // (a = L^T a_c reconstructs the full solution -- see
    // recover_munteanu_unsymmetric_solution, which does the same thing for
    // an actual solved a_c rather than building the whole operator).
    std::vector<std::vector<double>> Lt(static_cast<std::size_t>(num_edges),
                                         std::vector<double>(static_cast<std::size_t>(num_cotree), 0.0));
    for (int e = 0; e < num_edges; ++e) {
        if (!tc.is_tree_edge[static_cast<std::size_t>(e)]) {
            Lt[static_cast<std::size_t>(e)][static_cast<std::size_t>(F.cotree_local_index[static_cast<std::size_t>(e)])] = 1.0;
        }
    }
    std::vector<int> tree_edge_to_group(static_cast<std::size_t>(num_edges), -1);
    for (int g = 0; g < tc.num_groups; ++g) {
        if (tc.parent_group[static_cast<std::size_t>(g)] != -1) {
            tree_edge_to_group[static_cast<std::size_t>(tc.discovering_edge[static_cast<std::size_t>(g)])] = g;
        }
    }
    for (int e = 0; e < num_edges; ++e) {
        if (!tc.is_tree_edge[static_cast<std::size_t>(e)]) continue;
        const int g = tree_edge_to_group[static_cast<std::size_t>(e)];
        const int free_col = F.free_group_index[static_cast<std::size_t>(g)];
        for (int col = 0; col < num_cotree; ++col) {
            Lt[static_cast<std::size_t>(e)][static_cast<std::size_t>(col)] =
                -F.values[static_cast<std::size_t>(col)][static_cast<std::size_t>(free_col)];
        }
    }

    // Reduced matrix = (select cotree rows of M) * L^T -- an oblique
    // Petrov-Galerkin projection (test space = plain cotree selector, trial
    // space = L^T), computed here via dense multiplication at this phase's
    // test-mesh scale (see docs/ENGINEERING_STANDARDS.md).
    const auto M_dense = M.to_dense();
    std::vector<std::vector<double>> product(static_cast<std::size_t>(num_edges),
                                              std::vector<double>(static_cast<std::size_t>(num_cotree), 0.0));
    for (int r = 0; r < num_edges; ++r) {
        for (int k = 0; k < num_edges; ++k) {
            const double mrk = M_dense[static_cast<std::size_t>(r)][static_cast<std::size_t>(k)];
            if (mrk == 0.0) continue;
            for (int c = 0; c < num_cotree; ++c) {
                product[static_cast<std::size_t>(r)][static_cast<std::size_t>(c)] += mrk * Lt[static_cast<std::size_t>(k)][static_cast<std::size_t>(c)];
            }
        }
    }

    variant.reduced_matrix = SparseMatrix(num_cotree, num_cotree);
    for (int e = 0; e < num_edges; ++e) {
        if (tc.is_tree_edge[static_cast<std::size_t>(e)]) continue;
        const int row = F.cotree_local_index[static_cast<std::size_t>(e)];
        for (int c = 0; c < num_cotree; ++c) {
            const double val = product[static_cast<std::size_t>(e)][static_cast<std::size_t>(c)];
            if (val != 0.0) variant.reduced_matrix.add(row, c, val);
        }
    }
    variant.reduced_matrix.compress();
    return variant;
}

std::vector<double> recover_munteanu_unsymmetric_solution(const std::vector<double>& a_c, const GaugeVariant& variant,
                                                            const EssentialIncidenceMatrix& F,
                                                            const TreeCotreeResult& tc, int num_edges) {
    (void)variant;
    std::vector<double> a(static_cast<std::size_t>(num_edges), 0.0);
    for (int e = 0; e < num_edges; ++e) {
        const int c = F.cotree_local_index[static_cast<std::size_t>(e)];
        if (c != -1) a[static_cast<std::size_t>(e)] = a_c[static_cast<std::size_t>(c)];
    }
    for (int g = 0; g < tc.num_groups; ++g) {
        if (tc.parent_group[static_cast<std::size_t>(g)] == -1) continue;
        const int free_col = F.free_group_index[static_cast<std::size_t>(g)];
        double at = 0.0;
        for (int col = 0; col < F.num_cotree_edges; ++col) {
            at -= F.values[static_cast<std::size_t>(col)][static_cast<std::size_t>(free_col)] * a_c[static_cast<std::size_t>(col)];
        }
        a[static_cast<std::size_t>(tc.discovering_edge[static_cast<std::size_t>(g)])] = at;
    }
    return a;
}

namespace {

// B*v for B = A^T*A, via two sparse matvecs -- B itself is never formed, so
// this stays O(nnz) regardless of how dense A^T*A would be if it were ever
// written out (which, for a tree-cotree-reduced matrix, can be
// substantially denser than A itself -- see docs/TREE_COTREE_GAUGE.md
// Sec. 4-5 on Method D's fill-in). The two matvec primitives used to be
// hand-written here over the old COO triplets; they are now `Sparse<T>`'s
// own CSR members (docs/ROADMAP.md Phase 03.5, step 2), so this is the only
// piece left that is specific to the condition-number estimate.
std::vector<double> ata_matvec(const SparseMatrix& A, const std::vector<double>& v) {
    return A.matvec_transpose(A.matvec(v));
}

double dot(const std::vector<double>& a, const std::vector<double>& b) {
    double s = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) s += a[i] * b[i];
    return s;
}

double norm2(const std::vector<double>& v) { return std::sqrt(dot(v, v)); }

void normalize(std::vector<double>& v) {
    const double n = norm2(v);
    if (n > 1e-300) {
        for (auto& x : v) x /= n;
    }
}

// Solves B*x = b via unpreconditioned Conjugate Gradient (Hestenes &
// Stiefel, 1952 -- see docs/REFERENCES.md, "Matrix conditioning"), using
// only ata_matvec. B = A^T*A is symmetric positive semi-definite by
// construction, which is all CG needs; B is never formed or factored, so
// this replaces what used to be a fresh O(n^3) dense Gaussian elimination
// on every single inverse-power-iteration step (see
// docs/ENGINEERING_STANDARDS.md for why this was brought forward from its
// originally-planned Phase 04/05 slot) with an O(nnz)-per-iteration sparse
// solve. Returns false, with x left at its last iterate, if CG breaks down
// on a near-zero curvature direction p^T*B*p -- the standard symptom of B
// being (numerically) singular along that direction; the caller reports
// that as +infinity rather than a bogus finite condition number. Running
// out of cg_max_iterations without reaching cg_tol is NOT treated as
// failure: inverse power iteration only needs a decent improving direction
// each outer step ("inexact inverse iteration" is a standard technique,
// e.g. Golub & Van Loan, "Matrix Computations"), so an inexact inner solve
// is fine and the outer Rayleigh-quotient loop's own convergence check is
// what ultimately decides when lambda_min has settled.
bool cg_solve_ata(const SparseMatrix& A, const std::vector<double>& b, std::vector<double>& x,
                   int max_cg_iterations, double cg_tol) {
    const std::size_t n = b.size();
    x.assign(n, 0.0);
    std::vector<double> r = b;
    double rs_old = dot(r, r);
    const double b_norm = std::max(1.0, norm2(b));
    if (std::sqrt(rs_old) < cg_tol * b_norm) return true;
    std::vector<double> p = r;
    for (int iter = 0; iter < max_cg_iterations; ++iter) {
        const std::vector<double> Bp = ata_matvec(A, p);
        const double pBp = dot(p, Bp);
        if (std::abs(pBp) < 1e-300) return false;
        const double alpha = rs_old / pBp;
        for (std::size_t i = 0; i < n; ++i) {
            x[i] += alpha * p[i];
            r[i] -= alpha * Bp[i];
        }
        const double rs_new = dot(r, r);
        if (std::sqrt(rs_new) < cg_tol * b_norm) return true;
        const double beta = rs_new / rs_old;
        for (std::size_t i = 0; i < n; ++i) p[i] = r[i] + beta * p[i];
        rs_old = rs_new;
    }
    return true;
}

}  // namespace

double estimate_condition_number(const SparseMatrix& A, int max_iterations, double tol) {
    const int n = A.rows();
    if (n <= 1) return 1.0;

    // Largest eigenvalue of B = A^T*A (== largest squared singular value of
    // A) via plain power iteration. B is never formed -- ata_matvec is two
    // O(nnz) sparse passes over A, not an O(n^2) dense product.
    std::vector<double> v(static_cast<std::size_t>(n), 1.0 / std::sqrt(static_cast<double>(n)));
    double lambda_max = 0.0;
    for (int iter = 0; iter < max_iterations; ++iter) {
        std::vector<double> w = ata_matvec(A, v);
        normalize(w);
        const double next = dot(w, ata_matvec(A, w));
        v = w;
        if (std::abs(next - lambda_max) < tol * std::max(1.0, std::abs(next))) {
            lambda_max = next;
            break;
        }
        lambda_max = next;
    }

    // Smallest eigenvalue of B via inverse power iteration -- same
    // algorithm as before, but each "solve B*x = u" step now uses
    // cg_solve_ata (sparse, O(nnz) per CG iteration) instead of dense_solve
    // (dense O(n^3) Gaussian elimination from scratch every single outer
    // step).
    std::vector<double> u(static_cast<std::size_t>(n), 1.0 / std::sqrt(static_cast<double>(n)));
    double lambda_min = 0.0;
    bool singular = false;
    const int cg_max_iterations = std::max(50, std::min(n, 2000));
    const double cg_tol = 1e-10;
    for (int iter = 0; iter < max_iterations; ++iter) {
        std::vector<double> w;
        const bool ok = cg_solve_ata(A, u, w, cg_max_iterations, cg_tol);
        if (!ok) {
            singular = true;
            break;
        }
        normalize(w);
        const double next = dot(w, ata_matvec(A, w));
        u = w;
        if (std::abs(next - lambda_min) < tol * std::max(1.0, std::abs(next))) {
            lambda_min = next;
            break;
        }
        lambda_min = next;
    }

    if (singular || lambda_min <= 0.0) {
        return std::numeric_limits<double>::infinity();
    }
    return std::sqrt(lambda_max / lambda_min);
}

}  // namespace aphi_solver
