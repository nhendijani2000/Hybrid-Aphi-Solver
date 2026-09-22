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

    // The unknowns are the cotree edges that are not on a Dirichlet
    // (n x A = 0) surface: tree edges are the gauge, Dirichlet edges are
    // fixed at zero by the boundary condition. Every Dirichlet edge joins two
    // nodes of the same surface -- the same group -- so it has no row in the
    // group-level incidence matrix, and dropping it changes nothing in F
    // beyond removing its (empty) row.
    const int num_edges = mesh.num_edges();
    result.cotree_local_index.assign(static_cast<std::size_t>(num_edges), -1);
    int num_cotree = 0;
    for (int e = 0; e < num_edges; ++e) {
        if (!tc.is_tree_edge[static_cast<std::size_t>(e)] && !tc.is_dirichlet_edge[static_cast<std::size_t>(e)]) {
            result.cotree_local_index[static_cast<std::size_t>(e)] = num_cotree++;
        }
    }
    result.num_cotree_edges = num_cotree;

    // Per-group tree bookkeeping, one pass. parent_group[g] is always a
    // strictly smaller id (a parent is discovered before its child), so
    // depth can be filled in increasing g without a traversal.
    //
    // sign[g] is the +-1 picked up when descending through group g's
    // discovering edge: +1 when g's own node is that edge's canonical second
    // (higher-index) endpoint, -1 when it is the first -- since
    // build_gradient_matrix puts G(edge, first) = -1 and G(edge, second) =
    // +1, read group-wise.
    std::vector<int> depth(static_cast<std::size_t>(num_groups), 0);
    std::vector<int> sign(static_cast<std::size_t>(num_groups), 0);
    for (int g = 0; g < num_groups; ++g) {
        const int p = tc.parent_group[static_cast<std::size_t>(g)];
        if (p == -1) continue;  // root group: depth 0, never contributes a column
        depth[static_cast<std::size_t>(g)] = depth[static_cast<std::size_t>(p)] + 1;
        const auto& e = mesh.edges[static_cast<std::size_t>(tc.discovering_edge[static_cast<std::size_t>(g)])];
        const int v = (tc.node_group[static_cast<std::size_t>(e.first)] == g) ? e.first : e.second;
        sign[static_cast<std::size_t>(g)] = (v == e.second) ? 1 : -1;
    }

    // Each cotree edge's row of F is its fundamental cycle: walk both
    // endpoints' groups up the tree to their common ancestor, emitting one
    // entry per tree edge stepped over. Everything above the common ancestor
    // is shared by both endpoints and cancels, so it is never visited --
    // which is why this is O(cycle length) per row rather than O(V) per free
    // group. See EssentialIncidenceMatrix::F for the derivation.
    result.F = SparseMatrix(num_cotree, num_free);
    for (int e = 0; e < num_edges; ++e) {
        const int row = result.cotree_local_index[static_cast<std::size_t>(e)];
        if (row < 0) continue;  // tree or Dirichlet edge: not an unknown
        const auto& edge = mesh.edges[static_cast<std::size_t>(e)];
        int gi = tc.node_group[static_cast<std::size_t>(edge.first)];
        int gj = tc.node_group[static_cast<std::size_t>(edge.second)];

        auto emit = [&](int g, double s) {
            result.F.add(row, result.free_group_index[static_cast<std::size_t>(g)],
                          s * static_cast<double>(sign[static_cast<std::size_t>(g)]));
        };

        while (gi != gj) {
            const int di = depth[static_cast<std::size_t>(gi)];
            const int dj = depth[static_cast<std::size_t>(gj)];
            if (di > dj) {
                emit(gi, -1.0);
                gi = tc.parent_group[static_cast<std::size_t>(gi)];
            } else if (dj > di) {
                emit(gj, +1.0);
                gj = tc.parent_group[static_cast<std::size_t>(gj)];
            } else {
                // Equal depth and still distinct. If both are roots, the two
                // endpoints lie in different connected pieces of the mesh --
                // each piece has its own root -- so there is no common
                // ancestor, nothing cancels, and both paths are already fully
                // emitted. (An edge cannot actually join two pieces, so this
                // is a guard rather than a live case.)
                if (tc.parent_group[static_cast<std::size_t>(gi)] == -1 ||
                    tc.parent_group[static_cast<std::size_t>(gj)] == -1) {
                    break;
                }
                emit(gi, -1.0);
                gi = tc.parent_group[static_cast<std::size_t>(gi)];
                emit(gj, +1.0);
                gj = tc.parent_group[static_cast<std::size_t>(gj)];
            }
        }
    }
    result.F.compress();
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

GaugeIndexMap build_albanese_rubinacci_index_map(const std::vector<int>& a_dof_edge, int num_phi_dofs,
                                                   const TreeCotreeResult& tc) {
    const int num_a = static_cast<int>(a_dof_edge.size());
    GaugeIndexMap map;
    map.full_to_reduced.assign(static_cast<std::size_t>(num_a + num_phi_dofs), -1);
    map.reduced_to_full.reserve(static_cast<std::size_t>(num_a + num_phi_dofs));

    for (int k = 0; k < num_a; ++k) {
        const int edge = a_dof_edge[static_cast<std::size_t>(k)];
        // Eliminated: a tree edge (the gauge, a_t = 0) or an edge on an
        // n x A = 0 surface (the boundary condition). Both are fixed at zero.
        if (tc.is_tree_edge[static_cast<std::size_t>(edge)] || tc.is_dirichlet_edge[static_cast<std::size_t>(edge)]) {
            continue;
        }
        map.full_to_reduced[static_cast<std::size_t>(k)] = static_cast<int>(map.reduced_to_full.size());
        map.reduced_to_full.push_back(k);
    }
    // Phi DOFs are never gauged -- the gauge freedom removed here is A's
    // alone (docs/FORMULATION.md Sec 5.3) -- so every one of them survives.
    for (int p = 0; p < num_phi_dofs; ++p) {
        const int full = num_a + p;
        map.full_to_reduced[static_cast<std::size_t>(full)] = static_cast<int>(map.reduced_to_full.size());
        map.reduced_to_full.push_back(full);
    }

    map.reduced_size = static_cast<int>(map.reduced_to_full.size());
    return map;
}

GaugeVariant build_albanese_rubinacci_gauge(const SparseMatrix& M, const TreeCotreeResult& tc) {
    GaugeVariant variant;
    variant.name = "Albanese-Rubinacci";
    // The bare edge-indexed case: every matrix index IS a mesh edge, so the
    // A-DOF -> edge map is the identity and there are no Phi DOFs. Expressed
    // through the same index map and the same principal-submatrix primitive
    // that a coupled [a; Phi] system uses, so the reduction itself has one
    // implementation rather than two.
    const int num_edges = M.rows();
    std::vector<int> identity(static_cast<std::size_t>(num_edges));
    for (int e = 0; e < num_edges; ++e) identity[static_cast<std::size_t>(e)] = e;

    const GaugeIndexMap map = build_albanese_rubinacci_index_map(identity, 0, tc);
    variant.cotree_local_index = map.full_to_reduced;
    variant.reduced_matrix = M.principal_submatrix(map.full_to_reduced, map.reduced_size);
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

    // L^T (E x num_cotree), sparse: the identity on cotree rows, -F^T on
    // tree rows, so a = L^T a_c reconstructs the full edge solution. The
    // tree row belonging to group g is minus F's column for that group,
    // which is one row of F^T -- so the whole tree part is just F
    // transposed, scattered to the discovering edges.
    // Dirichlet edges -- surface tree edges included -- get no row at all:
    // they are zero by the boundary condition, so they contribute nothing to
    // the reconstructed field.
    SparseMatrix Lt(num_edges, num_cotree);
    for (int e = 0; e < num_edges; ++e) {
        const int c = F.cotree_local_index[static_cast<std::size_t>(e)];
        if (c >= 0) Lt.add(e, c, 1.0);
    }
    const SparseMatrix Ft = F.F.transposed();  // num_free_groups x num_cotree
    const auto& ft_row_ptr = Ft.row_ptr();
    const auto& ft_col = Ft.col_index();
    const auto& ft_val = Ft.values();
    for (int g = 0; g < tc.num_groups; ++g) {
        if (tc.parent_group[static_cast<std::size_t>(g)] == -1) continue;
        const int e = tc.discovering_edge[static_cast<std::size_t>(g)];
        const int free_row = F.free_group_index[static_cast<std::size_t>(g)];
        for (int k = ft_row_ptr[static_cast<std::size_t>(free_row)];
             k < ft_row_ptr[static_cast<std::size_t>(free_row) + 1]; ++k) {
            Lt.add(e, ft_col[static_cast<std::size_t>(k)], -ft_val[static_cast<std::size_t>(k)]);
        }
    }
    Lt.compress();

    // Reduced matrix = (select cotree rows of M) * L^T -- an oblique
    // Petrov-Galerkin projection (test space = plain cotree selector, trial
    // space = L^T). Both factors are sparse and the product goes through
    // Sparse::multiply (Gustavson), so no dense E x E intermediate is ever
    // formed. It used to densify M, which is ~80 GB at 10^5 edges and the
    // reason Method D could not run on a real mesh at all.
    SparseMatrix ScM(num_cotree, num_edges);
    const auto& m_row_ptr = M.row_ptr();
    const auto& m_col = M.col_index();
    const auto& m_val = M.values();
    for (int e = 0; e < num_edges; ++e) {
        const int row = F.cotree_local_index[static_cast<std::size_t>(e)];
        if (row == -1) continue;  // tree row: not in the test space
        for (int k = m_row_ptr[static_cast<std::size_t>(e)]; k < m_row_ptr[static_cast<std::size_t>(e) + 1]; ++k) {
            ScM.add(row, m_col[static_cast<std::size_t>(k)], m_val[static_cast<std::size_t>(k)]);
        }
    }
    ScM.compress();

    variant.reduced_matrix = ScM.multiply(Lt);
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
    // a_t = -F^T a_c, for every tree edge at once: one O(nnz) transposed
    // matvec instead of a per-group loop over every cotree column, which
    // was scanning the dense F column by column.
    const std::vector<double> ft_ac = F.F.matvec_transpose(a_c);
    for (int g = 0; g < tc.num_groups; ++g) {
        if (tc.parent_group[static_cast<std::size_t>(g)] == -1) continue;
        const int free_col = F.free_group_index[static_cast<std::size_t>(g)];
        a[static_cast<std::size_t>(tc.discovering_edge[static_cast<std::size_t>(g)])] =
            -ft_ac[static_cast<std::size_t>(free_col)];
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
