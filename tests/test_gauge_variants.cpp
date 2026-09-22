// Minimal, dependency-free test runner -- see the comment in
// tests/CMakeLists.txt for why no external framework is fetched here.
//
// Covers docs/ROADMAP.md Phase 03 steps 2-3: the Albanese-Rubinacci and
// Munteanu-unsymmetric tree-cotree gauge variants (I. Munteanu, "Tree-cotree
// condensation properties"), built on top of Phase 03 step 1's spanning-tree
// decomposition, plus the crude power-iteration condition-number estimate.

#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <iostream>
#include <string>
#include <vector>

#include "aphi_solver/complex_matrix.hpp"
#include "aphi_solver/gauge_variants.hpp"
#include "aphi_solver/incidence.hpp"
#include "aphi_solver/mesh.hpp"
#include "aphi_solver/tree_cotree.hpp"

using namespace aphi_solver;

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

bool nearly(double a, double b, double tol = 1e-6) { return std::abs(a - b) < tol; }
bool nearly(Complex a, Complex b, double tol = 1e-6) { return std::abs(a - b) < tol; }

Mesh make_single_tet() {
    Mesh m;
    m.nodes = {Vec3(0, 0, 0), Vec3(1, 0, 0), Vec3(0, 1, 0), Vec3(0, 0, 1)};
    m.tets = {{0, 1, 2, 3}};
    m.build_topology();
    return m;
}

Mesh make_two_tets_sharing_a_face() {
    Mesh m;
    m.nodes = {Vec3(0, 0, 0), Vec3(1, 0, 0), Vec3(0, 1, 0), Vec3(0, 0, 1), Vec3(1, 1, 1)};
    m.tets = {{0, 1, 2, 3}, {1, 2, 3, 4}};
    m.build_topology();
    return m;
}

// Transpose and multiply are members of the sparse type itself since the
// Phase 03.5 migration to Sparse<T> (docs/ROADMAP.md) -- no local copy and
// no free-function wrapper here anymore.

std::vector<double> dense_matvec(const std::vector<std::vector<double>>& A, const std::vector<double>& x) {
    std::vector<double> y(A.size(), 0.0);
    for (std::size_t i = 0; i < A.size(); ++i) {
        double s = 0.0;
        for (std::size_t j = 0; j < A[i].size(); ++j) s += A[i][j] * x[j];
        y[i] = s;
    }
    return y;
}

// Curl-curl "vacuum" test matrix M = C^T * C (unit reluctivity nu = 1) --
// no assembly pipeline exists yet (Phase 04), so this is the correct,
// physically legitimate stand-in curl-curl matrix for a vacuum test
// problem, not a placeholder invented for testing.
SparseMatrix build_test_M(const Mesh& mesh) {
    const SparseMatrix C = build_curl_matrix(mesh);
    return C.transposed().multiply(C);
}

// Independent verification of estimate_condition_number itself (Sept 2026,
// after it was rewritten from a dense O(n^3)-per-iteration implementation
// to a sparse, matrix-free Conjugate-Gradient-based one -- see
// docs/ENGINEERING_STANDARDS.md and docs/REFERENCES.md, "Matrix
// conditioning"): checked here against DIAGONAL matrices, whose singular
// values -- and hence whose exact condition number -- are known
// analytically (the diagonal entries themselves), independent of anything
// this project's gauge machinery computes. This is deliberately separate
// from the tree-cotree gauge checks above: those confirm the gauge
// reduction is correct; this confirms the *estimator* is correct, on
// inputs where the right answer is known by construction rather than by
// trusting the estimator's own output.
SparseMatrix make_diagonal_matrix(const std::vector<double>& diag_values) {
    const int n = static_cast<int>(diag_values.size());
    SparseMatrix A(n, n);
    for (std::size_t i = 0; i < diag_values.size(); ++i) {
        A.add(static_cast<int>(i), static_cast<int>(i), diag_values[i]);
    }
    A.compress();
    return A;
}

void run_condition_number_analytic_checks() {
    // Identity: every singular value is 1, so kappa = 1 exactly.
    {
        const SparseMatrix I = make_diagonal_matrix({1.0, 1.0, 1.0, 1.0, 1.0});
        const double kappa = estimate_condition_number(I);
        check(nearly(kappa, 1.0, 1e-6), "condition-number estimator: 5x5 identity gives kappa = 1");
    }

    // Small diagonal with a hand-picked, exactly-known condition number.
    {
        const SparseMatrix D = make_diagonal_matrix({1.0, 2.0, 4.0, 100.0});
        const double kappa = estimate_condition_number(D);
        check(nearly(kappa, 100.0, 1e-4), "condition-number estimator: diag(1,2,4,100) gives kappa = 100");
    }

    // Larger diagonal (n = 500, singular values spread geometrically from 1
    // to 1000, kappa = 1000 exactly) -- exercises the estimator at a scale
    // where the OLD dense implementation (O(n^2) memory just for A^T*A,
    // O(n^3) per inverse-power-iteration step) would already be a
    // noticeably heavy single test to run; the sparse rewrite handles it as
    // one more matvec-bound case. Timed and printed (not hard-asserted --
    // machine speed varies) as a visible sanity check that this stays fast.
    {
        const int n = 500;
        std::vector<double> diag_values(static_cast<std::size_t>(n));
        for (int i = 0; i < n; ++i) {
            // Geometric spread from 1 to 1000 across n values.
            diag_values[static_cast<std::size_t>(i)] = std::pow(1000.0, static_cast<double>(i) / static_cast<double>(n - 1));
        }
        const SparseMatrix D = make_diagonal_matrix(diag_values);
        const auto t0 = std::chrono::steady_clock::now();
        const double kappa = estimate_condition_number(D);
        const auto t1 = std::chrono::steady_clock::now();
        const double elapsed_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        check(nearly(kappa, 1000.0, 1.0), "condition-number estimator: 500x500 diagonal (kappa=1000 exactly) recovered");
        std::cout << "500x500 diagonal test: kappa estimate = " << kappa << " (exact: 1000), " << elapsed_ms << " ms\n";
    }
}

// docs/ROADMAP.md Phase 03.5, step 4: the gauge reduction must act on the
// A-DOF rows/columns of a larger coupled [a; Phi] system while leaving Phi
// alone, instead of assuming the matrix index space IS the mesh edge index
// space. Two things make this more than a reindexing exercise, and both are
// exercised below:
//
//  - A-DOFs may be a STRICT SUBSET of the mesh's edges (PEC tangential edges
//    are removed from the unknown set, docs/FORMULATION.md Sec 5.4), so the
//    A-DOF -> edge map is not the identity.
//  - The system matrix is COMPLEX, while the tree/cotree decomposition it is
//    reduced by is pure topology. The index map carries no scalar type at
//    all, which is what lets one reduction serve both.
void test_coupled_system_reduction() {
    const Mesh mesh = make_two_tets_sharing_a_face();
    const TreeCotreeResult tc = build_tree_cotree(mesh, std::vector<bool>(static_cast<std::size_t>(mesh.num_edges()), false));

    // Drop edge 0 from the unknown set, standing in for a PEC tangential
    // edge: A-DOFs are now a strict subset of the mesh's edges.
    std::vector<int> a_dof_edge;
    for (int e = 1; e < mesh.num_edges(); ++e) a_dof_edge.push_back(e);
    const int num_a = static_cast<int>(a_dof_edge.size());
    const int num_phi = 3;

    const GaugeIndexMap map = build_albanese_rubinacci_index_map(a_dof_edge, num_phi, tc);

    check(num_a < mesh.num_edges(), "coupled: A-DOFs are a strict subset of the mesh edges");

    int expected_cotree = 0;
    for (int e : a_dof_edge) {
        if (!tc.is_tree_edge[static_cast<std::size_t>(e)]) ++expected_cotree;
    }
    check(map.reduced_size == expected_cotree + num_phi,
          "coupled: reduced size == surviving A-DOFs + all Phi DOFs");

    bool phi_all_kept = true;
    for (int p = 0; p < num_phi; ++p) {
        if (map.full_to_reduced[static_cast<std::size_t>(num_a + p)] < 0) phi_all_kept = false;
    }
    check(phi_all_kept, "coupled: every Phi DOF survives (Phi is never gauged)");

    bool tree_a_all_dropped = true;
    for (int k = 0; k < num_a; ++k) {
        const bool is_tree = tc.is_tree_edge[static_cast<std::size_t>(a_dof_edge[static_cast<std::size_t>(k)])];
        const bool dropped = map.full_to_reduced[static_cast<std::size_t>(k)] < 0;
        if (is_tree != dropped) tree_a_all_dropped = false;
    }
    check(tree_a_all_dropped, "coupled: exactly the tree-edge A-DOFs are eliminated");

    // A complex coupled system with a recognisable value at every position,
    // so a misplaced row or column shows up as a wrong number rather than
    // just a wrong count.
    const int n_full = num_a + num_phi;
    SparseMatrixZ full(n_full, n_full);
    for (int r = 0; r < n_full; ++r) {
        for (int c = 0; c < n_full; ++c) {
            full.add(r, c, Complex(static_cast<double>(r + 1), static_cast<double>(c + 1)));
        }
    }
    full.compress();

    const SparseMatrixZ reduced = full.principal_submatrix(map.full_to_reduced, map.reduced_size);
    check(reduced.rows() == map.reduced_size && reduced.cols() == map.reduced_size,
          "coupled: reduced matrix is square with the mapped size");

    bool entries_match = true;
    for (int rr = 0; rr < map.reduced_size; ++rr) {
        for (int cc = 0; cc < map.reduced_size; ++cc) {
            const int r_full = map.reduced_to_full[static_cast<std::size_t>(rr)];
            const int c_full = map.reduced_to_full[static_cast<std::size_t>(cc)];
            if (!nearly(reduced.at(rr, cc), full.at(r_full, c_full))) entries_match = false;
        }
    }
    check(entries_match, "coupled: every surviving entry keeps its original value at its new position");

    // The Phi-Phi block must come through completely untouched -- it sits at
    // the bottom-right of both the full and the reduced system, since only
    // A-DOFs are ever removed and they all precede it.
    bool phi_block_intact = true;
    for (int p = 0; p < num_phi; ++p) {
        for (int q = 0; q < num_phi; ++q) {
            const int rr = map.full_to_reduced[static_cast<std::size_t>(num_a + p)];
            const int cc = map.full_to_reduced[static_cast<std::size_t>(num_a + q)];
            if (!nearly(reduced.at(rr, cc), full.at(num_a + p, num_a + q))) phi_block_intact = false;
        }
    }
    check(phi_block_intact, "coupled: the Phi-Phi block passes through unchanged");

    // Right-hand-side restriction and solution expansion round-trip, with
    // every eliminated entry coming back as exactly zero -- which for a
    // tree-edge A-DOF is the gauge condition a_t = 0 itself, not padding.
    std::vector<Complex> rhs_full(static_cast<std::size_t>(n_full));
    for (int i = 0; i < n_full; ++i) rhs_full[static_cast<std::size_t>(i)] = Complex(i + 1.0, -(i + 1.0));

    const std::vector<Complex> rhs_reduced = restrict_vector(rhs_full, map);
    check(static_cast<int>(rhs_reduced.size()) == map.reduced_size, "coupled: restricted RHS has the reduced length");

    const std::vector<Complex> expanded = expand_solution(rhs_reduced, map);
    check(expanded.size() == rhs_full.size(), "coupled: expanded solution has full length");

    bool round_trip_ok = true;
    bool eliminated_are_zero = true;
    for (int i = 0; i < n_full; ++i) {
        if (map.full_to_reduced[static_cast<std::size_t>(i)] >= 0) {
            if (!nearly(expanded[static_cast<std::size_t>(i)], rhs_full[static_cast<std::size_t>(i)])) round_trip_ok = false;
        } else if (expanded[static_cast<std::size_t>(i)] != Complex(0.0, 0.0)) {
            eliminated_are_zero = false;
        }
    }
    check(round_trip_ok, "coupled: restrict -> expand returns surviving entries unchanged");
    check(eliminated_are_zero, "coupled: eliminated tree-edge entries expand to exactly zero (a_t = 0)");
}

// --- Gauge correctness with n x A = 0 surfaces ------------------------------
//
// The boundary-first construction (`Claude outputs/
// tree_cotree_boundary_first_proposal.md`) exists because two simpler trees
// fail here, in opposite ways: a plain spanning tree ignoring the surfaces
// over-constrains (B wrong by 26-36 % on cube_4, with no error raised), and
// one root per surface under-constrains (a singular matrix). The two checks
// below catch both: nullity 0 rules out the second, and exact recovery of B
// for an arbitrary admissible field rules out the first.

Mesh make_kuhn_cube(int n) {
    Mesh m;
    const int s = n + 1;
    auto id = [s](int i, int j, int k) { return i + s * j + s * s * k; };
    for (int k = 0; k <= n; ++k)
        for (int j = 0; j <= n; ++j)
            for (int i = 0; i <= n; ++i) m.nodes.emplace_back(i, j, k);
    const int kuhn[6][4][3] = {{{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {1, 1, 1}}, {{0, 0, 0}, {1, 1, 0}, {0, 1, 0}, {1, 1, 1}},
                               {{0, 0, 0}, {0, 1, 0}, {0, 1, 1}, {1, 1, 1}}, {{0, 0, 0}, {0, 1, 1}, {0, 0, 1}, {1, 1, 1}},
                               {{0, 0, 0}, {0, 0, 1}, {1, 0, 1}, {1, 1, 1}}, {{0, 0, 0}, {1, 0, 1}, {1, 0, 0}, {1, 1, 1}}};
    for (int k = 0; k < n; ++k)
        for (int j = 0; j < n; ++j)
            for (int i = 0; i < n; ++i)
                for (const auto& t : kuhn) {
                    m.tets.push_back({id(i + t[0][0], j + t[0][1], k + t[0][2]), id(i + t[1][0], j + t[1][1], k + t[1][2]),
                                      id(i + t[2][0], j + t[2][1], k + t[2][2]), id(i + t[3][0], j + t[3][1], k + t[3][2])});
                }
    m.build_topology();
    return m;
}

std::vector<bool> plane_face_mask(const Mesh& m, int axis, double value) {
    std::vector<bool> mask(static_cast<std::size_t>(m.num_edges()), false);
    for (int f = 0; f < m.num_faces(); ++f) {
        if (!m.is_boundary_face(f)) continue;
        const auto& v = m.faces[static_cast<std::size_t>(f)];
        bool in_plane = true;
        for (int a = 0; a < 3; ++a) {
            const Vec3& p = m.nodes[static_cast<std::size_t>(v[a])];
            in_plane = in_plane && (axis == 0 ? p.x : (axis == 1 ? p.y : p.z)) == value;
        }
        if (!in_plane) continue;
        mask[static_cast<std::size_t>(m.find_edge(v[0], v[1]))] = true;
        mask[static_cast<std::size_t>(m.find_edge(v[1], v[2]))] = true;
        mask[static_cast<std::size_t>(m.find_edge(v[0], v[2]))] = true;
    }
    return mask;
}

// Rank deficiency by dense Gaussian elimination with complete pivoting. On
// these matrices the dropped pivots sit near 1e-14 and the kept ones above
// 0.05, so a relative threshold of 1e-10 separates them unambiguously.
int dense_nullity(std::vector<std::vector<double>> A) {
    const int n = static_cast<int>(A.size());
    if (n == 0) return 0;
    std::vector<double> pivots;
    for (int k = 0; k < n; ++k) {
        int pr = k, pc = k;
        double best = 0.0;
        for (int i = k; i < n; ++i)
            for (int j = k; j < n; ++j)
                if (std::abs(A[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)]) > best) {
                    best = std::abs(A[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)]);
                    pr = i;
                    pc = j;
                }
        pivots.push_back(best);
        if (best == 0.0) break;
        std::swap(A[static_cast<std::size_t>(k)], A[static_cast<std::size_t>(pr)]);
        for (auto& row : A) std::swap(row[static_cast<std::size_t>(k)], row[static_cast<std::size_t>(pc)]);
        for (int i = k + 1; i < n; ++i) {
            const double f = A[static_cast<std::size_t>(i)][static_cast<std::size_t>(k)] /
                             A[static_cast<std::size_t>(k)][static_cast<std::size_t>(k)];
            if (f == 0.0) continue;
            for (int j = k; j < n; ++j) {
                A[static_cast<std::size_t>(i)][static_cast<std::size_t>(j)] -=
                    f * A[static_cast<std::size_t>(k)][static_cast<std::size_t>(j)];
            }
        }
    }
    int nullity = n - static_cast<int>(pivots.size());
    for (double p : pivots) {
        if (p <= 1e-10 * pivots.front()) ++nullity;
    }
    return nullity;
}

void check_gauge_with_dirichlet(const Mesh& mesh, const std::vector<bool>& dirichlet, const std::string& label) {
    const int num_edges = mesh.num_edges();
    const TreeCotreeResult tc = build_tree_cotree(mesh, dirichlet);
    const SparseMatrix C = build_curl_matrix(mesh);
    const SparseMatrix M = C.transposed().multiply(C);

    int num_dirichlet = 0;
    for (bool d : dirichlet) num_dirichlet += d ? 1 : 0;
    const int expected_free = num_edges - num_dirichlet - tc.interior_tree_edge_count;

    // An arbitrary field that satisfies the boundary condition: zero on every
    // Dirichlet edge, anything elsewhere. j = M z is then a consistent
    // right-hand side, and any correct gauge must return a field with z's
    // curl -- the same B -- whatever A it picks.
    std::vector<double> z(static_cast<std::size_t>(num_edges), 0.0);
    for (int e = 0; e < num_edges; ++e) {
        if (!dirichlet[static_cast<std::size_t>(e)]) z[static_cast<std::size_t>(e)] = std::sin(1.7 * e + 0.3);
    }
    const std::vector<double> j = M.matvec(z);
    const std::vector<double> curl_z = C.matvec(z);
    auto relative_curl_error = [&](const std::vector<double>& a) {
        const std::vector<double> curl_a = C.matvec(a);
        double num = 0.0, den = 0.0;
        for (std::size_t f = 0; f < curl_z.size(); ++f) {
            num += (curl_a[f] - curl_z[f]) * (curl_a[f] - curl_z[f]);
            den += curl_z[f] * curl_z[f];
        }
        return std::sqrt(num / den);
    };

    // Method A.
    const GaugeVariant gauge_a = build_albanese_rubinacci_gauge(M, tc);
    check(gauge_a.reduced_matrix.rows() == expected_free,
          label + ": Method A keeps exactly edges - Dirichlet - interior tree unknowns");
    const auto a_dense = gauge_a.reduced_matrix.to_dense();
    const int nullity_a = dense_nullity(a_dense);
    check(nullity_a == 0, label + ": Method A reduced matrix has nullity 0 (gauge complete)");
    if (nullity_a != 0) return;  // singular: nothing further can be solved
    const std::vector<double> a_c = dense_solve(a_dense, select_cotree_entries(j, gauge_a.cotree_local_index));
    const std::vector<double> a_full = recover_albanese_rubinacci_solution(a_c, gauge_a, num_edges);
    check(relative_curl_error(a_full) < 1e-10,
          label + ": Method A recovers B exactly (gauge not over-constrained)");
    bool dirichlet_zero = true;
    for (int e = 0; e < num_edges; ++e) {
        if (dirichlet[static_cast<std::size_t>(e)] && a_full[static_cast<std::size_t>(e)] != 0.0) dirichlet_zero = false;
    }
    check(dirichlet_zero, label + ": Method A leaves every Dirichlet edge at exactly zero");

    // The B-recovery check must be able to FAIL: eliminate just one more free
    // edge -- a constraint that is not a gauge -- and B must come out wrong.
    // Without this, "error < 1e-10" could pass merely because the check is
    // blind. (The plain spanning tree this construction replaces
    // over-constrains the same way, with 37-43 extra edges on cube_4.)
    {
        GaugeIndexMap over;
        std::vector<int> identity(static_cast<std::size_t>(num_edges));
        for (int e = 0; e < num_edges; ++e) identity[static_cast<std::size_t>(e)] = e;
        over = build_albanese_rubinacci_index_map(identity, 0, tc);
        int extra = -1;
        for (int e = 0; e < num_edges && extra < 0; ++e) {
            if (over.full_to_reduced[static_cast<std::size_t>(e)] >= 0) extra = e;
        }
        std::vector<int> kept;
        for (int e = 0; e < num_edges; ++e) {
            if (over.full_to_reduced[static_cast<std::size_t>(e)] >= 0 && e != extra) kept.push_back(e);
        }
        std::vector<int> f2r(static_cast<std::size_t>(num_edges), -1);
        for (std::size_t k = 0; k < kept.size(); ++k) f2r[static_cast<std::size_t>(kept[k])] = static_cast<int>(k);
        const auto R = M.principal_submatrix(f2r, static_cast<int>(kept.size())).to_dense();
        std::vector<double> rhs(kept.size());
        for (std::size_t k = 0; k < kept.size(); ++k) rhs[k] = j[static_cast<std::size_t>(kept[k])];
        const std::vector<double> sol = dense_solve(R, rhs);
        std::vector<double> a_over(static_cast<std::size_t>(num_edges), 0.0);
        for (std::size_t k = 0; k < kept.size(); ++k) a_over[static_cast<std::size_t>(kept[k])] = sol[k];
        check(relative_curl_error(a_over) > 1e-3,
              label + ": negative control -- one extra eliminated edge visibly corrupts B");
    }

    // Method D.
    const EssentialIncidenceMatrix F = compute_essential_incidence_matrix(mesh, tc);
    check(F.num_cotree_edges == expected_free, label + ": Method D has the same unknowns as Method A");
    const GaugeVariant gauge_d = build_munteanu_unsymmetric_gauge(M, tc, F);
    const auto d_dense = gauge_d.reduced_matrix.to_dense();
    const int nullity_d = dense_nullity(d_dense);
    check(nullity_d == 0, label + ": Method D reduced matrix has nullity 0");
    if (nullity_d != 0) return;
    const std::vector<double> d_c = dense_solve(d_dense, select_cotree_entries(j, gauge_d.cotree_local_index));
    const std::vector<double> d_full = recover_munteanu_unsymmetric_solution(d_c, gauge_d, F, tc, num_edges);
    check(relative_curl_error(d_full) < 1e-10, label + ": Method D recovers B exactly");
}

void test_gauges_with_dirichlet_surfaces() {
    {
        const Mesh m = make_kuhn_cube(3);
        check_gauge_with_dirichlet(m, boundary_edge_mask(m), "cube n=3, n x A = 0 on the whole boundary");
    }
    {
        // Two separate surfaces: the case the previous one-root-per-body tree
        // left singular.
        const Mesh m = make_kuhn_cube(3);
        std::vector<bool> mask = plane_face_mask(m, 0, 0.0);
        const std::vector<bool> other = plane_face_mask(m, 0, 3.0);
        for (std::size_t e = 0; e < mask.size(); ++e) mask[e] = mask[e] || other[e];
        check_gauge_with_dirichlet(m, mask, "cube n=3, two separate Dirichlet faces");
    }
    {
        // One cell thick, both faces Dirichlet, interior edges joining them.
        const Mesh m = make_kuhn_cube(1);
        std::vector<bool> mask = plane_face_mask(m, 2, 0.0);
        const std::vector<bool> other = plane_face_mask(m, 2, 1.0);
        for (std::size_t e = 0; e < mask.size(); ++e) mask[e] = mask[e] || other[e];
        check_gauge_with_dirichlet(m, mask, "thin layer, both faces Dirichlet");
    }
}

}  // namespace

int main() {
    run_condition_number_analytic_checks();

    for (const auto& mesh_case : {std::pair<std::string, Mesh>{"single tet", make_single_tet()},
                                   std::pair<std::string, Mesh>{"two tets sharing a face", make_two_tets_sharing_a_face()}}) {
        const std::string& label = mesh_case.first;
        const Mesh& mesh = mesh_case.second;
        const int num_edges = mesh.num_edges();

        const TreeCotreeResult tc =
            build_tree_cotree(mesh, std::vector<bool>(static_cast<std::size_t>(num_edges), false));
        const EssentialIncidenceMatrix F = compute_essential_incidence_matrix(mesh, tc);

        check(F.num_free_groups == tc.interior_tree_edge_count,
              label + ": F's free-group count matches the interior tree edge count");
        check(F.num_cotree_edges == num_edges - tc.tree_edge_count, label + ": F's cotree count matches num_edges - tree_edge_count");

        const SparseMatrix M = build_test_M(mesh);
        const auto M_dense = M.to_dense();

        // --- Method A: reduced matrix must be EXACTLY M's principal
        //     submatrix on cotree rows/cols (Munteanu's own description:
        //     "eliminating rows and columns... which correspond to tree
        //     edges"). Checked directly against M_dense, independent of the
        //     implementation. ---
        const GaugeVariant gauge_a = build_albanese_rubinacci_gauge(M, tc);
        const auto a_dense = gauge_a.reduced_matrix.to_dense();
        bool submatrix_matches = true;
        for (int e1 = 0; e1 < num_edges; ++e1) {
            const int r = gauge_a.cotree_local_index[static_cast<std::size_t>(e1)];
            if (r == -1) continue;
            for (int e2 = 0; e2 < num_edges; ++e2) {
                const int c = gauge_a.cotree_local_index[static_cast<std::size_t>(e2)];
                if (c == -1) continue;
                if (!nearly(a_dense[static_cast<std::size_t>(r)][static_cast<std::size_t>(c)],
                            M_dense[static_cast<std::size_t>(e1)][static_cast<std::size_t>(e2)], 1e-9)) {
                    submatrix_matches = false;
                }
            }
        }
        check(submatrix_matches, label + ": Albanese-Rubinacci reduced matrix is exactly M's cotree principal submatrix");

        // --- Method D builds cleanly and has the right shape. ---
        const GaugeVariant gauge_d = build_munteanu_unsymmetric_gauge(M, tc, F);
        check(gauge_d.reduced_matrix.rows() == gauge_a.reduced_matrix.rows(),
              label + ": both gauge variants reduce to the same size");

        // --- The key physical cross-check: both gauge variants are valid
        //     gauge-fixings of the SAME curl-curl system, so for a
        //     manufactured, automatically-consistent right-hand side
        //     j = M*z, both should recover a field whose curl matches
        //     curl(z) exactly, even though the two recovered vector
        //     potentials a_A and a_D generally differ from each other and
        //     from z (different valid gauges, same physical B field). ---
        std::vector<double> z(static_cast<std::size_t>(num_edges));
        for (int e = 0; e < num_edges; ++e) z[static_cast<std::size_t>(e)] = 1.0 + 0.37 * static_cast<double>(e);
        const std::vector<double> j = dense_matvec(M_dense, z);

        const std::vector<double> j_c_a = select_cotree_entries(j, gauge_a.cotree_local_index);
        const std::vector<double> a_c_a = dense_solve(gauge_a.reduced_matrix.to_dense(), j_c_a);
        const std::vector<double> a_full_a = recover_albanese_rubinacci_solution(a_c_a, gauge_a, num_edges);

        const std::vector<double> j_c_d = select_cotree_entries(j, gauge_d.cotree_local_index);
        const std::vector<double> a_c_d = dense_solve(gauge_d.reduced_matrix.to_dense(), j_c_d);
        const std::vector<double> a_full_d = recover_munteanu_unsymmetric_solution(a_c_d, gauge_d, F, tc, num_edges);

        const SparseMatrix C = build_curl_matrix(mesh);
        const auto C_dense = C.to_dense();
        const std::vector<double> curl_a = dense_matvec(C_dense, a_full_a);
        const std::vector<double> curl_d = dense_matvec(C_dense, a_full_d);
        const std::vector<double> curl_z = dense_matvec(C_dense, z);

        bool curls_match_az = true, curls_match_ad = true;
        for (std::size_t f = 0; f < curl_z.size(); ++f) {
            if (!nearly(curl_a[f], curl_z[f], 1e-6)) curls_match_az = false;
            if (!nearly(curl_d[f], curl_z[f], 1e-6)) curls_match_ad = false;
        }
        check(curls_match_az, label + ": curl(a_AlbaneseRubinacci) matches curl(z) (same physical field)");
        check(curls_match_ad, label + ": curl(a_MunteanuUnsymmetric) matches curl(z) (same physical field)");

        // --- Condition-number estimates: both must be finite (the gauge
        //     genuinely removed the singularity) and >= 1 (a condition
        //     number can never be less than 1). Munteanu's own reported
        //     ordering (kappa_D <= kappa_A, i.e. the unsymmetric variant is
        //     at least as well-conditioned) is reported but not hard-
        //     asserted here -- our meshes are tiny (3x3 / 5x5 reduced
        //     systems) compared to the realistic problem sizes her own
        //     numerical tests used, so this is a methodology check, not a
        //     confirmation of her full ordering.
        const double kappa_a = estimate_condition_number(gauge_a.reduced_matrix);
        const double kappa_d = estimate_condition_number(gauge_d.reduced_matrix);
        check(std::isfinite(kappa_a) && kappa_a >= 1.0 - 1e-6, label + ": Albanese-Rubinacci condition number is finite and >= 1");
        check(std::isfinite(kappa_d) && kappa_d >= 1.0 - 1e-6, label + ": Munteanu-unsymmetric condition number is finite and >= 1");
        std::cout << label << ": kappa_A = " << kappa_a << ", kappa_D = " << kappa_d << "\n";
    }

    test_coupled_system_reduction();
    test_gauges_with_dirichlet_surfaces();

    std::cout << g_checks - g_failures << "/" << g_checks << " checks passed.\n";
    return g_failures == 0 ? 0 : 1;
}
