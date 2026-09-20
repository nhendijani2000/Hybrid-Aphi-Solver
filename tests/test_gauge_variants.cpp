// Minimal, dependency-free test runner -- see the comment in
// tests/CMakeLists.txt for why no external framework is fetched here.
//
// Covers docs/ROADMAP.md Phase 03 steps 2-3: the Albanese-Rubinacci and
// Munteanu-unsymmetric tree-cotree gauge variants (I. Munteanu, "Tree-cotree
// condensation properties"), built on top of Phase 03 step 1's spanning-tree
// decomposition, plus the crude power-iteration condition-number estimate.

#include <chrono>
#include <cmath>
#include <iostream>
#include <string>

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

// transpose(SparseMatrix) now lives in incidence.hpp/.cpp (promoted Sept
// 2026 once tools/compare_gauges.cpp needed the same operation) -- no local
// copy here anymore.

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
    return multiply(transpose(C), C);
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
    SparseMatrix A;
    A.rows = static_cast<int>(diag_values.size());
    A.cols = static_cast<int>(diag_values.size());
    for (std::size_t i = 0; i < diag_values.size(); ++i) {
        A.add(static_cast<int>(i), static_cast<int>(i), diag_values[i]);
    }
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

}  // namespace

int main() {
    run_condition_number_analytic_checks();

    for (const auto& mesh_case : {std::pair<std::string, Mesh>{"single tet", make_single_tet()},
                                   std::pair<std::string, Mesh>{"two tets sharing a face", make_two_tets_sharing_a_face()}}) {
        const std::string& label = mesh_case.first;
        const Mesh& mesh = mesh_case.second;
        const int num_edges = mesh.num_edges();

        std::vector<bool> is_pec(static_cast<std::size_t>(mesh.num_nodes()), false);
        const TreeCotreeResult tc = build_tree_cotree(mesh, is_pec);
        const EssentialIncidenceMatrix F = compute_essential_incidence_matrix(mesh, tc);

        check(F.num_free_groups == tc.tree_edge_count, label + ": F's free-group count matches the tree edge count");
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
        check(gauge_d.reduced_matrix.rows == gauge_a.reduced_matrix.rows,
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

    std::cout << g_checks - g_failures << "/" << g_checks << " checks passed.\n";
    return g_failures == 0 ? 0 : 1;
}
