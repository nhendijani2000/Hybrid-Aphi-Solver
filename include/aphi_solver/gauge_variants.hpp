#pragma once

#include <string>
#include <vector>

#include "aphi_solver/incidence.hpp"
#include "aphi_solver/mesh.hpp"
#include "aphi_solver/tree_cotree.hpp"

namespace aphi_solver {

/// Solves the dense linear system A x = b via Gaussian elimination with
/// partial pivoting. Generic (any size), used by the gauge-variant and
/// condition-number code below at Phase-02/03-scale problem sizes (small,
/// correctness-focused test meshes -- per docs/ENGINEERING_STANDARDS.md,
/// this is the right tool here; a real sparse factorization/iterative
/// solver is Phase 04/05's job, once assembly produces production-size
/// systems). Throws std::runtime_error if A is numerically singular.
std::vector<double> dense_solve(std::vector<std::vector<double>> A, std::vector<double> b);

/// The "essential incidence matrix" F = G_c * G_t^{-1} in I. Munteanu,
/// "Tree-cotree condensation properties" (Sec. III.A) -- G here being the
/// discrete gradient/node-edge incidence matrix restricted to the *free*
/// DOF groups from a TreeCotreeResult (a reference/PEC-body group's
/// potential is fixed, not a free unknown, so its column is simply dropped
/// rather than included), and split by rows into its cotree part G_c and
/// (square, invertible) tree part G_t after reordering edges cotree-then-
/// tree. F is what expresses each tree edge's DOF as a linear combination
/// of cotree DOFs (a_t = -F^T a_c), the basis of Munteanu's "unsymmetric"
/// gauge variant below.
///
/// Built via an O(V) tree back-substitution rather than a general dense
/// matrix inversion: G_t is exactly the spanning tree's own node-incidence
/// structure (TreeCotreeResult::parent_group / discovering_edge already
/// records it), so solving G_t x = b for any b is a single top-down walk
/// of the tree in group-discovery order, accumulating each group's value
/// from its parent's -- not an O(n^3) dense elimination. See
/// docs/ENGINEERING_STANDARDS.md.
struct EssentialIncidenceMatrix {
    int num_cotree_edges = 0;
    int num_free_groups = 0;

    /// Dense [cotree_local_index][free_group_index] values -- dense because
    /// F is not generally sparse even though G_c/G_t are (each row of F is
    /// a cumulative sum along a tree path, which can touch many free
    /// groups); fine at this phase's test-mesh scale.
    std::vector<std::vector<double>> values;

    /// Global edge index -> compact cotree-local row index, or -1 for a
    /// tree edge.
    std::vector<int> cotree_local_index;

    /// Group id -> compact free-group column index, or -1 for a reference
    /// group.
    std::vector<int> free_group_index;
};

EssentialIncidenceMatrix compute_essential_incidence_matrix(const Mesh& mesh, const TreeCotreeResult& tc);

/// One gauge variant's reduced system and the means to recover the full
/// edge-DOF solution from it.
struct GaugeVariant {
    std::string name;

    /// The reduced (num_cotree x num_cotree) system matrix.
    SparseMatrix reduced_matrix;

    /// Global edge index -> compact cotree-local index, or -1 for a tree
    /// edge (same convention as EssentialIncidenceMatrix::cotree_local_index,
    /// repeated here so a GaugeVariant is self-contained).
    std::vector<int> cotree_local_index;
};

/// Restricts a full E-length vector to its cotree entries, in
/// cotree_local_index order -- the right-hand-side reduction j_c = j[cotree]
/// that both gauge variants below share (both use the same plain
/// cotree-row test space, per Munteanu Sec. V).
std::vector<double> select_cotree_entries(const std::vector<double>& full, const std::vector<int>& cotree_local_index);

/// Which rows/columns of a system the tree-cotree gauge eliminates, and
/// where the survivors land. Deliberately holds nothing but indices: the
/// same map applies to a real or a complex matrix, so this is what lets one
/// implementation of the reduction serve both (`Sparse<T>::
/// principal_submatrix`). See docs/ROADMAP.md Phase 03.5, step 4.
struct GaugeIndexMap {
    /// full_to_reduced[i] == -1 if system index i is eliminated, otherwise
    /// its index in the reduced system.
    std::vector<int> full_to_reduced;

    /// The inverse: reduced_to_full[k] is the full-system index of reduced
    /// index k. Used to scatter a reduced solution back out.
    std::vector<int> reduced_to_full;

    int reduced_size = 0;
};

/// Builds the Albanese-Rubinacci index map for a coupled [a; Phi] system,
/// laid out as `assemble_sparse` produces it: A-DOFs first in
/// [0, a_dof_edge.size()), then Phi DOFs.
///
/// `a_dof_edge[k]` is the mesh edge that A-DOF k represents. In the
/// full-wave regime that is the identity (N_A == N_edges,
/// docs/FORMULATION.md Sec 5.2), but it becomes a genuine mapping once PEC
/// tangential edges are removed from the unknown set (Sec 5.4) or the
/// reduced low-frequency regime confines A to Omega_c -- which is precisely
/// why the gauge reduction cannot keep assuming its index space *is* the
/// mesh edge index space.
///
/// Every tree-edge A-DOF is eliminated (a_t = 0); every cotree A-DOF and
/// every Phi DOF is kept. Phi is never gauged: the gauge freedom being
/// removed is A's alone (docs/FORMULATION.md Sec 5.3).
GaugeIndexMap build_albanese_rubinacci_index_map(const std::vector<int>& a_dof_edge, int num_phi_dofs,
                                                   const TreeCotreeResult& tc);

/// Restricts a full-system right-hand side to the kept indices.
template <typename T>
std::vector<T> restrict_vector(const std::vector<T>& full, const GaugeIndexMap& map) {
    std::vector<T> out(static_cast<std::size_t>(map.reduced_size), T{});
    for (int k = 0; k < map.reduced_size; ++k) {
        out[static_cast<std::size_t>(k)] = full[static_cast<std::size_t>(map.reduced_to_full[static_cast<std::size_t>(k)])];
    }
    return out;
}

/// Scatters a reduced solution back to full-system length, leaving every
/// eliminated entry at exactly zero -- which for a tree-edge A-DOF is not a
/// placeholder but the gauge condition itself, a_t = 0.
template <typename T>
std::vector<T> expand_solution(const std::vector<T>& reduced, const GaugeIndexMap& map) {
    std::vector<T> out(map.full_to_reduced.size(), T{});
    for (int k = 0; k < map.reduced_size; ++k) {
        out[static_cast<std::size_t>(map.reduced_to_full[static_cast<std::size_t>(k)])] =
            reduced[static_cast<std::size_t>(k)];
    }
    return out;
}

/// Method A -- Albanese & Rubinacci's gauge (Munteanu's "Method A"):
/// literally sets every tree-edge DOF to zero, a_t = 0, and solves the
/// principal submatrix of the full curl-curl system M restricted to
/// cotree rows and columns (Munteanu, Sec. IV: "equivalent to eliminating
/// those rows and columns... which correspond to the tree edges"). This is
/// a genuine submatrix of M, so it preserves M's sparsity pattern exactly
/// (no fill-in) -- the simplest of the compared variants, though not the
/// best-conditioned.
GaugeVariant build_albanese_rubinacci_gauge(const SparseMatrix& M, const TreeCotreeResult& tc);

/// Recovers the full E-length solution from an Albanese-Rubinacci reduced
/// solution: cotree entries as solved, every tree entry exactly zero.
std::vector<double> recover_albanese_rubinacci_solution(const std::vector<double>& a_c, const GaugeVariant& variant,
                                                          int num_edges);

/// Method D -- Munteanu's "unsymmetric" gauge variant: rather than zeroing
/// the tree edges, expresses them as a_t = -F^T a_c via the essential
/// incidence matrix F above, and reduces the full system by an oblique
/// (Petrov-Galerkin) projection -- test with the plain cotree-row selector,
/// trial with a = L^T a_c (L = [I | -F]) -- derived here from Munteanu's
/// own stated general framework (Sec. V: a reduced system
/// W^T A V y = W^T f for test/trial bases W, V) together with her explicit
/// statement of which bases this method uses, rather than transcribed
/// directly from her equation (19): the OCR/extraction of that specific
/// equation reused the symbol "F" for two different things elsewhere in
/// the paper (both "number of mesh faces" and "the essential incidence
/// matrix"), which made copying it verbatim too risky to trust. The
/// resulting reduced matrix is generally NOT symmetric even though M is --
/// exactly the reason Munteanu calls this variant "unsymmetric", and a
/// useful independent check that this derivation lines up with her own
/// naming. The full step-by-step algebra (G_t invertibility, the
/// G^T a = 0 <=> a_t = -F^T a_c equivalence, a from-scratch transversality
/// proof that this gauge choice is valid, and the block-matrix reduction
/// itself) is written out in full in `docs/TREE_COTREE_GAUGE.md` -- this
/// comment summarizes it, but that document is the actual derivation.
GaugeVariant build_munteanu_unsymmetric_gauge(const SparseMatrix& M, const TreeCotreeResult& tc,
                                               const EssentialIncidenceMatrix& F);

/// Recovers the full E-length solution from a Munteanu-unsymmetric reduced
/// solution: cotree entries as solved, tree entries via a_t = -F^T a_c.
std::vector<double> recover_munteanu_unsymmetric_solution(const std::vector<double>& a_c, const GaugeVariant& variant,
                                                            const EssentialIncidenceMatrix& F,
                                                            const TreeCotreeResult& tc, int num_edges);

/// A crude condition-number ESTIMATE (largest singular value / smallest
/// singular value) via power iteration on A^T*A (largest eigenvalue) and
/// inverse power iteration on the same matrix (smallest eigenvalue) --
/// exactly the "even a crude power-iteration estimate" docs/ROADMAP.md
/// Phase 03 step 3 asks for, not a full SVD. Using A^T*A rather than A
/// directly is what makes this correct for BOTH gauge variants uniformly:
/// Method A's reduced matrix is symmetric (a principal submatrix of the SPD
/// M), but Method D's is explicitly not, so its eigenvalues aren't its
/// singular values and power iteration on A itself would not give a
/// meaningful condition number for it.
///
/// Sparse throughout: A^T*A is never formed (it can be substantially denser
/// than A -- see docs/TREE_COTREE_GAUGE.md Sec. 4-5 on Method D's fill-in),
/// every matvec is O(nnz), and the inverse-power-iteration step solves with
/// Conjugate Gradient (O(nnz) per CG iteration) instead of a dense O(n^3)
/// Gaussian elimination from scratch every step. This was originally a
/// dense implementation scoped to Phase 03's tiny validation meshes, and
/// was brought forward to sparse ahead of its planned Phase 04/05 slot
/// specifically so this estimator (and the gauge-comparison CLI built on
/// it) can run on realistic mesh sizes -- see docs/ENGINEERING_STANDARDS.md.
double estimate_condition_number(const SparseMatrix& A, int max_iterations = 500, double tol = 1e-12);

}  // namespace aphi_solver
