#pragma once

#include <vector>

#include "aphi_solver/mesh.hpp"

namespace aphi_solver {

/// A minimal sparse real matrix in COO (triplet) form: (row, col, value)
/// entries, duplicates summed. This is what docs/ROADMAP.md Phase 02 step 2
/// means by "sparse, entry-per-orientation, cheap to build" for the
/// incidence matrices below -- it is not meant to be a general-purpose
/// sparse-linear-algebra type (Phase 04/05 will need a real sparse solver
/// input format; this one exists to make G, C, and the CG=0 check cheap and
/// easy to verify, nothing more).
struct SparseMatrix {
    struct Entry {
        int row;
        int col;
        double value;
    };

    int rows = 0;
    int cols = 0;
    std::vector<Entry> entries;

    void add(int row, int col, double value) { entries.push_back({row, col, value}); }

    /// Sums duplicate (row, col) entries in place. G and C never produce
    /// duplicates by construction (see build_gradient_matrix /
    /// build_curl_matrix), but a product of two SparseMatrix values can, so
    /// `multiply` below calls this before returning.
    void coalesce();

    /// Dense conversion, for small test meshes / debugging only -- never use
    /// this on a production-size mesh.
    std::vector<std::vector<double>> to_dense() const;
};

/// Naive sparse-times-sparse product `A * B`. Not optimized (this project's
/// real assembly will need a proper sparse solver, per Phase 05); this exists
/// so the CG = 0 identity can be checked directly and cheaply on
/// Phase-02-scale test meshes.
SparseMatrix multiply(const SparseMatrix& a, const SparseMatrix& b);

/// Returns the transpose of `m` -- a new SparseMatrix with every (row, col,
/// value) entry remapped to (col, row, value). O(nnz), and coalesced before
/// returning (transposing alone can't create duplicates, but keeping this
/// consistent with `multiply`'s contract costs nothing and avoids surprises
/// for a caller who chains the two). Promoted here (Sept 2026) from a
/// test-local helper in tests/test_gauge_variants.cpp once a second
/// consumer -- tools/compare_gauges.cpp -- needed the same operation, per
/// this project's practice of sharing rather than duplicating once
/// something is used in more than one place.
SparseMatrix transpose(const SparseMatrix& m);

/// True iff every entry of `m` has magnitude <= tol (an all-zero matrix,
/// possibly with explicit but negligible entries from cancellation in
/// `multiply`).
bool is_zero_matrix(const SparseMatrix& m, double tol = 1e-9);

/// Builds the discrete gradient operator G : (node values) -> (edge values),
/// size num_edges x num_nodes. For edge e = (i, j) with i < j (Mesh's
/// canonical edge orientation, see mesh.hpp), row e has G(e, i) = -1 and
/// G(e, j) = +1, so (G * phi)(e) = phi(j) - phi(i) -- the standard
/// finite-difference "gradient along the edge" used to define tree-cotree's
/// node-edge incidence (docs/ROADMAP.md Phase 02 step 2; Munteanu's
/// convention, already cited in this project).
SparseMatrix build_gradient_matrix(const Mesh& mesh);

/// Builds the discrete curl operator C : (edge values) -> (face values),
/// size num_faces x num_edges. For face f = (a, b, c) with a < b < c (Mesh's
/// canonical face orientation), row f has C(f, edge(a,b)) = +1,
/// C(f, edge(b,c)) = +1, C(f, edge(a,c)) = -1 -- the signed sum of the
/// face's three boundary edges traversed a -> b -> c -> a. This is the
/// standard discrete Stokes / boundary operator for a triangulated surface
/// (e.g. Bossavit 1998 Ch. 5, already cited in docs/REFERENCES.md), not
/// specific to any one implementation.
SparseMatrix build_curl_matrix(const Mesh& mesh);

}  // namespace aphi_solver
