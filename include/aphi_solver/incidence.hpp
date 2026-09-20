#pragma once

#include "aphi_solver/mesh.hpp"
#include "aphi_solver/sparse_matrix.hpp"

namespace aphi_solver {

/// The real-valued sparse type the incidence operators use. Through Phases
/// 02-03 this was a separate COO struct with a `std::map`-based `coalesce()`
/// and free `multiply`/`transpose`/`is_zero_matrix` helpers; it is now an
/// alias for `Sparse<double>` (`sparse_matrix.hpp`), whose CSR form the
/// assembled system needs anyway -- see `docs/ROADMAP.md` Phase 03.5, step
/// 2. The old free functions became member functions in the move:
/// `multiply(a, b)` -> `a.multiply(b)`, `transpose(m)` -> `m.transposed()`,
/// `is_zero_matrix(m)` -> `m.is_zero()`, `coalesce()` -> `compress()`, and
/// the `rows`/`cols` fields are now `rows()`/`cols()` accessors.
using SparseMatrix = SparseMatrixD;

/// Builds the discrete gradient operator G : (node values) -> (edge values),
/// size num_edges x num_nodes. For edge e = (i, j) with i < j (Mesh's
/// canonical edge orientation, see mesh.hpp), row e has G(e, i) = -1 and
/// G(e, j) = +1, so (G * phi)(e) = phi(j) - phi(i) -- the standard
/// finite-difference "gradient along the edge" used to define tree-cotree's
/// node-edge incidence (docs/ROADMAP.md Phase 02 step 2; Munteanu's
/// convention, already cited in this project).
///
/// Returned already compressed, so the caller can use CSR operations
/// directly.
SparseMatrix build_gradient_matrix(const Mesh& mesh);

/// Builds the discrete curl operator C : (edge values) -> (face values),
/// size num_faces x num_edges. For face f = (a, b, c) with a < b < c (Mesh's
/// canonical face orientation), row f has C(f, edge(a,b)) = +1,
/// C(f, edge(b,c)) = +1, C(f, edge(a,c)) = -1 -- the signed sum of the
/// face's three boundary edges traversed a -> b -> c -> a. This is the
/// standard discrete Stokes / boundary operator for a triangulated surface
/// (e.g. Bossavit 1998 Ch. 5, already cited in docs/REFERENCES.md), not
/// specific to any one implementation.
///
/// Note that C and G are built from the mesh's *canonical* global edge and
/// face orientations, independently of any tet's local vertex ordering. The
/// per-tet basis functions in basis_functions.hpp are the place where local
/// ordering enters and must be reconciled with these conventions via
/// Mesh::tet_edge_signs -- see whitney_edge_value_global there.
///
/// Returned already compressed, so the caller can use CSR operations
/// directly.
SparseMatrix build_curl_matrix(const Mesh& mesh);

}  // namespace aphi_solver
