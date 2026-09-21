#pragma once

#include <stdexcept>
#include <string>

#include "aphi_solver/mesh.hpp"

namespace aphi_solver {

/// Thrown for any malformed/unsupported input to read_gmsh_msh.
class GmshReadError : public std::runtime_error {
public:
    explicit GmshReadError(const std::string& what) : std::runtime_error(what) {}
};

/// Reads an unstructured tetrahedral mesh from a Gmsh **ASCII** ".msh" file.
///
/// **Supported versions: 2.x and 4.1+.** 2.2 is the legacy format most
/// third-party exporters still emit; 4.1 is what current Gmsh writes by
/// default, so both are needed in practice. MSH 4.0 is rejected with a
/// message rather than guessed at -- it stores `$Nodes` with tags and
/// coordinates interleaved, where 4.1 separates them, and silently
/// misreading that would produce a plausible-looking but wrong mesh. Binary
/// files (`file_type != 0`) are rejected.
///
/// The two formats differ in more than layout, which is why this is a real
/// parser and not a reformat: in 2.2 an element's physical-group tag sits on
/// the element line, while in 4.x the element line carries only its geometric
/// *entity*, and the entity-to-physical-group mapping lives in `$Entities`.
/// A 4.x reader that skipped `$Entities` would return every element
/// untagged.
///
/// What is kept:
///  - Tetrahedra (element type 4) with their physical-group tag
///    (`Mesh::tet_tags`).
///  - Tagged surface triangles (element type 2) as
///    `Mesh::tagged_boundary_faces` -- how PEC walls, ports and the outer
///    truncation surface are marked.
///  - `$PhysicalNames`, if present, as `Mesh::physical_names`, so an input
///    file's declared region names can be checked against the mesh's own
///    rather than taken on trust.
/// Point, line and higher-order elements are read past and discarded.
///
/// Gmsh's node numbering is remapped to 0-based indices; node ids need not
/// be contiguous, 1-based, or ordered.
///
/// Before returning, calls `Mesh::build_topology()` (so edges, faces and
/// face->tet adjacency are ready) and then validates the result: every tet
/// must use four distinct nodes and be non-degenerate (tested as |V| / h^3
/// against the tet's own longest edge, so the check is scale-free), and no
/// face may be shared by more than two tets.
///
/// Throws GmshReadError on a missing file, an unsupported or missing
/// `$MeshFormat`, a malformed or unterminated section, a duplicate node id,
/// an element referencing an unknown node, or any validation failure above.
///
/// **Known gap (Sept 2026): this has never been run against a file produced
/// by Gmsh itself.** Every mesh the solver has read so far was synthesised
/// by this project's own tooling -- `tools/generate_cube_mesh.py`, the
/// fixtures in `tests/test_gmsh_reader.cpp`, or a hand conversion. Gmsh is
/// not installed on the development machine, so the format is implemented
/// from the specification.
///
/// The fixtures deliberately cover the structures real Gmsh output has that
/// a minimal spec-conforming file does not -- multi-block `$Nodes` (one
/// block per geometric entity), `parametric` coordinate lines carrying
/// extra u/u,v values, `$Entities` with trailing bounding-entity lists and
/// negative orientation tags, entities with zero and with multiple physical
/// tags, scientific-notation coordinates, and CRLF line endings. That
/// closes "does the implementation match its design". It cannot close "is
/// the design a correct reading of the format": only one real export can.
/// Exporting any mesh from Gmsh as `Version 2 ASCII` or `Version 4.1 ASCII`
/// and reading it is a worthwhile five-minute check before trusting this on
/// a real geometry -- it will either pass or fail loudly, since the
/// validation added here rejects a malformed result rather than proceeding.
Mesh read_gmsh_msh(const std::string& path);

}  // namespace aphi_solver
