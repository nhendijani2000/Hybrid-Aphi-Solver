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
Mesh read_gmsh_msh(const std::string& path);

}  // namespace aphi_solver
