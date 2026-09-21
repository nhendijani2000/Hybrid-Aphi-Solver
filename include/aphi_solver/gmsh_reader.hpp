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
/// **Validated against real Gmsh output (Sept 2026).** Both parsers were
/// written from the specification, which cannot by itself catch a
/// *misreading* of that specification -- so `tools/gmsh_validation_box.geo`
/// was meshed with Gmsh 4.13.1 and exported twice, as
/// `meshes/validation_box_v22.msh` and `meshes/validation_box_v41.msh`.
/// Both parse, and `tests/test_gmsh_reader.cpp` keeps them as permanent
/// fixtures: 338 nodes, 1110 tets, Euler characteristic 1, total volume
/// exactly the 1 mm box's 1e-9, no non-manifold face, CG = 0, every tet
/// tagged 7, all 180 tagged triangles resolving to real boundary faces, and
/// the two formats agreeing on every invariant.
///
/// That last point is the strongest part: a misreading would have to corrupt
/// 2.2 and 4.1 *identically* to escape the cross-check, despite the two
/// carrying physical tags by completely different mechanisms (on the element
/// line vs. through `$Entities`) and listing nodes in different order.
///
/// The hand-written fixtures alongside them remain useful for the structures
/// this particular geometry happens not to exercise -- `parametric`
/// coordinate lines, entities with zero or multiple physical tags, CRLF
/// endings, and the malformed inputs that must be rejected.
Mesh read_gmsh_msh(const std::string& path);

}  // namespace aphi_solver
