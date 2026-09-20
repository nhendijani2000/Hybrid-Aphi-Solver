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

/// Reads an unstructured tetrahedral mesh from a Gmsh ASCII ".msh" file,
/// format version 2.2 (the legacy format documented at gmsh.info -- still
/// the simplest and most widely produced by other meshers' Gmsh exporters;
/// per docs/ROADMAP.md Phase 02 step 1, "a reasonable first target").
///
/// Only volume elements (Gmsh element type 4, the 4-node tetrahedron) are
/// kept; surface/line/point elements in the file (used for boundary tagging)
/// are read past and discarded here -- boundary/material tagging is a later
/// phase's concern (Phase 03 PEC handling, Phase 06 material regions), not
/// Phase 02's mesh-ingestion step. Gmsh's own 1-based node numbering is
/// remapped to this codebase's 0-based indexing; node numbers need not be
/// contiguous in the file.
///
/// Throws GmshReadError on a missing file, an unsupported/missing
/// $MeshFormat version, or a malformed $Nodes/$Elements section. Calls
/// Mesh::build_topology() before returning, so the result's edges/faces are
/// ready to use immediately.
Mesh read_gmsh_msh(const std::string& path);

}  // namespace aphi_solver
