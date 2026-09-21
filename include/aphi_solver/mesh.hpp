#pragma once

#include <array>
#include <cstddef>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace aphi_solver {

/// A minimal 3-D point/vector type. Deliberately dependency-free, matching the
/// project's existing style (see complex_matrix.hpp).
struct Vec3 {
    double x = 0.0, y = 0.0, z = 0.0;

    Vec3() = default;
    Vec3(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}

    Vec3 operator+(const Vec3& o) const { return {x + o.x, y + o.y, z + o.z}; }
    Vec3 operator-(const Vec3& o) const { return {x - o.x, y - o.y, z - o.z}; }
    Vec3 operator*(double s) const { return {x * s, y * s, z * s}; }
    double dot(const Vec3& o) const { return x * o.x + y * o.y + z * o.z; }
    Vec3 cross(const Vec3& o) const {
        return {y * o.z - z * o.y, z * o.x - x * o.z, x * o.y - y * o.x};
    }
    double norm() const;
};

/// A tetrahedron, stored as 4 indices into Mesh::nodes. Local vertex order
/// 0..3 is arbitrary at read time (whatever the mesh file gives); everything
/// derived from it (local edge/face numbering below) is defined relative to
/// that order, per tet.
using TetVerts = std::array<int, 4>;

/// Local-vertex-pair numbering for a tet's 6 edges. This follows the same
/// edge order as the vertex labels 5..10 in J.-M. Jin, *The Finite Element
/// Method in Electromagnetics*, 3rd ed., 2014, Fig. 5.3 (already the source
/// for this project's P2 nodal Phi element, see docs/FORMULATION.md Sec 5.1) --
/// using the same order here means a tet's local edge index lines up directly
/// with which of Phi's 6 edge-midpoint nodes sits on that edge, which matters
/// once Phase 04 assembles both fields from the same per-tet loop.
/// In Jin's 1-indexed vertex labels this is (1,2),(1,3),(1,4),(2,3),(3,4),(2,4);
/// 0-indexed (used throughout this codebase) that is:
inline constexpr std::array<std::pair<int, int>, 6> kTetLocalEdgeVerts = {{
    {0, 1},  // local edge 0 -- Jin's edge-midpoint node 5
    {0, 2},  // local edge 1 -- node 6
    {0, 3},  // local edge 2 -- node 7
    {1, 2},  // local edge 3 -- node 8
    {2, 3},  // local edge 4 -- node 9
    {1, 3},  // local edge 5 -- node 10
}};

/// Local face numbering: face i is opposite local vertex i, with the other
/// three vertices listed in the cyclic order (i+1, i+2, i+3) mod 4. This is
/// the standard "face opposite vertex i" convention (e.g. Jin 2014's own
/// tetrahedron figures use it) and is purely a local-indexing choice, not
/// sourced from any one implementation.
inline constexpr std::array<std::array<int, 3>, 4> kTetLocalFaceVerts = {{
    {1, 2, 3},  // face 0, opposite vertex 0
    {2, 3, 0},  // face 1, opposite vertex 1
    {3, 0, 1},  // face 2, opposite vertex 2
    {0, 1, 2},  // face 3, opposite vertex 3
}};

/// A tetrahedral mesh with derived topology (globally unique edges and
/// faces, each in a fixed canonical orientation) needed to build the
/// discrete gradient (G) and curl (C) incidence matrices in incidence.hpp,
/// and to evaluate basis functions per tet in basis_functions.hpp.
///
/// Canonical orientation convention (used throughout, see docs/ROADMAP.md
/// Phase 02 step 2 and docs/FORMULATION.md): every global edge (i,j) is
/// stored with i < j, defining its positive direction as i -> j; every
/// global face (a,b,c) is stored with a < b < c, defining its positive
/// traversal as a -> b -> c -> a. These are the two facts incidence.hpp's
/// C and G are built from.
/// A boundary triangle that carried a tag in the mesh file (Gmsh elm-type
/// 2), used to mark PEC walls, ports and the outer truncation surface.
/// `nodes` is sorted ascending to match Mesh::faces' canonical orientation
/// (mesh.hpp's convention), so a tagged face can be looked up directly
/// against the derived topology via Mesh::find_face.
///
/// What a tag *means* -- copper, PEC, port 1 -- is the input file's job,
/// not the mesh's. Mesh carries integers only, which is what keeps
/// gmsh_reader.cpp format-agnostic and dependency-free.
struct TaggedFace {
    std::array<int, 3> nodes{};
    int tag = -1;
};

/// One entry of the mesh file's `$PhysicalNames` section: the human-readable
/// name a physical group was given in the .geo. Keyed by (dimension, tag)
/// rather than tag alone, because Gmsh's physical tags are only unique
/// *within* a dimension -- a volume group and a surface group can both be
/// tag 1, which is exactly the EDA case (region 1 = a conductor, boundary
/// 1 = a PEC wall).
///
/// This exists so the input file's declared region names can be checked
/// against the mesh's own instead of taken on trust: a JSON that says tag 1
/// is "via_conductor" while the mesh calls it "fr4_dielectric" is a
/// materials mix-up that would otherwise solve quietly and wrongly.
struct PhysicalName {
    int dimension = 0;  // 0 point, 1 curve, 2 surface, 3 volume
    int tag = -1;
    std::string name;
};

/// The tets sharing one global face. A conforming tetrahedral mesh gives
/// every face exactly 1 (boundary) or 2 (interior) -- `count` records the
/// true number even if it exceeds 2, so a non-manifold mesh can be detected
/// rather than silently truncated.
struct FaceTets {
    std::array<int, 2> tets{-1, -1};
    int count = 0;
};

class Mesh {
public:
    std::vector<Vec3> nodes;
    std::vector<TetVerts> tets;

    /// Physical-group tag of each tet, parallel to `tets` (Gmsh's
    /// convention: the first of an element's tags is its physical-group id).
    /// `build_topology` guarantees `tet_tags.size() == tets.size()`, filling
    /// -1 for a mesh built without tags (every hand-built test mesh, and any
    /// file whose elements carry no tags).
    std::vector<int> tet_tags;

    /// Tagged boundary triangles from the mesh file, in file order. Not
    /// derived from `tets` -- these come only from surface elements the mesh
    /// file actually contained, so this is empty for a mesh built in code or
    /// read from a file with no surface elements.
    std::vector<TaggedFace> tagged_boundary_faces;

    /// `$PhysicalNames` entries from the mesh file, if it had that section.
    /// Empty otherwise -- a mesh is perfectly usable without names.
    std::vector<PhysicalName> physical_names;

    /// For each global face, which tets contain it. Built by
    /// `build_topology` by inverting `tet_faces`, so it is always consistent
    /// with the derived topology rather than with whatever surface elements
    /// the file happened to carry.
    ///
    /// This is what makes "is this face on the boundary?" answerable at all
    /// (`is_boundary_face` below): a face of the *domain* boundary is one
    /// that exactly one tet touches. That is a geometric fact about the
    /// mesh, independent of tagging -- a mesh can have an untagged outer
    /// surface, and Phase 07's ABC needs the boundary regardless.
    std::vector<FaceTets> face_tets;

    /// Globally unique edges, each (i, j) with i < j. Index into this vector
    /// is the edge's global DOF/row index used by incidence.hpp and
    /// basis_functions.hpp.
    std::vector<std::pair<int, int>> edges;

    /// Globally unique faces, each (a, b, c) with a < b < c.
    std::vector<std::array<int, 3>> faces;

    /// For each tet, the 6 global edge indices in kTetLocalEdgeVerts order.
    std::vector<std::array<int, 6>> tet_edges;

    /// For each tet, the orientation of each of its 6 local edges relative to
    /// the corresponding global edge's canonical direction: +1 if the local
    /// direction (tets[t][vi] -> tets[t][vj]) agrees with the canonical
    /// low-index -> high-index direction, -1 if it is reversed.
    ///
    /// This is the `mEdgeSign` of the user's prior 3dedyaphi implementation
    /// (docs/FORMULATION.md Sec 5.1). It is NOT optional bookkeeping: a tet's
    /// local vertex order comes from the mesh file and is arbitrary, so on a
    /// real mesh roughly 58% of (tet, local edge) pairs are reversed (measured
    /// on meshes/cube_*.msh). The Whitney edge function built from local
    /// vertices vi, vj is the *negative* of the one built from the canonical
    /// pair whenever that happens, so any assembly that scatters a local edge
    /// contribution into global edge DOF tet_edges[t][le] must multiply it by
    /// this sign -- see whitney_edge_value_global in basis_functions.hpp.
    ///
    /// Only A's edge functions need this. Phi's P2 nodal functions do not: the
    /// edge-midpoint shape function 4*L_v0*L_v1 is symmetric in v0 and v1, so
    /// it is independent of which way the edge is traversed.
    ///
    /// signed char rather than int: 6 bytes per tet instead of 24, in a array
    /// read once per tet in the assembly loop right beside tet_edges.
    std::vector<std::array<signed char, 6>> tet_edge_signs;

    /// For each tet, the 4 global face indices in kTetLocalFaceVerts order.
    std::vector<std::array<int, 4>> tet_faces;

    /// Builds edges, faces, tet_edges, tet_faces from `nodes` and `tets`.
    /// Must be called once after `tets` is populated (the Gmsh reader in
    /// gmsh_reader.hpp calls this itself).
    void build_topology();

    /// Signed volume of tet t (positive iff vertices 0,1,2,3 are a
    /// right-handed / positively-oriented ordering). Zero or near-zero means
    /// a degenerate tet.
    double signed_tet_volume(int t) const;

    int num_nodes() const { return static_cast<int>(nodes.size()); }
    int num_edges() const { return static_cast<int>(edges.size()); }
    int num_faces() const { return static_cast<int>(faces.size()); }
    int num_tets() const { return static_cast<int>(tets.size()); }

    /// Returns the global edge index for the (unordered) vertex pair {i, j},
    /// or -1 if no such edge exists in this mesh. O(log num_edges) via an
    /// internal lookup map built by build_topology -- used by
    /// incidence.hpp's build_curl_matrix, so this needs to be more than a
    /// linear scan once meshes stop being test-sized.
    int find_edge(int i, int j) const;

    /// Returns the global face index for the (unordered) vertex triple
    /// {a, b, c}, or -1 if no such face exists in this mesh. Same cost and
    /// same intended use as find_edge: occasional point queries, not a
    /// per-entity traversal. Needed to resolve a `tagged_boundary_faces`
    /// entry against the topology derived from the tets -- a tagged triangle
    /// that resolves to -1 is not a face of any tet in the volume mesh.
    int find_face(int a, int b, int c) const;

    /// True iff global face `f` lies on the domain boundary, i.e. exactly
    /// one tet contains it.
    bool is_boundary_face(int f) const {
        return face_tets[static_cast<std::size_t>(f)].count == 1;
    }

    /// The name given to physical group `tag` of dimension `dimension`, or
    /// an empty string if the mesh carried no `$PhysicalNames` entry for it.
    /// A linear scan: there are a handful of physical groups even in a large
    /// model, and this is called once per region at load time, never per
    /// mesh entity.
    std::string physical_name(int dimension, int tag) const;

private:
    std::map<std::pair<int, int>, int> edge_lookup_;
    std::map<std::array<int, 3>, int> face_lookup_;
};

}  // namespace aphi_solver
