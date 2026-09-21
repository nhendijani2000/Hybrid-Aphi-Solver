// A minimal 3-D geometry whose only purpose is to produce a mesh written by
// Gmsh ITSELF, so this project's reader can be checked against real output
// rather than against fixtures written from the specification.
//
// Why this is worth doing: as of Sept 2026 no file produced by Gmsh has ever
// been read by this solver (see the "Known gap" note in
// include/aphi_solver/gmsh_reader.hpp). Gmsh is not installed on the
// development machine, so both the 2.2 and 4.1 parsers were implemented from
// the spec. The test fixtures cover the structures real output has --
// multi-block $Nodes, parametric coordinates, $Entities bounding lists,
// CRLF -- but a fixture written from one's own reading of a format cannot
// catch a misreading of that format.
//
// The physical groups below are deliberate, not decoration. They are what
// exercises the path most likely to be wrong:
//   - In 2.2, an element's physical tag sits on the element line.
//   - In 4.1, it does NOT: the element line names a geometric entity, and
//     the entity -> physical-group mapping lives in $Entities.
// A mesh with no physical groups would leave every element untagged and
// would barely test that path at all. Both a volume group and surface groups
// are defined so tet tags and boundary-face tags are both covered, and the
// explicit numeric tags make the expected values easy to check by eye.
//
// ---------------------------------------------------------------------------
// Generate BOTH formats from the same geometry (the cross-check is stronger
// than either file alone -- the two readers should agree exactly):
//
//   gmsh gmsh_validation_box.geo -3 -format msh22 -o validation_box_v22.msh
//   gmsh gmsh_validation_box.geo -3 -format msh41 -o validation_box_v41.msh
//
// From the GUI instead: Mesh -> 3D, then File -> Export, choose "Mesh - Gmsh
// MSH (*.msh)", and in the dialog pick version 2 ASCII / 4.1 ASCII with
// "Save all elements" OFF (so only the physical groups above are written)
// and binary OFF. This reader rejects binary files by design.
//
// Expected, for a correct parse of either file:
//   - tet tags        all 7
//   - tagged faces    tags 20 and 21 only
//   - physical names  (3,7)="conductor", (2,20)="pec_wall", (2,21)="port"
//   - every tet non-degenerate, no face shared by more than two tets
//     (the reader validates all of this and throws rather than proceeding)
// ---------------------------------------------------------------------------

SetFactory("OpenCASCADE");

// A 1 mm cube -- units are arbitrary here; what matters is that the reader's
// degeneracy check is scale-free, so a small geometry is a useful test of
// that too.
Box(1) = {0, 0, 0, 1e-3, 1e-3, 1e-3};

// Volume 1 is the box; surfaces 1..6 are its faces.
Physical Volume("conductor", 7) = {1};
Physical Surface("pec_wall", 20) = {1};
Physical Surface("port", 21) = {2};

// Coarse on purpose: a few hundred tets is plenty to exercise every parsing
// path, and keeps the file small enough to inspect by hand if something
// disagrees.
Mesh.MeshSizeMax = 0.35e-3;
Mesh.MeshSizeMin = 0.15e-3;
