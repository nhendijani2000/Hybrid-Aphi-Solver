// The Phase 01 EDA target problem, as a meshable geometry.
//
// docs/FORMULATION.md Sec. 7.1: a straight cylindrical via of radius
// a = 0.1 mm on-axis inside a cylindrical return/shield conductor of radius
// b = 0.5 mm, length L = 2 mm, filled with a lossy dielectric representative
// of an EDA package (eps_r = 4.3, tan(delta) = 0.02). This is the geometry
// Phase 04/05 must reproduce Z_0, L', C' and R_DC for:
//
//   Z_0  = (1/(2*pi)) * sqrt(mu/eps) * ln(b/a)
//   L'   = (mu/(2*pi)) * ln(b/a)
//   C'   = (2*pi*eps) / ln(b/a)
//   R_DC = 1 / (sigma_cond * pi * a^2)
//
// Both conductor and dielectric are meshed as VOLUMES rather than the
// dielectric alone with a PEC inner boundary: R_DC is a property of current
// flowing in the via's cross-section, so the conductor needs volume DOFs to
// carry it. The outer shield is a PEC boundary condition rather than a
// meshed volume -- the return path's own resistance is not one of the
// closed-form checks above, and meshing a shell would add elements that
// earn nothing at this stage.
//
// ---------------------------------------------------------------------------
//   gmsh coax_via.geo -3 -format msh41 -o ../meshes/coax_via.msh
//
// Entity numbering below was NOT guessed: the geometry was meshed once
// without physical groups and the tags read back out of the resulting 4.1
// $Entities section (boolean operations renumber entities, so guessing which
// surface is the outer wall is exactly how the wrong boundary gets tagged
// silently). If a future Gmsh renumbers them, the per-region volume check in
// the reader test catches it -- the meshed volume of each physical group is
// compared against pi*a^2*L and pi*(b^2-a^2)*L, which no mislabelling
// survives.
// ---------------------------------------------------------------------------

SetFactory("OpenCASCADE");

a = 0.1e-3;   // inner conductor (via) radius
b = 0.5e-3;   // shield inner radius
L = 2.0e-3;   // length

Cylinder(1) = {0, 0, 0, 0, 0, L, a};
Cylinder(2) = {0, 0, 0, 0, 0, L, b};

// Fragment so the two volumes share their r = a interface conformally.
// Without it each would mesh independently and the interface nodes would not
// match, which breaks tangential continuity across the material boundary --
// the very thing the A-Phi formulation needs to resolve there.
BooleanFragments{ Volume{1}; Delete; }{ Volume{2}; Delete; }

// --- Regions ---------------------------------------------------------------
// Volume 1: r <= a, the via.          Volume 2: a <= r <= b, the dielectric.
Physical Volume("via_conductor", 1) = {1};
Physical Volume("dielectric", 2) = {2};

// --- Boundaries ------------------------------------------------------------
// Surface 4 is the r = b cylindrical wall: the shield, a PEC boundary.
Physical Surface("shield_pec", 10) = {4};

// z = L: the via is shorted to the return conductor there (Sec. 7.1), so
// both end faces at that plane are one conducting termination.
Physical Surface("short_end", 11) = {2, 5};

// z = 0: the feed. The via's own end face and the dielectric's annular end
// face are tagged separately, because the lumped current source drives the
// former against the latter -- Phase 06 replaces this with a real port, but
// keeping them distinct now means that change is a reinterpretation of
// existing tags rather than a re-mesh.
Physical Surface("port_conductor", 20) = {3};
Physical Surface("port_return", 21) = {6};

// --- Mesh sizing -----------------------------------------------------------
// The via radius (0.1 mm) is the smallest feature and has to be resolved
// across, not just along: too coarse here and the conductor is a few
// elements wide, which will not carry a credible current distribution for
// R_DC. Refined on the conductor, coarser out in the dielectric where the
// fields vary slowly.
Field[1] = Distance;
Field[1].SurfacesList = {1};
Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = 0.035e-3;
Field[2].SizeMax = 0.13e-3;
Field[2].DistMin = 0.05e-3;
Field[2].DistMax = 0.30e-3;
Background Field = 2;

Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;
