// ---------------------------------------------------------------------------
// 08_MixedPort_Interior -- case 02's cylinder with its top cap INSIDE the box.
//
//     gmsh stub_rod.geo -3 -o stub_rod.msh
//
// WHY THIS GEOMETRY EXISTS
//
// Every port in cases 01-07 sits on the outer boundary, where `n x A = 0` pins
// the gauge function psi to zero and the terminal quantities are therefore
// gauge invariant -- measured to 1e-12 by 07_GaugeInvariance. That protection
// is what makes the suite pass, and docs/GAUGE_CHOICE.md Sec. 1.2 argues it
// does not extend to a port inside the domain.
//
// This is the minimal change that puts a port there. The wire is case 02's,
// unchanged in radius and material, but it stops at z = Lw instead of running
// the full height of the box, so:
//
//     wire_bottom  at z = 0   is on the box's bottom face -- ORDINARY boundary
//                             terminal, the 0 V reference;
//     wire_top     at z = Lw  is INSIDE the air -- conductor below, air above.
//                             A single-potential terminal there is the MIXED
//                             MATERIAL INTERIOR PORT, which is the one
//                             configuration tree-cotree has no rule for.
//
// Current enters the interior cap, flows down the rod, and leaves through the
// grounded cap on the boundary. That path is why this geometry works where an
// earlier attempt did not: a rod floating entirely inside the domain with a cut
// across it has NO return path, the problem is ill-posed, and the fields come
// out as garbage for a reason that has nothing to do with the gauge. The
// conductor-only control caught that, and the case was deleted. Here one
// terminal is still on the boundary, so the conduction path is closed and
// exactly ONE thing has changed against case 02/03.
//
// The control for this case is 03_Cylinder_1A_50Hz itself: same wire, same
// material, same drive, both caps on the boundary. Any difference is the
// interior mixed port and nothing else.
//
// GEOMETRY: a = 1.5 mm, box 40 mm cube, rod from z = 0 to z = 28 mm, so the
// interior cap sits 12 mm clear of the top face. Mesh sizing is case 02's,
// field for field.
// ---------------------------------------------------------------------------
SetFactory("Built-in");

a  = 1.5;     // wire radius, mm -- circumradius of the polygon
W  = 40.0;    // box side, mm
L  = 40.0;    // box height, mm
Lw = 28.0;    // WIRE height -- the one number that differs from case 02
N  = 40;      // sides of the wire polygon

lc_skin = 0.25;
lc_far  = 6.0;
d_far   = 4.0;
lc_core = 0.4;

// --- wire cross-section: a regular N-gon inscribed in radius a --------------
For i In {0:N-1}
  theta = 2*Pi*i/N;
  Point(100+i) = {a*Cos(theta), a*Sin(theta), 0, lc_skin};
EndFor
For i In {0:N-2}
  Line(100+i) = {100+i, 101+i};
EndFor
Line(100+N-1) = {100+N-1, 100};

Curve Loop(1) = {100:100+N-1};
Plane Surface(1) = {1};

// --- box cross-section: a square with the wire removed ----------------------
h = W/2;
Point(1) = {-h, -h, 0, lc_far};
Point(2) = { h, -h, 0, lc_far};
Point(3) = { h,  h, 0, lc_far};
Point(4) = {-h,  h, 0, lc_far};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Curve Loop(2) = {1, 2, 3, 4};
Plane Surface(2) = {2, 1};

// --- stage 1: wire and the air beside it, z = 0 .. Lw -----------------------
// Both in ONE call so the wire's lateral surface is created once and shared,
// which is what keeps the mesh conformal across the wire/air interface.
//
// Return layout per extruded surface: [0] = top, [1] = volume, then one entry
// per boundary curve. Surface 1 has N curves, so surface 2's entries start at
// N + 2.
ext1[] = Extrude {0, 0, Lw} { Surface{1, 2}; };

wire_volume   = ext1[1];
wire_top      = ext1[0];        // z = Lw -- INTERIOR once stage 2 sits on it
air_lower     = ext1[N + 3];
annulus_top   = ext1[N + 2];    // z = Lw, the air beside the wire

lateral[] = {};
For i In {0:N-1}
  lateral[] += ext1[2 + i];
EndFor

// --- stage 2: all air above the wire, z = Lw .. L ---------------------------
// Extruding BOTH top surfaces -- the wire's cap and the annulus -- continues
// the full square cross-section upward as air. Because stage 2 is built on
// stage 1's surfaces, the mesh stays conformal and `wire_top` becomes a genuine
// interior face with a tetrahedron on each side: conductor below, air above.
// That two-sided-ness is exactly what the port binding requires of an interior
// terminal, and what makes the material contrast ACROSS the port face real.
ext2[] = Extrude {0, 0, L - Lw} { Surface{wire_top, annulus_top}; };

air_above_cap = ext2[1];
air_above_ann = ext2[N + 3];

Physical Volume("wire", 1) = {wire_volume};
Physical Volume("air", 2)  = {air_lower, air_above_cap, air_above_ann};

Physical Surface("wire_bottom", 10) = {1};         // on the boundary: 0 V
Physical Surface("wire_top", 11)    = {wire_top};  // INTERIOR: the mixed port

// --- graded sizing: case 02's, unchanged ------------------------------------
Field[1] = Distance;
Field[1].SurfacesList = {lateral[]};
Field[1].Sampling = 100;

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = lc_skin;
Field[2].SizeMax = lc_far;
Field[2].DistMin = 0.0;
Field[2].DistMax = d_far;

Field[3] = MathEval;
Field[3].F = Sprintf("%g", lc_core);

Field[4] = Restrict;
Field[4].InField = 3;
Field[4].VolumesList = {wire_volume};

Field[5] = Min;
Field[5].FieldsList = {2, 4};

Background Field = 5;

Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;

Mesh.Algorithm3D = 1;
Mesh.Optimize = 1;
Mesh.OptimizeNetgen = 1;
Mesh.MshFileVersion = 4.1;
