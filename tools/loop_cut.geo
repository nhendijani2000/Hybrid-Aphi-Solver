// A conducting ring driven through an internal cut -- the second test geometry.
// See examples/loop_internal_port.aphi.
//
// Meshed with **Gmsh 4.13.1**, as meshes/cylinder_box.msh was. A different
// version can renumber entities and size elements differently, which would
// silently change every count a test asserts. Regenerate with:
//
//     gmsh -3 tools/loop_cut.geo -o meshes/loop_cut.msh
//
// ---------------------------------------------------------------------------
// The shape: a washer, not a torus.
//
// An annulus in the x-y plane extruded along z is a closed conducting loop --
// current circulates azimuthally around the hole on the z axis. That lets the
// whole thing be built the way tools/cylinder_box.geo builds its wire: a 2D
// cross-section, then one Extrude. A revolved circular torus would need the
// OpenCASCADE kernel and a boolean to get the cut, for no gain here.
//
// Why the ring is POLYGONAL, inside and out: the same reason the cylinder's
// wire is (see that file). Planar lateral faces mean every surface triangle
// on them has a horizontal normal whatever the triangulation does, so the
// natural condition dPhi/dn = 0 holds there exactly. N must be EVEN, so that
// theta = 0 and theta = pi are both vertices and the two halves below meet
// exactly.
//
// ---------------------------------------------------------------------------
// How the internal cut is made.
//
// A port surface needs a tet on each side, so the cut cannot simply be drawn
// inside a volume -- it has to BE the interface between two volumes. So the
// annulus is built as two half-annuli, y >= 0 and y <= 0, sharing the two
// radial lines at theta = 0 and theta = pi. Extruding both in one call makes
// the surfaces over those shared lines once and shares them, which is what
// makes the mesh conformal across the cut and gives it tets on both sides.
//
// `loop_cut` is the one at theta = 0, the half-plane y = 0, x > 0 -- which is
// where examples/loop_internal_port.aphi says it is, and why the hint `+y`
// there reads as "counterclockwise seen from +z".
//
// The ring sits fully inside the air box, touching no outer face: the loop is
// floating, and its own cut is its only potential reference.
// ---------------------------------------------------------------------------

SetFactory("Built-in");

Ri = 0.6;    // ring inner radius, mm (circumradius of the inner polygon)
Ro = 1.0;    // ring outer radius, mm
t  = 0.3;    // ring thickness along z, mm
N  = 24;     // sides of the ring polygons -- MUST be even

Wxy = 4.0;   // air box side in x and y, mm
Wz  = 2.0;   // air box height, mm

lc_ring = 0.11;   // element size on the ring
lc_box  = 0.55;   // element size on the box wall

half = N/2;

// --- the two polygons -----------------------------------------------------
For i In {0:N-1}
  theta = 2*Pi*i/N;
  Point(1000 + i) = {Ro*Cos(theta), Ro*Sin(theta), -t/2, lc_ring};
  Point(2000 + i) = {Ri*Cos(theta), Ri*Sin(theta), -t/2, lc_ring};
EndFor

For i In {0:N-2}
  Line(1000 + i) = {1000 + i, 1001 + i};
  Line(2000 + i) = {2000 + i, 2001 + i};
EndFor
Line(1000 + N - 1) = {1000 + N - 1, 1000};
Line(2000 + N - 1) = {2000 + N - 1, 2000};

// The two radial lines the halves share. These are what become the cut
// surfaces once extruded.
Line(3000) = {2000, 1000};                // theta = 0,  x > 0  -- the port
Line(3001) = {2000 + half, 1000 + half};  // theta = pi, x < 0

// --- upper half-annulus, theta 0 -> pi ------------------------------------
// Boundary order fixes the order Extrude returns the lateral surfaces in,
// which is how the cut is identified below without guessing.
upper[] = {3000};
For i In {0:half-1}
  upper[] += {1000 + i};
EndFor
upper[] += {-3001};
For i In {half-1:0:-1}
  upper[] += {-(2000 + i)};
EndFor
Curve Loop(1) = upper[];
Plane Surface(1) = {1};

// --- lower half-annulus, theta pi -> 2pi ----------------------------------
lower[] = {3001};
For i In {half:N-1}
  lower[] += {1000 + i};
EndFor
lower[] += {-3000};
For i In {N-1:half:-1}
  lower[] += {-(2000 + i)};
EndFor
Curve Loop(2) = lower[];
Plane Surface(2) = {2};

// --- extrude both halves together ----------------------------------------
// ONE call, so the surfaces over the two shared radial lines are made once
// and shared -- the same reason cylinder_box.geo extrudes wire and air
// together. Without that the cut would be two coincident faces with a tet
// behind only one of them, and the port would have no plus side.
//
// Per extruded surface the return is [0] = top, [1] = volume, then one entry
// per boundary curve in the order given to Curve Loop. Surface 1 has N + 2
// curves, so surface 2's block starts at N + 4.
ring[] = Extrude {0, 0, t} { Surface{1, 2}; };

vol_upper = ring[1];
vol_lower = ring[N + 5];

// First curve of `upper` is 3000, so the first lateral surface of surface 1
// is the cut at theta = 0.
cut_port  = ring[2];
// First curve of `lower` is 3001, likewise.
cut_other = ring[N + 6];

// --- the air box ----------------------------------------------------------
hx = Wxy/2;
hz = Wz/2;
Point(1) = {-hx, -hx, -hz, lc_box};
Point(2) = { hx, -hx, -hz, lc_box};
Point(3) = { hx,  hx, -hz, lc_box};
Point(4) = {-hx,  hx, -hz, lc_box};
Point(5) = {-hx, -hx,  hz, lc_box};
Point(6) = { hx, -hx,  hz, lc_box};
Point(7) = { hx,  hx,  hz, lc_box};
Point(8) = {-hx,  hx,  hz, lc_box};

Line(1) = {1, 2}; Line(2) = {2, 3}; Line(3) = {3, 4}; Line(4) = {4, 1};
Line(5) = {5, 6}; Line(6) = {6, 7}; Line(7) = {7, 8}; Line(8) = {8, 5};
Line(9) = {1, 5}; Line(10) = {2, 6}; Line(11) = {3, 7}; Line(12) = {4, 8};

Curve Loop(11) = {1, 2, 3, 4};        Plane Surface(11) = {11};   // z = -hz
Curve Loop(12) = {5, 6, 7, 8};        Plane Surface(12) = {12};   // z = +hz
Curve Loop(13) = {1, 10, -5, -9};     Plane Surface(13) = {13};
Curve Loop(14) = {2, 11, -6, -10};    Plane Surface(14) = {14};
Curve Loop(15) = {3, 12, -7, -11};    Plane Surface(15) = {15};
Curve Loop(16) = {4, 9, -8, -12};     Plane Surface(16) = {16};

Surface Loop(21) = {11, 12, 13, 14, 15, 16};

// The ring's outer skin, WITHOUT the two internal cuts: CombinedBoundary
// drops surfaces shared between the volumes it is given, which is exactly
// what a cavity wall must not include.
skin[] = CombinedBoundary { Volume{vol_upper, vol_lower}; };
Surface Loop(22) = skin[];

// Air is the box with the ring as a cavity, sharing the ring's own surfaces,
// so the mesh is conformal across the interface -- the 3D analogue of
// `Plane Surface(2) = {2, 1}` in cylinder_box.geo.
Volume(31) = {21, 22};

// --- physical groups ------------------------------------------------------
Physical Volume("ring", 1) = {vol_upper, vol_lower};
Physical Volume("air", 2)  = {31};

// The port: one complete cross-section of the ring, interior to it, with ring
// tets on both sides.
Physical Surface("loop_cut", 10) = {cut_port};

// The box's own faces need no tag: the outer boundary is found from the
// topology (a face with one adjacent tet), which is what `outer =
// flux_tangential` uses. `cut_other` is deliberately untagged -- it is an
// artefact of building the ring in halves, not a port.

Mesh.MshFileVersion = 4.1;
