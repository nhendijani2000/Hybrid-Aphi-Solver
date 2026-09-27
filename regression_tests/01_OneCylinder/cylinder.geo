// One solid cylinder in a box, sized so that 50 Hz is already the interesting
// frequency and a sweep upward walks into strong skin effect.
//
//     gmsh -3 cylinder.geo -o cylinder.msh
//
// ---------------------------------------------------------------------------
// WHY THIS SIZE
//
// Skin depth in copper is delta = sqrt(2/(w mu sigma)) = 9.35 mm at 50 Hz. The
// previous geometry (a = 0.2 mm) had a/delta = 0.02, so there was no skin
// effect to see at any frequency a low-frequency solver cares about. Here
// a = 10 mm gives a/delta = 1.07 at 50 Hz: the current is already starting to
// crowd, and every decade up makes it stronger.
//
//     f        delta     a/delta   elem per delta   |J(0)/J(a)| exact
//     50 Hz    9.346 mm    1.07        4.67             0.9265
//     200 Hz   4.673 mm    2.14        2.34             0.5055
//     500 Hz   2.955 mm    3.38        1.48             0.1826   <- marginal
//     1 kHz    2.090 mm    4.79        1.04             0.0538   <- not resolved
// So this mesh is honest at 50 Hz (measured: 0.948 against the exact 0.9414 for
// the same bin centres, 0.7 %) and marginal by 200 Hz. Going higher needs lc_skin
// reduced, which costs elements fast -- and the direct solver is the binding
// constraint, not gmsh: 49276 unknowns already take 93 s to factor and 751 MB.
//
// The box is 200 mm across, W/a = 20, and the wire spans its full height so
// both caps lie on the boundary and the problem stays a clean two-terminal
// series drive.
//
// ---------------------------------------------------------------------------
// WHY THE WIRE IS A POLYGON, NOT A CIRCLE
//
// The DC resistance should come out exact, and that rests on Phi = V*z/L being
// reproducible by the element space. It is linear, so it lies in P1 subset P2 --
// but it must also satisfy dPhi/dn = 0 on the lateral surface, which needs that
// surface vertical. Mesh a true cylinder and its boundary triangles have
// vertices at different heights and angles, giving a normal with a z-component:
// grad(Phi).n is then not zero and the linear field is no longer the discrete
// solution. A regular N-gon makes every lateral face planar and vertical by
// construction.
//
// The price is that the meshed cross-section is the polygon's area, not pi*a^2.
// With N = 36 it is 0.51 % below -- and R must be compared against the polygon.
//
//     A_poly = (N/2) a^2 sin(2 pi / N) = 311.87 mm^2   (pi a^2 = 314.16)
// ---------------------------------------------------------------------------

SetFactory("Built-in");

a = 10.0;    // wire radius, mm -- circumradius of the polygon
W = 200.0;   // box side, mm -- W/a = 20, so the outer wall is genuinely far
L = 40.0;    // height, mm (wire and box alike)
N = 36;      // sides of the wire polygon

// Mesh sizing. These are the field parameters, NOT point sizes: the field
// below overrides point sizes entirely (see CharacteristicLengthExtendFromBoundary).
lc_skin = 2.0;    // at the conductor surface, where the skin depth lives
lc_far  = 16.0;   // out in the air, where nothing happens
d_far   = 45.0;   // distance over which one grows into the other

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

// The inner loop is the wire's own boundary, so the two surfaces share those
// curves and the mesh is conformal across the wire/air interface.
Plane Surface(2) = {2, 1};

// --- extrude both together --------------------------------------------------
// Both in ONE call, so the wire's lateral surface is created once and shared.
// No Layers{}: a structured extrusion of triangles gives prisms, and this
// project's reader takes tetrahedra only.
//
// Return layout per extruded surface: [0] = top, [1] = volume, then one entry
// per lateral face. Surface 1 has N boundary curves, so the air surface's
// entries start at N+2.
ext[] = Extrude {0, 0, L} { Surface{1, 2}; };

wire_volume = ext[1];
wire_top    = ext[0];
air_volume  = ext[N + 3];

// The wire's N lateral faces, which is what the refinement is measured from.
lateral[] = {};
For i In {0:N-1}
  lateral[] += ext[2 + i];
EndFor

Physical Volume("wire", 1) = {wire_volume};
Physical Volume("air", 2)  = {air_volume};
Physical Surface("wire_bottom", 10) = {1};
Physical Surface("wire_top", 11)    = {wire_top};

// The rest of the box needs no tag: the outer boundary is found from the
// topology (a face with one adjacent tet), which is what `outer =
// flux_tangential` uses.

// --- graded sizing: fine at the conductor surface, coarse far away ----------
// A Distance field measured from the wire's lateral faces, mapped through a
// Threshold. This refines on BOTH sides of the surface, which is what is
// wanted: the skin effect is inside, and the 1/r field outside wants some
// resolution too, but neither needs it more than a few skin depths away.
//
// Distance from the lateral surface reaches a = 10 mm at the axis, so the core
// of the conductor is meshed at roughly lc_skin + (10/d_far)(lc_far - lc_skin)
// = 3.8 mm -- coarse, and correctly so, since nothing happens there.
Field[1] = Distance;
Field[1].SurfacesList = {lateral[]};
Field[1].Sampling = 100;

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = lc_skin;
Field[2].SizeMax = lc_far;
Field[2].DistMin = 0.0;
Field[2].DistMax = d_far;

Background Field = 2;

// Without this the point sizes above compete with the field and win near the
// geometry, which would defeat the grading entirely.
Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;

Mesh.Algorithm3D = 1;          // Delaunay, which respects a background field well
Mesh.Optimize = 1;
Mesh.OptimizeNetgen = 0;
Mesh.MshFileVersion = 4.1;

// ---------------------------------------------------------------------------
// TO GO HIGHER IN FREQUENCY
//
// Halving lc_skin buys one factor of 2 in resolved delta, i.e. 4x in frequency,
// and costs roughly 8x the elements in the refined shell. Before reaching for
// it, check whether the answer being sought needs the field resolved at all --
// R and L converge considerably faster than the current profile does.
// ---------------------------------------------------------------------------
