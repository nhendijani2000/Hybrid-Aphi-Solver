// ---------------------------------------------------------------------------
// 07_GaugeInvariance -- a short fat copper rod, meshed coarsely on purpose.
//
//     gmsh gauge.geo -3 -o gauge.msh
//
// WHY THIS GEOMETRY. The case asks whether Phi and A move when the spanning
// tree changes while E, B, H and J do not, so Phi has to MOVE BY SOMETHING
// MEASURABLE or the test proves nothing. The gauge-sensitive fraction of Phi
// is about omega*L/|Z| = sin(arg Z) -- see 06_TwoWires_Quadrature's report
// section 7 -- so the case has to sit well away from the resistive limit.
//
// 01_OneCylinder's proportions do that: a = 10 mm, L = 40 mm, W = 200 mm at
// 50 Hz gives omega*L/R = 2.89, where changing the tree moves Phi by 15 % of
// the applied volt. A long thin wire at the same frequency would give
// omega*L/R ~ 0.07 and Phi would barely move, and the test would pass for the
// wrong reason. verify_gauge.py MEASURES the ratio and asserts the regime
// rather than trusting this comment.
//
// WHY IT IS COARSE. Five solves run per check -- three spanning trees under a
// voltage drive and two under a current drive -- so element count is the
// budget. Nothing here is compared against an analytic field: the comparison is
// between two solves of the SAME problem, so discretisation error cancels
// exactly and a coarse mesh costs nothing. The skin effect is 01_OneCylinder's
// job, and the Bessel profile is 04's.
// ---------------------------------------------------------------------------

SetFactory("Built-in");

a = 10.0;    // wire radius, mm -- circumradius of the polygon
W = 200.0;   // box side, mm
L = 40.0;    // height, mm (wire and box alike)
N = 20;      // sides of the wire polygon -- 36 in case 01, coarser here

// Sizing. Deliberately loose: see the header. The polygon facet is
// 2a sin(pi/N) = 3.13 mm at N = 20, and the design rule from case 01 is that
// the volume size stays COMPARABLE to the facet width, so 3.5 mm is the floor
// here rather than a number worth tuning.
lc_skin = 3.5;
lc_core = 3.5;
lc_far  = 30.0;
d_far   = 45.0;

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

// Shared curves, so the mesh is conformal across the wire/air interface.
Plane Surface(2) = {2, 1};

// --- extrude both together --------------------------------------------------
// Return layout per extruded surface: [0] = top, [1] = volume, then one entry
// per lateral face. Surface 1 has N boundary curves, so the air surface's
// entries start at N+2.
ext[] = Extrude {0, 0, L} { Surface{1, 2}; };

wire_volume = ext[1];
wire_top    = ext[0];
air_volume  = ext[N + 3];

lateral[] = {};
For i In {0:N-1}
  lateral[] += ext[2 + i];
EndFor

Physical Volume("wire", 1) = {wire_volume};
Physical Volume("air", 2)  = {air_volume};

// BOTH END FACES MUST BE TAGGED, and this is the point of the case rather than
// a detail. The ports live on these two surfaces, they lie IN the air box's own
// end faces because the wire spans the full height, and the outer boundary
// carries n x A = 0. That is what pins the gauge function psi to zero on the
// whole outer boundary -- so the terminals stay exact in every gauge while the
// interior of Phi does not. verify_gauge.py measures psi on both.
Physical Surface("wire_bottom", 10) = {1};
Physical Surface("wire_top", 11)    = {wire_top};

// --- graded sizing ----------------------------------------------------------
// A Distance field from the wire's lateral faces through a Threshold, with the
// interior capped by Restrict+Min. The distance field is UNSIGNED and reads a
// at the axis, so without the cap the core coarsens as if it were far away --
// see case 01's geo for the full account of why that matters there. It matters
// less here, where no analytic profile is being compared, but the construction
// is kept so the two meshes differ only in size.
Field[1] = Distance;
Field[1].SurfacesList = {lateral[]};
Field[1].Sampling = 60;

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = lc_skin;
Field[2].SizeMax = lc_far;
Field[2].DistMin = 0.0;
Field[2].DistMax = d_far;

Field[3] = MathEval;
Field[3].F = Sprintf("%g", lc_core);   // a bare "lc_core" stores the LITERAL STRING

Field[4] = Restrict;
Field[4].InField = 3;
Field[4].VolumesList = {wire_volume};

Field[5] = Min;
Field[5].FieldsList = {2, 4};

Background Field = 5;
