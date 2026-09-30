// ---------------------------------------------------------------------------
// 04_Cylinder_SkinDepth -- the same copper rod as 01_OneCylinder, run where the
// SKIN EFFECT DOMINATES rather than where it is a 2.6 % correction.
//
//     gmsh cylinder_skin.geo -3 -o cylinder_skin.msh
//
// WHY THIS CASE EXISTS
//
// 01_OneCylinder runs a = 10 mm at 50 Hz, so a/delta = 1.07 and R_ac/R_dc =
// 1.027. The skin effect is real but it is only 2.6 % of the answer, so
// matching Kelvin to 0.04 % validates the EFFECT ITSELF to only about 1.4 %.
// This case runs the same rod at 393 Hz:
//
//                       01 @ 50 Hz      04 @ 393 Hz
//     a/delta              1.070           3.000
//     |J(0)/J(a)| exact    0.9265          0.2522      <- 4x surface-to-axis
//     R_ac/R_dc exact      1.0267          1.7680      <- 77 % of R, not 2.6 %
//
// At a/delta = 3 the skin effect is 77 % of R, so the same 0.1 % agreement
// validates the effect to 0.13 %. That is the point: not a higher frequency for
// its own sake, but a regime where the quantity being checked is most of the
// answer instead of a correction to it.
//
// WHY 393 Hz AND NOT MORE. At a/delta = 5 (1092 Hz) the exact |J(0)/J(a)| is
// 0.054, so the core current is 5 % of the surface value and a relative error
// there is being taken against nearly zero -- the reference stops being
// well conditioned. It also needs N = 128 and about 1.7 M unknowns. a/delta = 3
// is the most demanding point that still has a meaningful exact reference.
//
// ---------------------------------------------------------------------------
// HOW THE MESH WAS SIZED, which is the whole content of this file
//
// Measured on 01's mesh (N = 36, h = 1.65 mm at the rim), worst-case error in
// the radial |J| profile against Bessel evaluated at the same radii:
//
//     f        delta      a/delta   delta/h   worst band err   R vs Kelvin
//     50 Hz    9.346 mm    1.070     5.66        0.26 %         -0.038 %
//    175 Hz    4.996 mm    2.002     3.03        2.10 %         -0.050 %
//    300 Hz    3.815 mm    2.621     2.31        4.73 %         +0.386 %
//
// So the error is set by delta/h, and delta/h >= 4 is what holds it near 1 %.
// At 393 Hz, delta = 3.334 mm, so h must be about 0.833 mm.
//
// WHERE the refinement has to go was not obvious and was measured, not assumed.
// At 300 Hz the per-band error looks small in the annulus (0.1 % at r > 0.8a)
// and large in the core (4.7 % at r < 0.3a), which reads as "refine the core".
// That is the wrong reading. |J| in the core is FLAT there -- 0.3934, 0.3929,
// 0.3949 across r/a = 0 to 0.3 -- so the core is not failing to resolve local
// variation. It is inheriting accumulated error from the region where the field
// actually decays, r in [a - 2 delta, a], and delta/h was 1.45 to 2.31 across
// the WHOLE conductor at that frequency. The annulus's own relative error looks
// small only because |J| is large there.
//
// At a/delta = 3 the decay region r > a - 2 delta = 3.33 mm is 89 % of the
// cross-section, so there is nothing to be gained by grading inside the
// conductor: lc_core = lc_skin and the wire is meshed uniformly. Keeping the
// deep core coarse would save about 10 % of the elements and reintroduce the
// error this case exists to measure.
//
// AND N HAD TO RISE WITH IT. The facet is 2a sin(pi/N), and the design rule is
// that the volume element size stays COMPARABLE to the facet width; a volume
// size well below the facet makes slivers. At N = 36 the facet is 1.743 mm, so
// h = 0.833 mm would be less than half of it. N = 76 gives a facet of
// 0.827 mm, matched to h. This is the same rule that killed the N = 96
// experiment in case 02 from the other side, where the facet was 0.098 mm and
// the volume size 0.25 mm -- the ratio is what matters, in either direction.
// ---------------------------------------------------------------------------

SetFactory("Built-in");

a = 10.0;    // wire radius, mm -- circumradius of the polygon
W = 200.0;   // box side, mm -- W/a = 20, so the outer wall is genuinely far
L = 40.0;    // height, mm (wire and box alike)
N = 76;      // sides of the wire polygon; facet = 2a sin(pi/N) = 0.827 mm

// The target frequency and its skin depth, kept here because every size below
// is derived from them. delta = sqrt(2/(omega mu sigma)) for copper.
f_target = 393.0;   // Hz
delta    = 3.334;   // mm at 393 Hz in copper, sigma = 5.8e7

lc_skin = 0.833;    // = delta/4, at and just outside the conductor surface
lc_core = 0.833;    // = lc_skin: the decay region is 89 % of the cross-section
lc_far  = 16.0;     // out in the air, where nothing happens
d_hold  = 1.667;    // = delta/2, distance the fine size is HELD outward
d_far   = 45.0;     // distance over which it then grows to lc_far

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

// --- sizing -----------------------------------------------------------------
// TWO fields, because the conductor and the air want different things and the
// Distance field cannot tell them apart -- it is UNSIGNED, so on the axis it
// reads 10 mm, the same as a point 10 mm out in the air. That is the bug that
// left 01_OneCylinder's core at 4.4 mm elements.
//
// Field[2] grades the AIR: it HOLDS lc_skin for the first d_hold = delta/2
// outward, then grows to lc_far. DistMin is the important part. With
// DistMin = 0 -- which is what 01 had -- the fine size decays immediately and
// there is no resolved shell at all, only a fine surface.
//
// Field[4] caps the CONDUCTOR at lc_core uniformly. Restrict returns a huge
// size outside its volume, which is exactly what Min wants, so it binds only
// inside the wire and leaves the air grading untouched.
Field[1] = Distance;
Field[1].SurfacesList = {lateral[]};
Field[1].Sampling = 200;          // 100 was enough at N=36; N=76 has finer facets

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = lc_skin;
Field[2].SizeMax = lc_far;
Field[2].DistMin = d_hold;        // HOLD the fine size this far out, then grow
Field[2].DistMax = d_far;

Field[3] = MathEval;
Field[3].F = Sprintf("%g", lc_core);   // a bare "lc_core" stores the LITERAL STRING

Field[4] = Restrict;
Field[4].InField = 3;
Field[4].VolumesList = {wire_volume};

Field[5] = Min;
Field[5].FieldsList = {2, 4};

Background Field = 5;

// Without this the point sizes above compete with the field and win near the
// wire, giving a mesh that ignores the grading.
Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;

Mesh.Algorithm3D = 1;          // Delaunay, which respects a background field well
Mesh.Optimize = 1;
Mesh.OptimizeNetgen = 1;       // element SHAPE at fixed size; 24 % less B scatter
                               // and a 26 % faster factorization on case 02
Mesh.MshFileVersion = 4.1;
