// One solid cylinder in a box, sized to match the Ansys Maxwell A-Phi voltage
// example -- a slender conductor at 50 Hz, where the potential is readable.
//
//     gmsh -3 cylinder.geo -o cylinder.msh
//
// ---------------------------------------------------------------------------
// WHERE THESE DIMENSIONS COME FROM
//
// They are inferred from the Ansys plot, not measured -- it carries no scale
// bar. Two things fix them.
//
// The plot shows a clean continuous gradient, which means the case is
// resistance dominated. omega*L/R scales as a^2, and we have it measured at
// three points: 0.0358 (the Ansys loop, a = 1.12 mm, clean), 0.0593 (our 1 Hz
// run, clean) and 2.889 (a 10 mm rod at 50 Hz, ruined). Taking 0.1 as a
// generous ceiling for clean gives
//
//     a < 10 mm * sqrt(0.1/2.889) = 1.86 mm
//
// The plotted range is 0 to 3.2836e-04 V over three rods, so about 1.1e-04 V
// each. With a = 1.5 mm and l = 40 mm that is R = 9.76e-05 ohm and hence
// I = 1.13 A -- landing on a 1 A excitation, the conventional setup and the
// one the loop reference also used. That self-consistency is why 1.5 mm.
//
// Length and box size follow from the picture proportions (rods about 1/12 of
// the domain diameter, spanning most of its height) and are the weakest
// numbers here: proportions off by 2x move them by 2x.
//
// a/delta = 0.16 at 50 Hz, so there is NO skin effect to resolve. The surface
// refinement is for the field gradient outside the conductor, not a boundary
// layer. Contrast 01_OneCylinder, which is the opposite case.
//
// N = 24 rather than 36: at a = 1.5 mm a 36-gon has 0.26 mm facets, and every
// facet edge must carry an element whatever lc_skin says, which drove the mesh
// to 18696 nodes. 24 facets give 0.39 mm edges and a 1.1 % area deficit
// against pi a^2 -- still the polygon, not the circle, that R is compared to.
//
// N = 96 WAS TRIED AND REVERTED (2026-09-28). The 24-gon's 15-degree facets were
// a candidate explanation for the azimuthal lobing in |B| at the conductor
// surface. N = 96 makes the cross-section round to 0.07 % in area and costs only
// 59 % more nodes -- but it does NOT smooth B, and it introduces a worse
// artefact. Measured, N=24 -> N=96, at identical volume sizing:
//
//     azimuthal sd/mean   0.029-0.057  ->  0.027-0.063     unchanged
//     unknowns            114390       ->  183150
//     factorise           290 s, 2.1 GB -> 795 s, 4.2 GB
//
// The reason is element QUALITY, not count. The 96 boundary nodes are forced by
// the polygon, but the interior only supports lc_core = 0.5 mm, so gmsh fans
// dozens of slivers from each interior node out to the dense boundary. Those
// slivers put feathery radial spikes into |B| at r = a, replacing the lobing.
//
// THE RULE: facet(N) = 2a sin(pi/N) must be comparable to the VOLUME size lc.
//
//     lc = 0.5 mm -> matched N ~ 19        N = 24 -> facet 0.39 mm, matched
//     lc = 0.7 mm -> matched N ~ 13        N = 96 -> facet 0.098 mm, 5x too fine
//
// Using N = 96 honestly needs lc_skin ~ 0.1 mm, about 125x the elements in the
// surface shell. The pairing here is lc_skin = 0.25 with N = 40 (facet 0.235
// against a 0.25 mm volume size); at lc_skin = 0.35 it was N = 32, at 0.7 N = 24.
// ---------------------------------------------------------------------------
// frequency and a sweep upward walks into strong skin effect.
//
//     gmsh -3 cylinder.geo -o cylinder.msh
//
SetFactory("Built-in");

a = 1.5;     // wire radius, mm -- circumradius of the polygon
W = 40.0;    // box side, mm
L = 40.0;    // height, mm (wire and box alike)
N = 40;      // sides of the wire polygon -- matched to lc_skin, see the note above

// Mesh sizing. These are the field parameters, NOT point sizes: the field
// below overrides point sizes entirely (see CharacteristicLengthExtendFromBoundary).
lc_skin = 0.25;   // at the conductor surface -- the lever that smooths B
lc_far  = 6.0;    // out in the air, where nothing happens
d_far   = 4.0;    // distance over which one grows into the other

// The conductor core is NOT covered by the distance field above: distance from
// the lateral face reaches only a = 1.5 mm on the axis, so the core is sized at
// lc_skin + (1.5/d_far)(lc_far - lc_skin) = 1.2 mm against a 1.5 mm radius --
// about ONE element spanning the core. MEASURED consequence: only 178 tet
// centroids land in the annulus 0.60-0.90 mm over the middle 40 % of the
// length, and |B| there came out 11.5 % high against the exact mu0*I*r/2*pi*a^2
// while every well-resolved band sat within 3 %. lc_core caps the size inside
// the wire and fixes it. At a/delta = 0.16 there is no boundary layer, so a
// UNIFORM wire mesh is the right target; lc_core < lc_skin simply makes it so.
lc_core = 0.4;    // inside the conductor, everywhere

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
// Distance from the lateral surface reaches only a = 1.5 mm at the axis, so the
// core would be meshed at roughly lc_skin + (1.5/d_far)(lc_far - lc_skin)
// = 1.2 mm -- far too coarse. Field 4 below overrides it.
Field[1] = Distance;
Field[1].SurfacesList = {lateral[]};
Field[1].Sampling = 100;

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = lc_skin;
Field[2].SizeMax = lc_far;
Field[2].DistMin = 0.0;
Field[2].DistMax = d_far;

// The threshold field alone leaves the core coarse (see lc_core above), so cap
// it inside the wire and take the smaller of the two everywhere. Restrict
// returns a huge size outside its volume, which is exactly what Min wants.
Field[3] = MathEval;
Field[3].F = Sprintf("%g", lc_core);

Field[4] = Restrict;
Field[4].InField = 3;
Field[4].VolumesList = {wire_volume};

Field[5] = Min;
Field[5].FieldsList = {2, 4};

Background Field = 5;

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
