// ---------------------------------------------------------------------------
// 05_Loop_1A_50Hz -- a conducting RING driven by an internal 1 A current source.
//
//     gmsh loop.geo -3 -o loop.msh
//
// WHY THIS CASE EXISTS
//
// Cases 01-04 are all a straight rod: simply connected, and the tree-cotree
// gauge has nothing hard to do. A closed ring is MULTIPLY CONNECTED, which is
// where the gauge is actually tested -- the cotree has to span the loop, and a
// current circulates with no terminal on any outer boundary to attach an
// electrode to. A ring can only be driven by CUTTING it.
//
// GEOMETRY. Loop radius R = 6.5 mm, wire radius a = 0.8 mm, copper, inside a
// cylinder of air 30 mm in radius and 30 mm tall, the ring at mid height.
//
// Those dimensions are not arbitrary. They come from two PUBLISHED closed forms
// solved simultaneously against a reference impedance measured for this
// geometry:
//
//     R_dc = 2 pi R / (sigma A)
//     L    = mu0 R [ ln(8R/a) - 2 + 1/4 ]      round wire, uniform current
//
// The 1/4 is the internal inductance of a round wire, mu0/(8 pi) per unit
// length, which applies while a << delta -- true here, a/delta = 0.086 at
// 50 Hz. R = 6.5 mm with a = 0.8 mm reproduces the reference R to 1.3 % and L
// to 0.27 %, and agrees independently with the dimensions measured off the
// reference model's own scale bar: centreline diameter 1.23 cm, tube 1.4-2.0 mm.
//
// THE REFERENCE WAS RUN AT 100 Hz, NOT 50. That was not stated anywhere; it
// falls out of the same two equations. Im(Z)/omega is 19.75 nH at 100 Hz and
// 39.50 nH at 50 Hz, and no geometry consistent with the reported R gives
// 39.50 nH -- the closed form gives 19.74 nH, matching the 100 Hz reading to
// 0.05 %. It does not matter for what this case checks: a/delta is 0.086 at
// 50 Hz and 0.121 at 100 Hz, so there is no skin effect either way and R and L
// are frequency independent between them.
//
// ---------------------------------------------------------------------------
// HOW THE CUT IS MADE, which is the whole difficulty.
//
// A port surface needs a tet on each side, so the cut cannot be a surface drawn
// inside a volume -- it has to BE the interface between two volumes. So the
// ring is revolved in TWO HALVES from one cross-section, and the two copies of
// that cross-section, at theta = 0 and theta = pi, become real interior faces
// with tets on both sides. The one at theta = 0 is tagged `loop_cut`.
//
// NOTE IT NEEDS TWO CUTS. A ring cut in one place is still one connected
// volume: the current can go round the other way. Only cutting at both theta=0
// and theta=pi gives two half-rings.
//
// WHY THE CROSS-SECTION IS A POLYGON AND NOT A CIRCLE. Partly the reason every
// other conductor in this project is polygonal -- planar lateral faces mean
// every surface triangle on them has a normal in the plane, so dPhi/dn = 0
// holds exactly. But here it is also the only thing that MESHES. An OCC
// `Torus` carrying an internal cut and fragmented against the air fails with
// "Invalid boundary mesh (overlapping facets)": a curved torus face ends up
// overlapping itself. That was tried eight ways -- partial tori welded with
// Coherence, cuts by oversized discs, by exact-size discs, by thin boxes, cuts
// moved off the parametric seam, the seam rotated away from the cuts, the air
// fragmented before and after the cut -- and every one of them failed the same
// way. Revolving a polygon has none of that trouble.
// ---------------------------------------------------------------------------

SetFactory("OpenCASCADE");

R  = 6.5;     // loop radius, mm -- centreline of the tube
a  = 0.8;     // tube circumradius, mm
M  = 20;      // sides of the tube cross-section
Rd = 60.0;    // air cylinder radius, mm -- 2x, to test the boundary effect on L
Hd = 60.0;    // air cylinder height, mm

// Mesh sizes. lc_ring has to resolve a cross-section 1.6 mm across, so 0.25 mm
// gives about 20 elements around the tube and 6 across it.
//
// AND IT IS TIED TO M BY THE FACET RULE. The facet is 2a sin(pi/M) = 0.250 mm
// at M = 20, matched to lc_ring. A volume size well ABOVE the facet width makes
// slivers -- that is what killed the N = 96 experiment in case 02 -- so M and
// lc_ring move together.
lc_ring = 0.25;
lc_far  = 10.00;
d_far   = 40.0;

// --- the tube cross-section, in the x-z plane centred at (R, 0, 0) ----------
For i In {0 : M-1}
  th = 2*Pi*i/M;
  Point(100+i) = {R + a*Cos(th), 0, a*Sin(th)};
EndFor
For i In {0 : M-2}
  Line(100+i) = {100+i, 101+i};
EndFor
Line(100+M-1) = {100+M-1, 100};
Curve Loop(1) = {100 : 100+M-1};
Plane Surface(1) = {1};

// --- revolve it both ways, so the two cross-sections are shared faces -------
half_pos() = Extrude { {0,0,1}, {0,0,0},  Pi } { Surface{1}; };
half_neg() = Extrude { {0,0,1}, {0,0,0}, -Pi } { Surface{1}; };
// Extrude returns [top surface, volume, lateral faces...]; entry 1 is the volume.
vp = half_pos(1);
vn = half_neg(1);

// --- the air, welded to the ring --------------------------------------------
Cylinder(500) = {0, 0, -Hd/2, 0, 0, Hd, Rd};
BooleanFragments{ Volume{vp, vn, 500}; Delete; }{}

// --- find the pieces again --------------------------------------------------
// The boolean renumbers everything, so the volumes are identified by WHERE THEY
// ARE. `In BoundingBox` returns only entities lying entirely inside the box,
// which separates the halves cleanly: the theta in [0,pi] half reaches y = -a
// at worst, never y = -(R+a), so a box stopping just below -a excludes the
// other half.
eps = 1e-6;
ring_pos() = Volume In BoundingBox{ -(R+a)-eps, -a-eps,     -a-eps,
                                     (R+a)+eps,  (R+a)+eps,  a+eps };
ring_neg() = Volume In BoundingBox{ -(R+a)-eps, -(R+a)-eps, -a-eps,
                                     (R+a)+eps,  a+eps,      a+eps };
all_vol()  = Volume "*";

If (#ring_pos() != 1 || #ring_neg() != 1)
  Error("expected one volume per ring half, got %g and %g",
        #ring_pos(), #ring_neg());
EndIf

air() = {};
For i In {0 : #all_vol()-1}
  v = all_vol(i);
  If (v != ring_pos(0) && v != ring_neg(0))
    air() += v;
  EndIf
EndFor

// The cut: the cross-section at theta = 0, the y = 0 plane on the +x side. The
// other shared cross-section is at theta = pi, on the -x side, so a box around
// +x picks out one and only one.
cut() = Surface In BoundingBox{ R-a-eps, -eps, -a-eps,
                                R+a+eps,  eps,  a+eps };
If (#cut() != 1)
  Error("expected exactly one cut surface at theta = 0, got %g", #cut());
EndIf

Physical Volume("ring", 1) = { ring_pos(0), ring_neg(0) };
Physical Volume("air", 2)  = air();
Physical Surface("loop_cut", 20) = cut();

// The outer boundary needs no tag: it is found from the topology -- a face with
// one adjacent tet -- which is what `outer = flux_tangential` uses.

// --- sizing -----------------------------------------------------------------
// Graded off the ring's own surface. There is no unsigned-distance trap to
// avoid here, unlike the straight-rod cases: the ring is a thin tube, so
// "distance from its surface" is small everywhere inside it and the field
// cannot coarsen an interior that is only 1.6 mm across.
ring_surf() = Boundary{ Volume{ ring_pos(0), ring_neg(0) }; };

Field[1] = Distance;
Field[1].SurfacesList = { ring_surf() };
Field[1].Sampling = 200;

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = lc_ring;
Field[2].SizeMax = lc_far;
Field[2].DistMin = a;        // hold the fine size one tube radius out
Field[2].DistMax = d_far;

Background Field = 2;

Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;

Mesh.Algorithm3D = 1;          // Delaunay, which respects a background field well
Mesh.Optimize = 1;
Mesh.OptimizeNetgen = 1;       // element SHAPE at fixed size
Mesh.MshFileVersion = 4.1;
