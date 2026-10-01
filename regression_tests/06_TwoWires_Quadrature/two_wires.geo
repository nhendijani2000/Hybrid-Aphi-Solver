// ---------------------------------------------------------------------------
// 06_TwoWires_Quadrature -- two parallel copper wires driven 90 degrees apart.
//
//     gmsh two_wires.geo -3 -o two_wires.msh
//
// WHY TWO WIRES, AND WHY THE QUADRATURE. This case exists to demonstrate the
// [postprocess] options, and the three that are hardest to understand --
// axial_ratio, peak against complex_magnitude, and phase -- are all about
// ELLIPTICAL POLARIZATION. Every other case in the suite is linearly polarised
// almost everywhere: case 05's median axial ratio is 9e-06. A single wire gives
// E along one axis and B purely azimuthal, the polarization ellipse collapses to
// a line, and those options have nothing to show.
//
// Two sources a quarter cycle apart do have something to show. Working the
// superposition of two infinite wires analytically, with I1 = I and I2 = jI:
//
//     b/a = 0      exactly, everywhere on the line of centres (y = 0)
//     b/a = 1.000  at (0, +/- d/2) -- PERFECTLY CIRCULAR, where
//                  complex_magnitude overstates the true peak by the full
//                  sqrt(2) = 41 %
//     b/a > 0.3 over 22 % of the plane
//
// So the axial_ratio figure has real structure: two bright spots on the
// perpendicular bisector, a null along the axis joining the wires.
//
// GEOMETRY. a = 2 mm, centres d = 10 mm apart, 20 mm long, in a 30 x 30 x 20 mm
// box of air. At 2.5 kHz in copper delta = 1.32 mm, so a/delta = 1.5 -- mild
// skin effect, enough that the phase varies with radius inside each conductor
// and `display = phase` has structure there rather than being flat.
//
// THE MESH IS DELIBERATELY COARSE. This is a demonstration you are meant to
// edit and re-run, so it is sized for a ten-second solve rather than for
// accuracy: about 2 elements per skin depth, where case 04 uses 8. The physics
// it asserts (the port currents, and B between the wires against the two-wire
// superposition) is insensitive to that; do not read a skin-depth profile off
// it.
// ---------------------------------------------------------------------------

SetFactory("OpenCASCADE");

a  = 2.0;      // wire radius, mm
d  = 10.0;     // centre-to-centre separation, mm
Lz = 20.0;     // length, mm
W  = 30.0;     // air box half-width is W/2, mm
H  = 30.0;     // air box half-height is H/2, mm

lc_wire = 0.6;
lc_far  = 3.0;
d_far   = 10.0;

// --- the two wires, running the full length so their ends are on the boundary
Cylinder(1) = {-d/2, 0, 0,  0, 0, Lz,  a};
Cylinder(2) = { d/2, 0, 0,  0, 0, Lz,  a};

// --- the air box, same extent in z so the wire ends sit in its end faces -----
Box(3) = {-W/2, -H/2, 0,  W, H, Lz};

BooleanFragments{ Volume{1, 2, 3}; Delete; }{}

// --- find the pieces by where they are --------------------------------------
// After fragmenting, the volume tags are reassigned; locate each by a bounding
// box around where it must be rather than by guessing a number.
eps = 1e-6;
w1() = Volume In BoundingBox{-d/2-a-eps, -a-eps, -eps,  -d/2+a+eps, a+eps, Lz+eps};
w2() = Volume In BoundingBox{ d/2-a-eps, -a-eps, -eps,   d/2+a+eps, a+eps, Lz+eps};
all() = Volume "*";

air() = {};
For i In {0 : #all()-1}
  isw = 0;
  If (all(i) == w1(0)) isw = 1; EndIf
  If (all(i) == w2(0)) isw = 1; EndIf
  If (isw == 0) air() += all(i); EndIf
EndFor

Physical Volume("wire1", 1) = { w1() };
Physical Volume("wire2", 2) = { w2() };
Physical Volume("air",   3) = { air() };

// --- the four terminals, one pair per wire ----------------------------------
// Each wire is independently driven, so each needs its own potential reference:
// a current source at z = 0 and a 0 V terminal at z = Lz.
b1() = Surface In BoundingBox{-d/2-a-eps, -a-eps, -eps,  -d/2+a+eps, a+eps, eps};
t1() = Surface In BoundingBox{-d/2-a-eps, -a-eps, Lz-eps, -d/2+a+eps, a+eps, Lz+eps};
b2() = Surface In BoundingBox{ d/2-a-eps, -a-eps, -eps,   d/2+a+eps, a+eps, eps};
t2() = Surface In BoundingBox{ d/2-a-eps, -a-eps, Lz-eps,  d/2+a+eps, a+eps, Lz+eps};

Physical Surface("wire1_bottom", 11) = { b1() };
Physical Surface("wire1_top",    12) = { t1() };
Physical Surface("wire2_bottom", 21) = { b2() };
Physical Surface("wire2_top",    22) = { t2() };

// --- sizing: fine on the wires, coarsening outward --------------------------
Field[1] = Distance;
Field[1].SurfacesList = { CombinedBoundary{ Volume{ w1(), w2() }; } };
Field[1].Sampling = 100;

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = lc_wire;
Field[2].SizeMax = lc_far;
Field[2].DistMin = a;
Field[2].DistMax = d_far;

Background Field = 2;
Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;

Mesh.Algorithm3D = 1;
Mesh.Optimize = 1;
Mesh.OptimizeNetgen = 1;
Mesh.MshFileVersion = 4.1;
