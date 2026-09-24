// Cylinder in a square box -- the first end-to-end test geometry.
// See `Claude outputs/input_file_plan.md` Sec. 5 and examples/cylinder_box.aphi.
//
// Meshed with **Gmsh 4.13.1** -- the same version that produced
// meshes/validation_box_v*.msh and meshes/coax_via.msh. A different version
// can renumber entities and size elements differently, which would silently
// change every count the tests assert. Regenerate with:
//
//     gmsh -3 tools/cylinder_box.geo -o meshes/cylinder_box.msh
//
// ---------------------------------------------------------------------------
// Why the wire is a POLYGON, not a circle.
//
// The DC resistance is meant to come out exact to round-off, and that rests
// on the exact solution Phi = V*z/L being reproducible by the element space.
// It is linear, so it lies in P1 subset P2 -- but it also has to satisfy the
// natural condition dPhi/dn = 0 on the conductor's lateral surface, which
// needs that surface to be vertical.
//
// Mesh a true cylinder and its boundary triangles generally have three
// vertices at three different heights and angles, giving a normal with a
// z-component: grad(Phi).n is then not zero and the linear field is no
// longer the discrete solution. Make the cross-section a regular N-gon and
// the lateral faces are planar and vertical *by construction*, so every
// surface triangle lying in them has a horizontal normal no matter how the
// volume is triangulated.
//
// The price is that the meshed cross-section is the polygon's area, not
// pi*a^2 -- which is why the test asserts R against the measured A_mesh.
// ---------------------------------------------------------------------------

SetFactory("Built-in");

a  = 0.2;    // wire radius, mm (circumradius of the polygon)
W  = 2.0;    // box side, mm
L  = 1.0;    // height, mm
N  = 24;     // sides of the wire polygon

lc_wire = 0.055;   // element size on the wire
lc_box  = 0.30;    // element size on the box wall

// --- wire cross-section: a regular N-gon inscribed in radius a ------------
For i In {0:N-1}
  theta = 2*Pi*i/N;
  Point(100+i) = {a*Cos(theta), a*Sin(theta), 0, lc_wire};
EndFor
For i In {0:N-2}
  Line(100+i) = {100+i, 101+i};
EndFor
Line(100+N-1) = {100+N-1, 100};

Curve Loop(1) = {100:100+N-1};
Plane Surface(1) = {1};

// --- box cross-section: a square with the wire removed --------------------
h = W/2;
Point(1) = {-h, -h, 0, lc_box};
Point(2) = { h, -h, 0, lc_box};
Point(3) = { h,  h, 0, lc_box};
Point(4) = {-h,  h, 0, lc_box};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Curve Loop(2) = {1, 2, 3, 4};

// The inner loop is the wire's own boundary, so the two surfaces share those
// curves and the mesh is conformal across the wire/air interface.
Plane Surface(2) = {2, 1};

// --- extrude both together ------------------------------------------------
// Both in ONE call, so the wire's lateral surface is created once and shared.
// No Layers{}: a structured extrusion of triangles gives prisms, which this
// project's reader does not take (it reads tetrahedra only). Extruding the
// geometry and letting the 3D algorithm fill it gives tets, and the vertical
// lateral faces -- the property the exactness argument needs -- come from the
// polygon, not from the mesh structure.
//
// Return layout, per extruded surface: [0] = top, [1] = volume, then one
// entry per lateral face. Surface 1 has N boundary curves, so the air
// surface's entries start at N+2.
ext[] = Extrude {0, 0, L} { Surface{1, 2}; };

wire_volume = ext[1];
wire_top    = ext[0];
air_volume  = ext[N + 3];

Physical Volume("wire", 1) = {wire_volume};
Physical Volume("air", 2)  = {air_volume};

// The wire's two end caps. Both lie on the box boundary (the wire spans the
// full height), so they are the two boundary ports of examples/cylinder_box.aphi.
Physical Surface("wire_bottom", 10) = {1};
Physical Surface("wire_top", 11)    = {wire_top};

// The rest of the box needs no tag: the outer boundary is found from the
// topology (a face with one adjacent tet), which is what `outer =
// flux_tangential` uses.

Mesh.MshFileVersion = 4.1;
