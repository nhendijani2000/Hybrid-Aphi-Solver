// A copper cube fully embedded in a dielectric box, for compare.py.
//
// WHY NOT tiny.geo. tiny.geo's conductor slab runs the full length of the
// domain, so its conductor/dielectric interface meets the outer boundary and
// the mesh is coarse enough that EVERY interior node ends up on that interface.
// That makes the one measurement compare.py exists for -- psi strictly inside
// the insulator against psi only on the interface, Jochum section 5.5 --
// impossible to take: the two sets are identical and the comparison is vacuous.
//
// Here the conductor is fully enclosed, so the interface is closed and there are
// three distinct node populations: strictly inside the conductor, on the
// interface, and strictly inside the insulator. Still small enough that
// condition numbers come from a dense SVD and are exact.
SetFactory("OpenCASCADE");
a = 3.0;      // conductor cube
W = 10.0;     // dielectric box
Box(1) = {-a/2, -a/2, -a/2,  a, a, a};
Box(2) = {-W/2, -W/2, -W/2,  W, W, W};
BooleanFragments{ Volume{1,2}; Delete; }{}

eps = 1e-6;
inner[] = Volume In BoundingBox{-a/2-eps,-a/2-eps,-a/2-eps,
                                 a/2+eps, a/2+eps, a/2+eps};
all[] = Volume "*";
outer[] = {};
For i In {0 : #all[]-1}
  isin = 0;
  For j In {0 : #inner[]-1}
    If (all[i] == inner[j]); isin = 1; EndIf
  EndFor
  If (isin == 0); outer[] += all[i]; EndIf
EndFor
Physical Volume("conductor", 1) = {inner[]};
Physical Volume("dielectric", 2) = {outer[]};

Mesh.MeshSizeMin = 1.3;
Mesh.MeshSizeMax = 2.2;
