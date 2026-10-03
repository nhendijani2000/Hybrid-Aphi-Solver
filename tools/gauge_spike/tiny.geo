// A deliberately tiny two-material box for the gauge spike: a copper slab
// inside a dielectric box, so eps varies and K_NE's material weighting is
// actually exercised. A few hundred tets, small enough for dense linear
// algebra so condition numbers can be computed exactly rather than estimated.
SetFactory("OpenCASCADE");
a = 2.0; W = 8.0; L = 8.0;
Box(1) = {-a/2, -a/2, -L/2,  a, a, L};
Box(2) = {-W/2, -W/2, -L/2,  W, W, L};
BooleanFragments{ Volume{1,2}; Delete; }{}
eps = 1e-6;
inner[] = Volume In BoundingBox{-a/2-eps,-a/2-eps,-L/2-eps, a/2+eps,a/2+eps,L/2+eps};
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
Mesh.MeshSizeMin = 2.0;
Mesh.MeshSizeMax = 3.0;
