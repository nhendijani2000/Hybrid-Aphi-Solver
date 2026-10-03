// A copper stub standing on the domain floor, for spikeC.py.
//
// WHY NOT embedded.geo. That one has the conductor fully enclosed, which is
// right for measuring OPERATORS (compare.py) and wrong for measuring a PORT:
// current injected at a terminal on an enclosed conductor can only leave
// through sigma = 0 dielectric, so at 50 Hz there is no return path, the
// problem is ill-posed, and the answer is garbage for reasons that have nothing
// to do with the gauge. That is exactly the error that killed the first attempt
// at 08_MixedPort_Interior.
//
// Here the copper stands ON the box floor, so:
//
//     z = 0    copper meets the outer boundary  -> Phi = 0 reference
//     z = h    copper's top face, INTERIOR      -> the mixed terminal,
//                                                  conductor below, air above
//
// which is 08_MixedPort_Interior's configuration in miniature, with a closed
// conduction path. Deliberately boxy and coarse: a few hundred tets so the
// dense linear algebra stays exact.
SetFactory("OpenCASCADE");
a = 3.0;      // copper cross-section, mm
h = 6.0;      // copper height -- stops short of the box top
W = 10.0;     // box side, mm

Box(1) = {-a/2, -a/2, 0,     a, a, h};
Box(2) = {-W/2, -W/2, 0,     W, W, W};
BooleanFragments{ Volume{1,2}; Delete; }{}

eps = 1e-6;
inner[] = Volume In BoundingBox{-a/2-eps,-a/2-eps,-eps, a/2+eps,a/2+eps,h+eps};
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
