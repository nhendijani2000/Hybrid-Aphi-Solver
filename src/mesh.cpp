#include "aphi_solver/mesh.hpp"

#include <algorithm>
#include <cmath>

namespace aphi_solver {

double Vec3::norm() const { return std::sqrt(x * x + y * y + z * z); }

namespace {

// Canonical key for an edge: the vertex pair sorted ascending.
std::pair<int, int> edge_key(int i, int j) {
    return (i < j) ? std::make_pair(i, j) : std::make_pair(j, i);
}

// Canonical key for a face: the vertex triple sorted ascending.
std::array<int, 3> face_key(int a, int b, int c) {
    std::array<int, 3> k{a, b, c};
    std::sort(k.begin(), k.end());
    return k;
}

}  // namespace

void Mesh::build_topology() {
    edges.clear();
    faces.clear();
    tet_edges.assign(tets.size(), {});
    tet_edge_signs.assign(tets.size(), {});
    tet_faces.assign(tets.size(), {});

    // Keep tet_tags parallel to tets. The Gmsh reader fills it as it reads,
    // so sizes already agree there; a mesh built directly in code (every
    // test fixture) has no tags at all, and -1 is what "untagged" means.
    if (tet_tags.size() != tets.size()) {
        tet_tags.assign(tets.size(), -1);
    }

    // Deduplication is done by SORTING, not by a lookup table.
    //
    // The obvious implementation keeps a std::map from key to index and
    // probes it once per (tet, local entity) -- 6T edge probes and 4T face
    // probes, each a red-black tree descent through scattered nodes. That
    // was measured as the single largest cost in the whole read: 23.8 ms of
    // a 53 ms load on a 16 040-tet mesh, more than the ASCII parsing.
    //
    // Instead: emit one (key, slot) record per local entity, sort so equal
    // keys land adjacent, and make one linear scan that numbers the unique
    // keys and fills tet_edges/tet_faces. A sort and a sequential scan are
    // both cache-friendly and branch-predictable, which a tree is not.
    //
    // It also leaves `edges` and `faces` in ascending key order, so
    // find_edge/find_face can binary-search them directly and the two lookup
    // maps disappear -- about 3 MB of tree nodes on this mesh, as well as
    // their time.
    //
    // Measured on meshes/cylinder_box.msh: build_topology 23.8 -> 11.7 ms,
    // the whole read 53.2 -> 38.3 ms, boundary_edge_mask 0.28 -> 0.15 ms
    // (it calls find_edge about 98 000 times). A 2x win, not the 3-5x a
    // first guess suggested: the sort is still O(n log n) and
    // comparison-based, so what this buys is locality, not a better order.
    //
    // The edge and face numbering changes as a result -- it is now sorted
    // rather than discovery order. Nothing depends on it: indices are DOF
    // labels, and the counts the tests pin (tree size, free unknowns) are
    // topological invariants.
    const std::size_t nt = tets.size();

    struct EdgeRef {
        int lo, hi, slot;  // slot encodes tet * 6 + local edge
    };
    std::vector<EdgeRef> eref(nt * 6);
    for (std::size_t t = 0; t < nt; ++t) {
        const TetVerts& tv = tets[t];
        for (int le = 0; le < 6; ++le) {
            const int v0 = tv[static_cast<std::size_t>(kTetLocalEdgeVerts[static_cast<std::size_t>(le)].first)];
            const int v1 = tv[static_cast<std::size_t>(kTetLocalEdgeVerts[static_cast<std::size_t>(le)].second)];
            EdgeRef& r = eref[t * 6 + static_cast<std::size_t>(le)];
            r.lo = v0 < v1 ? v0 : v1;
            r.hi = v0 < v1 ? v1 : v0;
            r.slot = static_cast<int>(t * 6 + static_cast<std::size_t>(le));
            // +1 when the local traversal v0 -> v1 already matches the stored
            // canonical direction (low index -> high index), -1 when it is the
            // reverse -- see Mesh::tet_edge_signs.
            tet_edge_signs[t][static_cast<std::size_t>(le)] =
                (v0 < v1) ? static_cast<signed char>(1) : static_cast<signed char>(-1);
        }
    }
    std::sort(eref.begin(), eref.end(), [](const EdgeRef& a, const EdgeRef& b) {
        return a.lo != b.lo ? a.lo < b.lo : a.hi < b.hi;
    });
    for (std::size_t i = 0; i < eref.size(); ++i) {
        if (i == 0 || eref[i].lo != eref[i - 1].lo || eref[i].hi != eref[i - 1].hi) {
            edges.emplace_back(eref[i].lo, eref[i].hi);
        }
        const std::size_t slot = static_cast<std::size_t>(eref[i].slot);
        tet_edges[slot / 6][slot % 6] = static_cast<int>(edges.size()) - 1;
    }

    struct FaceRef {
        int a, b, c, slot;  // slot encodes tet * 4 + local face
    };
    std::vector<FaceRef> fref(nt * 4);
    for (std::size_t t = 0; t < nt; ++t) {
        const TetVerts& tv = tets[t];
        for (int lf = 0; lf < 4; ++lf) {
            const auto key = face_key(tv[static_cast<std::size_t>(kTetLocalFaceVerts[static_cast<std::size_t>(lf)][0])],
                                      tv[static_cast<std::size_t>(kTetLocalFaceVerts[static_cast<std::size_t>(lf)][1])],
                                      tv[static_cast<std::size_t>(kTetLocalFaceVerts[static_cast<std::size_t>(lf)][2])]);
            FaceRef& r = fref[t * 4 + static_cast<std::size_t>(lf)];
            r.a = key[0];
            r.b = key[1];
            r.c = key[2];
            r.slot = static_cast<int>(t * 4 + static_cast<std::size_t>(lf));
        }
    }
    std::sort(fref.begin(), fref.end(), [](const FaceRef& x, const FaceRef& y) {
        if (x.a != y.a) return x.a < y.a;
        if (x.b != y.b) return x.b < y.b;
        return x.c < y.c;
    });
    for (std::size_t i = 0; i < fref.size(); ++i) {
        if (i == 0 || fref[i].a != fref[i - 1].a || fref[i].b != fref[i - 1].b ||
            fref[i].c != fref[i - 1].c) {
            faces.push_back({fref[i].a, fref[i].b, fref[i].c});
        }
        const std::size_t slot = static_cast<std::size_t>(fref[i].slot);
        tet_faces[slot / 4][slot % 4] = static_cast<int>(faces.size()) - 1;
    }

    // Invert tet_faces into face -> tets. One pass over the same 4-per-tet
    // data, no extra lookups: tet_faces already holds the global face index.
    face_tets.assign(faces.size(), FaceTets{});
    for (std::size_t t = 0; t < tets.size(); ++t) {
        for (int lf = 0; lf < 4; ++lf) {
            FaceTets& ft = face_tets[static_cast<std::size_t>(tet_faces[t][static_cast<std::size_t>(lf)])];
            if (ft.count < 2) ft.tets[static_cast<std::size_t>(ft.count)] = static_cast<int>(t);
            ++ft.count;  // counted even past 2, so a non-manifold mesh is detectable
        }
    }
}

std::string Mesh::physical_name(int dimension, int tag) const {
    for (const PhysicalName& p : physical_names) {
        if (p.dimension == dimension && p.tag == tag) return p.name;
    }
    return std::string();
}

double Mesh::signed_tet_volume(int t) const {
    const TetVerts& tv = tets[static_cast<std::size_t>(t)];
    const Vec3& p0 = nodes[static_cast<std::size_t>(tv[0])];
    const Vec3& p1 = nodes[static_cast<std::size_t>(tv[1])];
    const Vec3& p2 = nodes[static_cast<std::size_t>(tv[2])];
    const Vec3& p3 = nodes[static_cast<std::size_t>(tv[3])];
    const Vec3 a = p1 - p0;
    const Vec3 b = p2 - p0;
    const Vec3 c = p3 - p0;
    return a.cross(b).dot(c) / 6.0;
}

int Mesh::find_edge(int i, int j) const {
    // `edges` is left in ascending key order by build_topology, so a binary
    // search over it replaces what used to be a separate std::map -- faster
    // to probe (contiguous, cache-friendly) and one less structure to hold.
    const auto key = edge_key(i, j);
    const auto it = std::lower_bound(edges.begin(), edges.end(), key);
    return (it == edges.end() || *it != key) ? -1 : static_cast<int>(it - edges.begin());
}

int Mesh::find_face(int a, int b, int c) const {
    const auto key = face_key(a, b, c);
    const auto it = std::lower_bound(faces.begin(), faces.end(), key);
    return (it == faces.end() || *it != key) ? -1 : static_cast<int>(it - faces.begin());
}

}  // namespace aphi_solver
