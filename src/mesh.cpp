#include "aphi_solver/mesh.hpp"

#include <algorithm>
#include <cmath>
#include <map>

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
    edge_lookup_.clear();

    std::map<std::pair<int, int>, int>& edge_index = edge_lookup_;
    std::map<std::array<int, 3>, int> face_index;

    for (std::size_t t = 0; t < tets.size(); ++t) {
        const TetVerts& tv = tets[t];

        // Edges (6 per tet), in kTetLocalEdgeVerts order (Jin 2014 Fig. 5.3).
        for (int le = 0; le < 6; ++le) {
            const int v0 = tv[kTetLocalEdgeVerts[le].first];
            const int v1 = tv[kTetLocalEdgeVerts[le].second];
            const auto key = edge_key(v0, v1);
            auto it = edge_index.find(key);
            int gidx;
            if (it == edge_index.end()) {
                gidx = static_cast<int>(edges.size());
                edges.push_back(key);
                edge_index.emplace(key, gidx);
            } else {
                gidx = it->second;
            }
            tet_edges[t][le] = gidx;
            // +1 when the local traversal v0 -> v1 already matches the stored
            // canonical direction (low index -> high index), -1 when it is the
            // reverse -- see Mesh::tet_edge_signs.
            tet_edge_signs[t][le] = (v0 < v1) ? static_cast<signed char>(1) : static_cast<signed char>(-1);
        }

        // Faces (4 per tet), in kTetLocalFaceVerts order (opposite-vertex).
        for (int lf = 0; lf < 4; ++lf) {
            const int v0 = tv[kTetLocalFaceVerts[lf][0]];
            const int v1 = tv[kTetLocalFaceVerts[lf][1]];
            const int v2 = tv[kTetLocalFaceVerts[lf][2]];
            const auto key = face_key(v0, v1, v2);
            auto it = face_index.find(key);
            int gidx;
            if (it == face_index.end()) {
                gidx = static_cast<int>(faces.size());
                faces.push_back(key);
                face_index.emplace(key, gidx);
            } else {
                gidx = it->second;
            }
            tet_faces[t][lf] = gidx;
        }
    }
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
    const auto key = edge_key(i, j);
    auto it = edge_lookup_.find(key);
    return (it == edge_lookup_.end()) ? -1 : it->second;
}

}  // namespace aphi_solver
