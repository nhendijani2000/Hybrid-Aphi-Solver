#include "aphi_solver/gmsh_reader.hpp"

#include <fstream>
#include <sstream>
#include <unordered_map>

namespace aphi_solver {

namespace {

constexpr int kGmshTetType = 4;  // Gmsh elm-type 4 == 4-node tetrahedron.

std::string trim(const std::string& s) {
    const auto a = s.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) return "";
    const auto b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b - a + 1);
}

}  // namespace

Mesh read_gmsh_msh(const std::string& path) {
    std::ifstream in(path);
    if (!in) {
        throw GmshReadError("read_gmsh_msh: could not open file: " + path);
    }

    Mesh mesh;
    std::unordered_map<long long, int> gmsh_id_to_local;  // Gmsh node id -> Mesh::nodes index

    std::string line;
    bool saw_mesh_format = false;
    bool saw_nodes = false;
    bool saw_elements = false;

    while (std::getline(in, line)) {
        const std::string tag = trim(line);

        if (tag == "$MeshFormat") {
            if (!std::getline(in, line)) {
                throw GmshReadError("read_gmsh_msh: truncated $MeshFormat section");
            }
            std::istringstream fmt(line);
            double version = 0.0;
            int file_type = 0, data_size = 0;
            fmt >> version >> file_type >> data_size;
            if (version < 2.0 || version >= 3.0) {
                throw GmshReadError("read_gmsh_msh: only Gmsh format 2.x is supported, got " +
                                     std::to_string(version));
            }
            if (file_type != 0) {
                throw GmshReadError("read_gmsh_msh: only ASCII .msh files are supported (file_type 0)");
            }
            saw_mesh_format = true;
            continue;
        }

        if (tag == "$Nodes") {
            if (!std::getline(in, line)) throw GmshReadError("read_gmsh_msh: truncated $Nodes section");
            long long num_nodes = 0;
            {
                std::istringstream count_line(line);
                count_line >> num_nodes;
            }
            mesh.nodes.reserve(static_cast<std::size_t>(num_nodes));
            for (long long k = 0; k < num_nodes; ++k) {
                if (!std::getline(in, line)) throw GmshReadError("read_gmsh_msh: $Nodes section too short");
                std::istringstream node_line(line);
                long long gmsh_id;
                double x, y, z;
                node_line >> gmsh_id >> x >> y >> z;
                const int local_idx = static_cast<int>(mesh.nodes.size());
                mesh.nodes.emplace_back(x, y, z);
                gmsh_id_to_local[gmsh_id] = local_idx;
            }
            if (!std::getline(in, line) || trim(line) != "$EndNodes") {
                throw GmshReadError("read_gmsh_msh: expected $EndNodes");
            }
            saw_nodes = true;
            continue;
        }

        if (tag == "$Elements") {
            if (!std::getline(in, line)) throw GmshReadError("read_gmsh_msh: truncated $Elements section");
            long long num_elements = 0;
            {
                std::istringstream count_line(line);
                count_line >> num_elements;
            }
            for (long long k = 0; k < num_elements; ++k) {
                if (!std::getline(in, line)) throw GmshReadError("read_gmsh_msh: $Elements section too short");
                std::istringstream elem_line(line);
                long long elm_number, elm_type, num_tags;
                elem_line >> elm_number >> elm_type >> num_tags;
                for (long long tg = 0; tg < num_tags; ++tg) {
                    long long discard_tag;
                    elem_line >> discard_tag;
                }
                if (elm_type == kGmshTetType) {
                    TetVerts tv{};
                    for (int v = 0; v < 4; ++v) {
                        long long node_id;
                        elem_line >> node_id;
                        auto it = gmsh_id_to_local.find(node_id);
                        if (it == gmsh_id_to_local.end()) {
                            throw GmshReadError("read_gmsh_msh: element references unknown node id " +
                                                 std::to_string(node_id));
                        }
                        tv[static_cast<std::size_t>(v)] = it->second;
                    }
                    mesh.tets.push_back(tv);
                }
                // Non-tet elements (points, lines, triangles used for boundary
                // tagging) are intentionally not stored here -- see the header
                // comment for why.
            }
            if (!std::getline(in, line) || trim(line) != "$EndElements") {
                throw GmshReadError("read_gmsh_msh: expected $EndElements");
            }
            saw_elements = true;
            continue;
        }
    }

    if (!saw_mesh_format) throw GmshReadError("read_gmsh_msh: missing $MeshFormat section");
    if (!saw_nodes) throw GmshReadError("read_gmsh_msh: missing $Nodes section");
    if (!saw_elements) throw GmshReadError("read_gmsh_msh: missing $Elements section");
    if (mesh.tets.empty()) {
        throw GmshReadError("read_gmsh_msh: no tetrahedral (elm-type 4) elements found in " + path);
    }

    mesh.build_topology();
    return mesh;
}

}  // namespace aphi_solver
