#include "aphi_solver/gmsh_reader.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <map>
#include <sstream>
#include <unordered_map>
#include <utility>
#include <vector>

namespace aphi_solver {

namespace {

constexpr int kGmshTetType = 4;       // Gmsh elm-type 4 == 4-node tetrahedron.
constexpr int kGmshTriangleType = 2;  // Gmsh elm-type 2 == 3-node triangle.

std::string trim(const std::string& s) {
    const auto a = s.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) return "";
    const auto b = s.find_last_not_of(" \t\r\n");
    return s.substr(a, b - a + 1);
}

std::string next_line(std::istream& in, const char* what) {
    std::string line;
    if (!std::getline(in, line)) {
        throw GmshReadError(std::string("read_gmsh_msh: truncated ") + what + " section");
    }
    return line;
}

/// Reads past a section this reader does not consume, up to its `$End...`
/// line. Done explicitly rather than by ignoring unrecognized lines, so a
/// section whose *contents* happen to start with '$' cannot be mistaken for
/// a new section header.
void skip_section(std::istream& in, const std::string& open_tag) {
    const std::string end_tag = "$End" + open_tag.substr(1);
    std::string line;
    while (std::getline(in, line)) {
        if (trim(line) == end_tag) return;
    }
    throw GmshReadError("read_gmsh_msh: unterminated section " + open_tag);
}

void expect_end(std::istream& in, const char* end_tag) {
    std::string line;
    if (!std::getline(in, line) || trim(line) != end_tag) {
        throw GmshReadError(std::string("read_gmsh_msh: expected ") + end_tag);
    }
}

/// Accumulates the mesh while parsing, and owns the Gmsh-id -> 0-based
/// remapping both format versions need. Node ids need not be contiguous or
/// 1-based in either format.
struct Accumulator {
    Mesh mesh;
    std::unordered_map<long long, int> id_to_local;

    void add_node(long long gmsh_id, double x, double y, double z) {
        const auto inserted = id_to_local.emplace(gmsh_id, static_cast<int>(mesh.nodes.size()));
        if (!inserted.second) {
            throw GmshReadError("read_gmsh_msh: duplicate node id " + std::to_string(gmsh_id));
        }
        mesh.nodes.emplace_back(x, y, z);
    }

    int local(long long gmsh_id) const {
        auto it = id_to_local.find(gmsh_id);
        if (it == id_to_local.end()) {
            throw GmshReadError("read_gmsh_msh: element references unknown node id " +
                                 std::to_string(gmsh_id));
        }
        return it->second;
    }

    /// Consumes an element's node list from `rest` and stores it if it is a
    /// type this solver uses. Shared by both format versions: they differ in
    /// where the element's type and physical tag come from, not in what an
    /// element means once you have them.
    void add_element(long long elm_type, int physical_tag, std::istream& rest) {
        if (elm_type == kGmshTetType) {
            TetVerts tv{};
            for (int v = 0; v < 4; ++v) {
                long long id = 0;
                rest >> id;
                tv[static_cast<std::size_t>(v)] = local(id);
            }
            mesh.tets.push_back(tv);
            mesh.tet_tags.push_back(physical_tag);
        } else if (elm_type == kGmshTriangleType) {
            // Surface elements are how Gmsh marks PEC walls, ports and the
            // outer truncation boundary. Stored sorted, matching
            // Mesh::faces' canonical ascending orientation, so
            // Mesh::find_face resolves them directly. No meaning is attached
            // to the tag here -- that is the input file's job.
            TaggedFace face;
            for (int v = 0; v < 3; ++v) {
                long long id = 0;
                rest >> id;
                face.nodes[static_cast<std::size_t>(v)] = local(id);
            }
            std::sort(face.nodes.begin(), face.nodes.end());
            face.tag = physical_tag;
            mesh.tagged_boundary_faces.push_back(face);
        }
        // Other element types (points, lines, higher-order cells) are
        // skipped: nothing in this solver consumes them.
    }
};

/// `$PhysicalNames` -- identical in 2.2 and 4.x:
///     numNames
///     dimension physicalTag "name"
/// The name is quoted and may contain spaces, so it is taken between the
/// first and last quote rather than read as a whitespace-delimited token.
void parse_physical_names(std::istream& in, Mesh& mesh) {
    long long count = 0;
    {
        std::istringstream header(next_line(in, "$PhysicalNames"));
        header >> count;
    }
    for (long long k = 0; k < count; ++k) {
        const std::string line = next_line(in, "$PhysicalNames");
        std::istringstream iss(line);
        PhysicalName entry;
        long long tag = 0;
        iss >> entry.dimension >> tag;
        entry.tag = static_cast<int>(tag);

        std::string rest;
        std::getline(iss, rest);
        const auto q1 = rest.find('"');
        const auto q2 = rest.rfind('"');
        entry.name = (q1 != std::string::npos && q2 != std::string::npos && q2 > q1)
                          ? rest.substr(q1 + 1, q2 - q1 - 1)
                          : trim(rest);
        mesh.physical_names.push_back(entry);
    }
    expect_end(in, "$EndPhysicalNames");
}

// ---------------------------------------------------------------- format 2.2

void parse_nodes_v2(std::istream& in, Accumulator& acc) {
    long long num_nodes = 0;
    {
        std::istringstream header(next_line(in, "$Nodes"));
        header >> num_nodes;
    }
    acc.mesh.nodes.reserve(static_cast<std::size_t>(num_nodes));
    for (long long k = 0; k < num_nodes; ++k) {
        std::istringstream iss(next_line(in, "$Nodes"));
        long long id = 0;
        double x = 0, y = 0, z = 0;
        iss >> id >> x >> y >> z;
        acc.add_node(id, x, y, z);
    }
    expect_end(in, "$EndNodes");
}

void parse_elements_v2(std::istream& in, Accumulator& acc) {
    long long num_elements = 0;
    {
        std::istringstream header(next_line(in, "$Elements"));
        header >> num_elements;
    }
    for (long long k = 0; k < num_elements; ++k) {
        std::istringstream iss(next_line(in, "$Elements"));
        long long elm_number = 0, elm_type = 0, num_tags = 0;
        iss >> elm_number >> elm_type >> num_tags;

        // In 2.2 the tags sit on the element line itself, and by convention
        // the first is the physical-group id (the second, the elementary
        // geometrical entity, is unused here). -1 means no tags at all.
        int physical_tag = -1;
        for (long long tg = 0; tg < num_tags; ++tg) {
            long long value = 0;
            iss >> value;
            if (tg == 0) physical_tag = static_cast<int>(value);
        }
        acc.add_element(elm_type, physical_tag, iss);
    }
    expect_end(in, "$EndElements");
}

// ---------------------------------------------------------------- format 4.1

using EntityPhysical = std::map<std::pair<int, int>, int>;  // (dim, entityTag) -> physical tag

/// `$Entities` -- the reason 4.x needs more than a reformat of 2.2. In this
/// format an element line carries NO physical tag; it carries the geometric
/// entity it belongs to, and the entity-to-physical-group mapping lives
/// here. Without this section, every element in a 4.x file would come back
/// untagged.
///
///     numPoints numCurves numSurfaces numVolumes
///     pointTag   x y z                numPhysicalTags tags...
///     curveTag   minX..maxZ           numPhysicalTags tags... numBounding...
///     surfaceTag minX..maxZ           numPhysicalTags tags... numBounding...
///     volumeTag  minX..maxZ           numPhysicalTags tags... numBounding...
///
/// Points carry three coordinates; the other three kinds carry a six-value
/// bounding box. Trailing bounding-entity lists are not needed here.
void parse_entities_v4(std::istream& in, EntityPhysical& entity_physical) {
    long long counts[4] = {0, 0, 0, 0};
    {
        std::istringstream header(next_line(in, "$Entities"));
        header >> counts[0] >> counts[1] >> counts[2] >> counts[3];
    }
    for (int dim = 0; dim < 4; ++dim) {
        const int coord_values = (dim == 0) ? 3 : 6;  // point xyz vs bounding box
        for (long long k = 0; k < counts[dim]; ++k) {
            std::istringstream iss(next_line(in, "$Entities"));
            long long entity_tag = 0;
            iss >> entity_tag;
            for (int c = 0; c < coord_values; ++c) {
                double ignored = 0;
                iss >> ignored;
            }
            long long num_physical = 0;
            iss >> num_physical;
            int first_physical = -1;
            for (long long p = 0; p < num_physical; ++p) {
                long long value = 0;
                iss >> value;
                if (p == 0) first_physical = static_cast<int>(value);
            }
            entity_physical[{dim, static_cast<int>(entity_tag)}] = first_physical;
        }
    }
    expect_end(in, "$EndEntities");
}

/// `$Nodes` (4.1):
///     numEntityBlocks numNodes minNodeTag maxNodeTag
///     entityDim entityTag parametric numNodesInBlock
///       nodeTag        x numNodesInBlock
///       x y z          x numNodesInBlock
///
/// The split of tags from coordinates is what distinguishes 4.1 from 4.0,
/// and why 4.0 is rejected rather than guessed at.
void parse_nodes_v4(std::istream& in, Accumulator& acc) {
    long long num_blocks = 0, num_nodes = 0, min_tag = 0, max_tag = 0;
    {
        std::istringstream header(next_line(in, "$Nodes"));
        header >> num_blocks >> num_nodes >> min_tag >> max_tag;
    }
    acc.mesh.nodes.reserve(static_cast<std::size_t>(num_nodes));

    for (long long b = 0; b < num_blocks; ++b) {
        long long entity_dim = 0, entity_tag = 0, parametric = 0, block_nodes = 0;
        {
            std::istringstream header(next_line(in, "$Nodes"));
            header >> entity_dim >> entity_tag >> parametric >> block_nodes;
        }
        std::vector<long long> ids(static_cast<std::size_t>(block_nodes));
        for (long long k = 0; k < block_nodes; ++k) {
            std::istringstream iss(next_line(in, "$Nodes"));
            iss >> ids[static_cast<std::size_t>(k)];
        }
        for (long long k = 0; k < block_nodes; ++k) {
            std::istringstream iss(next_line(in, "$Nodes"));
            double x = 0, y = 0, z = 0;
            iss >> x >> y >> z;
            acc.add_node(ids[static_cast<std::size_t>(k)], x, y, z);
        }
    }
    expect_end(in, "$EndNodes");
}

/// `$Elements` (4.1):
///     numEntityBlocks numElements minElementTag maxElementTag
///     entityDim entityTag elementType numElementsInBlock
///       elementTag nodeTag...   x numElementsInBlock
///
/// The element type is per *block*, and the physical tag comes from the
/// block's entity via `$Entities`.
void parse_elements_v4(std::istream& in, Accumulator& acc, const EntityPhysical& entity_physical) {
    long long num_blocks = 0, num_elements = 0, min_tag = 0, max_tag = 0;
    {
        std::istringstream header(next_line(in, "$Elements"));
        header >> num_blocks >> num_elements >> min_tag >> max_tag;
    }
    for (long long b = 0; b < num_blocks; ++b) {
        long long entity_dim = 0, entity_tag = 0, elm_type = 0, block_elements = 0;
        {
            std::istringstream header(next_line(in, "$Elements"));
            header >> entity_dim >> entity_tag >> elm_type >> block_elements;
        }
        int physical_tag = -1;
        const auto it = entity_physical.find({static_cast<int>(entity_dim), static_cast<int>(entity_tag)});
        if (it != entity_physical.end()) physical_tag = it->second;

        for (long long k = 0; k < block_elements; ++k) {
            std::istringstream iss(next_line(in, "$Elements"));
            long long element_tag = 0;
            iss >> element_tag;
            acc.add_element(elm_type, physical_tag, iss);
        }
    }
    expect_end(in, "$EndElements");
}

// ------------------------------------------------------------------ validate

/// Structural checks that a malformed mesh should fail loudly at load time
/// rather than somewhere deep in assembly. Everything here is cheap
/// (O(tets) / O(faces)) and none of it depends on the file format.
void validate_mesh(const Mesh& mesh, const std::string& path) {
    if (mesh.tets.empty()) {
        throw GmshReadError("read_gmsh_msh: no tetrahedral (elm-type 4) elements found in " + path);
    }

    for (int t = 0; t < mesh.num_tets(); ++t) {
        const TetVerts& tv = mesh.tets[static_cast<std::size_t>(t)];
        for (int a = 0; a < 4; ++a) {
            for (int b = a + 1; b < 4; ++b) {
                if (tv[static_cast<std::size_t>(a)] == tv[static_cast<std::size_t>(b)]) {
                    throw GmshReadError("read_gmsh_msh: tet " + std::to_string(t) +
                                         " uses the same node twice (node index " +
                                         std::to_string(tv[static_cast<std::size_t>(a)]) + ")");
                }
            }
        }

        // Degeneracy is tested against the tet's own size, not an absolute
        // volume: a mesh in metres with micron features has legitimately
        // tiny tets, and any fixed threshold would either reject those or
        // miss a flat tet in a coarse mesh. |V| / h^3 is dimensionless.
        double h2 = 0.0;
        for (int a = 0; a < 4; ++a) {
            for (int b = a + 1; b < 4; ++b) {
                const Vec3 d = mesh.nodes[static_cast<std::size_t>(tv[static_cast<std::size_t>(a)])] -
                                mesh.nodes[static_cast<std::size_t>(tv[static_cast<std::size_t>(b)])];
                h2 = std::max(h2, d.dot(d));
            }
        }
        const double h = std::sqrt(h2);
        const double volume = std::abs(mesh.signed_tet_volume(t));
        if (volume <= 1e-12 * h * h * h) {
            throw GmshReadError("read_gmsh_msh: tet " + std::to_string(t) +
                                 " is degenerate (volume " + std::to_string(volume) +
                                 " vs. longest edge " + std::to_string(h) + ")");
        }
    }

    for (int f = 0; f < mesh.num_faces(); ++f) {
        const int count = mesh.face_tets[static_cast<std::size_t>(f)].count;
        if (count > 2) {
            throw GmshReadError("read_gmsh_msh: face " + std::to_string(f) + " is shared by " +
                                 std::to_string(count) +
                                 " tets; a conforming tetrahedral mesh allows at most 2");
        }
    }
}

}  // namespace

Mesh read_gmsh_msh(const std::string& path) {
    std::ifstream in(path);
    if (!in) {
        throw GmshReadError("read_gmsh_msh: could not open file: " + path);
    }

    Accumulator acc;
    EntityPhysical entity_physical;
    int major = 0, minor = 0;
    bool saw_format = false, saw_nodes = false, saw_elements = false;

    std::string line;
    while (std::getline(in, line)) {
        const std::string section = trim(line);
        if (section.empty() || section[0] != '$') continue;

        if (section == "$MeshFormat") {
            std::istringstream fmt(next_line(in, "$MeshFormat"));
            std::string version_str;
            int file_type = 0, data_size = 0;
            fmt >> version_str >> file_type >> data_size;
            {
                std::istringstream vs(version_str);
                char dot = '.';
                vs >> major >> dot >> minor;
            }
            if (file_type != 0) {
                throw GmshReadError("read_gmsh_msh: only ASCII .msh files are supported (file_type 0)");
            }
            const bool supported = (major == 2) || (major == 4 && minor >= 1);
            if (!supported) {
                throw GmshReadError(
                    "read_gmsh_msh: unsupported .msh version " + version_str +
                    " (supported: 2.x, and 4.1 or later 4.x; MSH 4.0 stores $Nodes differently from "
                    "4.1 -- re-export as 4.1 or 2.2)");
            }
            expect_end(in, "$EndMeshFormat");
            saw_format = true;
        } else if (section == "$PhysicalNames") {
            parse_physical_names(in, acc.mesh);
        } else if (section == "$Entities") {
            // Gmsh writes $Entities before $Elements, which is what lets a
            // single sequential pass resolve each element's physical group.
            if (major == 4) {
                parse_entities_v4(in, entity_physical);
            } else {
                skip_section(in, section);
            }
        } else if (section == "$Nodes") {
            if (!saw_format) throw GmshReadError("read_gmsh_msh: $Nodes before $MeshFormat");
            if (major == 2) {
                parse_nodes_v2(in, acc);
            } else {
                parse_nodes_v4(in, acc);
            }
            saw_nodes = true;
        } else if (section == "$Elements") {
            if (!saw_format) throw GmshReadError("read_gmsh_msh: $Elements before $MeshFormat");
            if (major == 2) {
                parse_elements_v2(in, acc);
            } else {
                parse_elements_v4(in, acc, entity_physical);
            }
            saw_elements = true;
        } else {
            skip_section(in, section);
        }
    }

    if (!saw_format) throw GmshReadError("read_gmsh_msh: missing $MeshFormat section");
    if (!saw_nodes) throw GmshReadError("read_gmsh_msh: missing $Nodes section");
    if (!saw_elements) throw GmshReadError("read_gmsh_msh: missing $Elements section");

    acc.mesh.build_topology();
    validate_mesh(acc.mesh, path);
    return acc.mesh;
}

}  // namespace aphi_solver
