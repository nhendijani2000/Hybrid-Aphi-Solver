#include "aphi_solver/problem_binding.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <string>

namespace aphi_solver {

int BoundPort::side_of_tet(int t) const {
    const auto it = std::lower_bound(touching_tets.begin(), touching_tets.end(), t);
    if (it == touching_tets.end() || *it != t) return 0;
    return tet_side[static_cast<std::size_t>(it - touching_tets.begin())];
}

void scale_mesh_to_metres(Mesh& mesh, LengthUnit unit) {
    const double s = length_scale(unit);
    if (s == 1.0) return;
    for (Vec3& p : mesh.nodes) {
        p.x *= s;
        p.y *= s;
        p.z *= s;
    }
}

namespace {

// Binding errors have no line of their own -- they are about the mesh --
// but they are always *caused* by something written in the input file, so
// they carry that line. "<input>" as the path keeps the message honest when
// the Problem was built in code rather than parsed.
[[noreturn]] void fail(int line, const std::string& message) {
    throw InputError("<input>", line, message);
}

/// Resolves a name written in the input file to a physical-group tag of the
/// given dimension. A name that is all digits is taken as the tag itself --
/// meshes without a $PhysicalNames section can only be addressed that way.
int resolve_tag(const Mesh& mesh, const std::string& name, int dimension, int line,
                const std::string& what) {
    const bool numeric = !name.empty() &&
                         name.find_first_not_of("0123456789") == std::string::npos;
    if (numeric) return std::stoi(name);

    for (const PhysicalName& p : mesh.physical_names) {
        if (p.dimension == dimension && p.name == name) return p.tag;
    }

    // List what the mesh does offer: a name that is merely misspelled is by
    // far the likeliest cause, and the fix is then obvious.
    std::string available;
    for (const PhysicalName& p : mesh.physical_names) {
        if (p.dimension != dimension) continue;
        if (!available.empty()) available += ", ";
        available += "'" + p.name + "'";
    }
    if (available.empty()) available = "(the mesh names none)";
    fail(line, what + " '" + name + "' is not a " +
                   (dimension == 3 ? "Physical Volume" : "Physical Surface") +
                   " in the mesh. Available: " + available);
}

std::string join_names(const std::vector<std::string>& v) {
    std::string out;
    for (std::size_t i = 0; i < v.size(); ++i) {
        if (i) out += ", ";
        out += "'" + v[i] + "'";
    }
    return out;
}

std::vector<int> sorted_unique(std::vector<int> v) {
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
    return v;
}

/// The unit normal of face f, from its nodes in ascending order. Its SIGN is
/// arbitrary -- it follows the node numbering, not any geometric intent --
/// which is why every use below takes |n.d| or fixes the sign deliberately.
Vec3 face_normal(const Mesh& mesh, int f) {
    const std::array<int, 3>& v = mesh.faces[static_cast<std::size_t>(f)];
    const Vec3& a = mesh.nodes[static_cast<std::size_t>(v[0])];
    const Vec3& b = mesh.nodes[static_cast<std::size_t>(v[1])];
    const Vec3& c = mesh.nodes[static_cast<std::size_t>(v[2])];
    const Vec3 n = (b - a).cross(c - a);
    const double len = n.norm();
    if (len == 0.0) return Vec3(0.0, 0.0, 0.0);  // degenerate; caught elsewhere
    return n * (1.0 / len);
}

Vec3 face_centroid(const Mesh& mesh, int f) {
    const std::array<int, 3>& v = mesh.faces[static_cast<std::size_t>(f)];
    const Vec3 sum = mesh.nodes[static_cast<std::size_t>(v[0])] +
                     mesh.nodes[static_cast<std::size_t>(v[1])] +
                     mesh.nodes[static_cast<std::size_t>(v[2])];
    return sum * (1.0 / 3.0);
}

/// The vertex of tet t that is NOT on face f. Used to decide which side of
/// the face the tet lies on -- a test that needs no winding and no
/// tolerance, only that the tet is non-degenerate.
int opposite_vertex(const Mesh& mesh, int t, int f) {
    const std::array<int, 3>& fv = mesh.faces[static_cast<std::size_t>(f)];
    for (int v : mesh.tets[static_cast<std::size_t>(t)]) {
        if (v != fv[0] && v != fv[1] && v != fv[2]) return v;
    }
    return -1;
}

/// Collects the faces of one physical surface tag, resolved to global face
/// indices, plus the node and edge sets they span.
BoundSurface resolve_surface(const Mesh& mesh, const std::string& name, int tag, int line) {
    BoundSurface s;
    s.name = name;
    s.tag = tag;

    std::vector<int> nodes, edges;
    for (const TaggedFace& tf : mesh.tagged_boundary_faces) {
        if (tf.tag != tag) continue;
        const int f = mesh.find_face(tf.nodes[0], tf.nodes[1], tf.nodes[2]);
        if (f < 0) {
            fail(line, "surface '" + name + "' has a triangle that is not a face of any "
                       "tetrahedron -- the surface mesh and the volume mesh disagree");
        }
        s.faces.push_back(f);
        for (int v : tf.nodes) nodes.push_back(v);
        edges.push_back(mesh.find_edge(tf.nodes[0], tf.nodes[1]));
        edges.push_back(mesh.find_edge(tf.nodes[1], tf.nodes[2]));
        edges.push_back(mesh.find_edge(tf.nodes[0], tf.nodes[2]));
    }
    if (s.faces.empty()) {
        fail(line, "surface '" + name + "' (tag " + std::to_string(tag) +
                       ") has no triangles in the mesh");
    }
    s.faces = sorted_unique(std::move(s.faces));
    s.nodes = sorted_unique(std::move(nodes));
    s.edges = sorted_unique(std::move(edges));
    return s;
}

/// Edges used by exactly one of `faces` -- the surface's rim.
std::vector<int> rim_of(const Mesh& mesh, const std::vector<int>& faces) {
    std::map<int, int> use_count;
    for (int f : faces) {
        const std::array<int, 3>& v = mesh.faces[static_cast<std::size_t>(f)];
        ++use_count[mesh.find_edge(v[0], v[1])];
        ++use_count[mesh.find_edge(v[1], v[2])];
        ++use_count[mesh.find_edge(v[0], v[2])];
    }
    std::vector<int> rim;
    for (const auto& kv : use_count) {
        if (kv.second == 1) rim.push_back(kv.first);
    }
    return rim;
}

/// Union-find over the conductor tets, joining two that share a face. Used
/// for the DC checks: every conducting piece that carries a port needs
/// exactly one potential reference.
class Components {
public:
    explicit Components(int n) : parent_(static_cast<std::size_t>(n)) {
        for (int i = 0; i < n; ++i) parent_[static_cast<std::size_t>(i)] = i;
    }
    int find(int x) {
        while (parent_[static_cast<std::size_t>(x)] != x) {
            parent_[static_cast<std::size_t>(x)] =
                parent_[static_cast<std::size_t>(parent_[static_cast<std::size_t>(x)])];
            x = parent_[static_cast<std::size_t>(x)];
        }
        return x;
    }
    void unite(int a, int b) {
        a = find(a);
        b = find(b);
        if (a != b) parent_[static_cast<std::size_t>(a)] = b;
    }

private:
    std::vector<int> parent_;
};

}  // namespace

BoundProblem bind_to_mesh(const Problem& problem, const Mesh& mesh) {
    BoundProblem out;
    out.problem = problem;
    const int num_tets = mesh.num_tets();

    // ---- bodies -----------------------------------------------------------
    out.body_of_tet.assign(static_cast<std::size_t>(num_tets), -1);

    for (const Body& b : problem.bodies) {
        BoundBody bb;
        bb.name = b.name;
        bb.volume = b.volume;
        bb.tag = resolve_tag(mesh, b.volume, 3, b.line, "body '" + b.name + "': volume");
        bb.sigma = b.sigma;
        bb.eps_r = b.eps_r;
        bb.mu_r = b.mu_r;

        const int index = static_cast<int>(out.bodies.size());
        for (int t = 0; t < num_tets; ++t) {
            if (mesh.tet_tags[static_cast<std::size_t>(t)] != bb.tag) continue;
            if (out.body_of_tet[static_cast<std::size_t>(t)] >= 0) {
                fail(b.line, "volume '" + bb.volume + "' is already claimed by body '" +
                                 out.bodies[static_cast<std::size_t>(
                                                out.body_of_tet[static_cast<std::size_t>(t)])]
                                     .name +
                                 "'");
            }
            out.body_of_tet[static_cast<std::size_t>(t)] = index;
            bb.tets.push_back(t);
        }
        if (bb.tets.empty()) {
            fail(b.line, "body '" + bb.name + "': volume '" + bb.volume + "' (tag " +
                             std::to_string(bb.tag) + ") has no tetrahedra in the mesh");
        }
        out.bodies.push_back(std::move(bb));
    }

    // Every tet must belong to some body: an unclaimed region would be
    // assembled with no material at all, which is a silent wrong answer
    // rather than a visible failure.
    for (int t = 0; t < num_tets; ++t) {
        if (out.body_of_tet[static_cast<std::size_t>(t)] < 0) {
            const int tag = mesh.tet_tags[static_cast<std::size_t>(t)];
            const std::string name = mesh.physical_name(3, tag);
            fail(0, "the mesh has a region no [Body ...] claims: tag " + std::to_string(tag) +
                        (name.empty() ? "" : " ('" + name + "')") +
                        ". Every volume needs a material, or it would assemble as nothing");
        }
    }

    // ---- where Phi lives --------------------------------------------------
    const bool conductors_only = phi_on_conductors_only(problem);
    out.phi_tet.assign(static_cast<std::size_t>(num_tets), true);
    if (conductors_only) {
        for (int t = 0; t < num_tets; ++t) {
            out.phi_tet[static_cast<std::size_t>(t)] =
                out.bodies[static_cast<std::size_t>(out.body_of_tet[static_cast<std::size_t>(t)])]
                    .is_conductor();
        }
    }

    // ---- ports ------------------------------------------------------------
    for (const Port& p : problem.ports) {
        BoundPort bp;
        bp.name = p.name;
        bp.type = p.type;
        bp.amplitude = port_amplitude(p);

        // A port may name several surfaces; they merge into one.
        BoundSurface merged;
        merged.name = p.surface.front();
        std::vector<int> faces, nodes, edges;
        for (const std::string& name : p.surface) {
            const int tag = resolve_tag(mesh, name, 2, p.line, "port '" + p.name + "': surface");
            const BoundSurface s = resolve_surface(mesh, name, tag, p.line);
            merged.tag = s.tag;
            faces.insert(faces.end(), s.faces.begin(), s.faces.end());
            nodes.insert(nodes.end(), s.nodes.begin(), s.nodes.end());
            edges.insert(edges.end(), s.edges.begin(), s.edges.end());
        }
        merged.faces = sorted_unique(std::move(faces));
        merged.nodes = sorted_unique(std::move(nodes));
        merged.edges = sorted_unique(std::move(edges));
        bp.surface = std::move(merged);

        // Boundary vs interior, checked against what the file declared. A
        // surface believed to be on the boundary but embedded in the volume
        // is exactly the mistake worth catching loudly.
        int on_boundary = 0;
        for (int f : bp.surface.faces) {
            if (mesh.is_boundary_face(f)) ++on_boundary;
        }
        const int n_faces = static_cast<int>(bp.surface.faces.size());
        if (aphi_solver::is_internal(p.type)) {
            if (on_boundary > 0) {
                fail(p.line, "port '" + p.name + "' is declared internal, but " +
                                 std::to_string(on_boundary) + " of its " +
                                 std::to_string(n_faces) +
                                 " faces lie on the domain boundary. An internal cut needs a "
                                 "tetrahedron on each side");
            }
        } else if (on_boundary != n_faces) {
            fail(p.line, "port '" + p.name + "' is declared a boundary port, but " +
                             std::to_string(n_faces - on_boundary) + " of its " +
                             std::to_string(n_faces) +
                             " faces are interior to the mesh");
        }

        // Orientation.
        const Vec3 n0 = face_normal(mesh, bp.surface.faces.front());
        if (aphi_solver::is_internal(p.type)) {
            // The cut must be planar, or "the" normal does not exist. The
            // magnitude is what is tested: each face's normal takes its sign
            // from node numbering, and the sign is not needed -- the plus
            // side comes from the opposite-vertex test below.
            for (int f : bp.surface.faces) {
                const double alignment = std::abs(face_normal(mesh, f).dot(n0));
                if (std::abs(alignment - 1.0) > 1e-8) {
                    fail(p.line, "port '" + p.name +
                                     "' is not planar: its faces do not share a normal "
                                     "(worst |n.n0| = " +
                                     std::to_string(alignment) +
                                     "). Only planar cuts are supported");
                }
            }
            if (p.current_direction.has_value()) {
                bp.direction = resolve_cut_direction(*p.current_direction, n0);
                bp.direction_from_hint = true;
            } else {
                bp.direction = canonical_orientation(n0);
                bp.direction_from_hint = false;
            }

            // Which adjacent tet of each face is on the plus side.
            bp.plus_side_tet.reserve(bp.surface.faces.size());
            for (int f : bp.surface.faces) {
                const FaceTets& ft = mesh.face_tets[static_cast<std::size_t>(f)];
                const Vec3 c = face_centroid(mesh, f);
                int plus = -1;
                for (int i = 0; i < ft.count; ++i) {
                    const int t = ft.tets[static_cast<std::size_t>(i)];
                    const int apex = opposite_vertex(mesh, t, f);
                    if (apex < 0) continue;
                    if ((mesh.nodes[static_cast<std::size_t>(apex)] - c).dot(bp.direction) > 0.0) {
                        plus = t;
                    }
                }
                if (plus < 0) {
                    fail(p.line, "port '" + p.name +
                                     "': a cut face has no tetrahedron on its plus side");
                }
                bp.plus_side_tet.push_back(plus);
            }
            bp.rim_edges = rim_of(mesh, bp.surface.faces);
        } else {
            // A boundary terminal: positive current enters the domain, so d
            // is the inward normal. The one adjacent tet is the inside, so
            // its opposite vertex says which way that is.
            const int f = bp.surface.faces.front();
            const FaceTets& ft = mesh.face_tets[static_cast<std::size_t>(f)];
            const int apex = opposite_vertex(mesh, ft.tets[0], f);
            const Vec3 inward = mesh.nodes[static_cast<std::size_t>(apex)] - face_centroid(mesh, f);
            bp.direction = inward.dot(n0) > 0.0 ? n0 : n0 * -1.0;
            bp.direction_from_hint = false;
            bp.plus_side_tet.assign(bp.surface.faces.size(), -1);
        }

        out.ports.push_back(std::move(bp));
    }

    // Two terminals sharing a node are a short circuit, which would make the
    // ports one port without saying so.
    for (std::size_t i = 0; i < out.ports.size(); ++i) {
        for (std::size_t j = i + 1; j < out.ports.size(); ++j) {
            std::vector<int> shared;
            std::set_intersection(out.ports[i].surface.nodes.begin(),
                                  out.ports[i].surface.nodes.end(),
                                  out.ports[j].surface.nodes.begin(),
                                  out.ports[j].surface.nodes.end(), std::back_inserter(shared));
            if (!shared.empty()) {
                fail(problem.ports[j].line,
                     "ports '" + out.ports[i].name + "' and '" + out.ports[j].name + "' share " +
                         std::to_string(shared.size()) +
                         " mesh node(s): they are short-circuited together");
            }
        }
    }

    // ---- checks that need Phi's support ----------------------------------
    // An internal cut's rim must lie on the boundary of Phi's support. At DC
    // that is the conductor surface, and a rim inside the conductor means
    // the cut does not fully span it -- current flows around the uncut part
    // and the port is partly shorted, with no symptom but a wrong answer.
    // At full wave Phi reaches into the dielectric, so the rim is always
    // interior: that is the delta gap, deliberate and merely warned about.
    std::vector<std::vector<int>> tets_of_edge(static_cast<std::size_t>(mesh.num_edges()));
    for (int t = 0; t < num_tets; ++t) {
        for (int e : mesh.tet_edges[static_cast<std::size_t>(t)]) {
            tets_of_edge[static_cast<std::size_t>(e)].push_back(t);
        }
    }
    std::vector<bool> edge_on_boundary(static_cast<std::size_t>(mesh.num_edges()), false);
    for (int f = 0; f < mesh.num_faces(); ++f) {
        if (!mesh.is_boundary_face(f)) continue;
        const std::array<int, 3>& v = mesh.faces[static_cast<std::size_t>(f)];
        edge_on_boundary[static_cast<std::size_t>(mesh.find_edge(v[0], v[1]))] = true;
        edge_on_boundary[static_cast<std::size_t>(mesh.find_edge(v[1], v[2]))] = true;
        edge_on_boundary[static_cast<std::size_t>(mesh.find_edge(v[0], v[2]))] = true;
    }

    for (std::size_t i = 0; i < out.ports.size(); ++i) {
        const BoundPort& bp = out.ports[i];
        if (!bp.is_internal()) continue;
        const int line = problem.ports[i].line;

        int interior_rim = 0;
        for (int e : bp.rim_edges) {
            if (edge_on_boundary[static_cast<std::size_t>(e)]) continue;
            bool leaves_support = false;
            for (int t : tets_of_edge[static_cast<std::size_t>(e)]) {
                if (!out.phi_tet[static_cast<std::size_t>(t)]) leaves_support = true;
            }
            if (!leaves_support) ++interior_rim;
        }

        if (interior_rim > 0) {
            if (conductors_only) {
                fail(line, "port '" + bp.name + "': the cut does not fully span its conductor -- " +
                               std::to_string(interior_rim) +
                               " of its rim edges are surrounded by conductor. Current would flow "
                               "around the uncut part and the port would be partly shorted");
            }
            out.warnings.push_back(
                "port '" + bp.name + "' is a delta gap: its rim lies inside the region where Phi "
                "lives, so the field there is singular and the gap capacitance depends on mesh "
                "refinement. This is the standard idealization, not an error.");
        }
    }

    // ---- conduction paths -------------------------------------------------
    // Built from sigma alone, in every regime: a path is a property of the
    // material, not of the formulation. Two bodies that touch are one path.
    {
        Components comp(num_tets);
        std::vector<bool> conducting(static_cast<std::size_t>(num_tets), false);
        for (int t = 0; t < num_tets; ++t) {
            conducting[static_cast<std::size_t>(t)] =
                out.bodies[static_cast<std::size_t>(out.body_of_tet[static_cast<std::size_t>(t)])]
                    .is_conductor();
        }
        for (int f = 0; f < mesh.num_faces(); ++f) {
            const FaceTets& ft = mesh.face_tets[static_cast<std::size_t>(f)];
            if (ft.count != 2) continue;
            const int a = ft.tets[0], b = ft.tets[1];
            if (conducting[static_cast<std::size_t>(a)] && conducting[static_cast<std::size_t>(b)]) {
                comp.unite(a, b);
            }
        }

        out.path_of_tet.assign(static_cast<std::size_t>(num_tets), -1);
        std::map<int, int> root_to_path;
        for (int t = 0; t < num_tets; ++t) {
            if (!conducting[static_cast<std::size_t>(t)]) continue;
            const int root = comp.find(t);
            auto it = root_to_path.find(root);
            if (it == root_to_path.end()) {
                it = root_to_path.emplace(root, static_cast<int>(out.conduction_paths.size())).first;
                out.conduction_paths.emplace_back();
            }
            out.path_of_tet[static_cast<std::size_t>(t)] = it->second;
            out.conduction_paths[static_cast<std::size_t>(it->second)].tets.push_back(t);
        }
        for (ConductionPath& path : out.conduction_paths) {
            std::vector<int> bodies;
            for (int t : path.tets) bodies.push_back(out.body_of_tet[static_cast<std::size_t>(t)]);
            path.bodies = sorted_unique(std::move(bodies));
        }

        // Which paths each port touches, and whether it fixes a potential
        // there. A voltage port does; so does an internal current port,
        // which holds its minus side at 0 V. A boundary current port does
        // not -- it prescribes a flux, not a level.
        for (std::size_t i = 0; i < out.ports.size(); ++i) {
            const BoundPort& bp = out.ports[i];
            const bool is_reference = !is_current_driven(bp.type) || bp.is_internal();
            std::set<int> paths;
            for (int f : bp.surface.faces) {
                const FaceTets& ft = mesh.face_tets[static_cast<std::size_t>(f)];
                for (int k = 0; k < ft.count; ++k) {
                    const int p = out.path_of_tet[static_cast<std::size_t>(
                        ft.tets[static_cast<std::size_t>(k)])];
                    if (p >= 0) paths.insert(p);
                }
            }
            if (paths.empty() && conductors_only) {
                fail(problem.ports[i].line,
                     "port '" + bp.name +
                         "' touches no conductor, so at DC it can carry no current");
            }
            for (int p : paths) {
                ConductionPath& path = out.conduction_paths[static_cast<std::size_t>(p)];
                path.ports.push_back(bp.name);
                if (is_reference) ++path.reference_count;
            }
        }

        // At DC (and in the reduced variant) Phi lives only on the
        // conductors, so each path is its own Phi problem and needs its own
        // reference. At full wave (sigma + j*omega*eps) couples everything
        // through the dielectric, the domain is one Phi problem, and a
        // single global reference -- already required at parse time --
        // suffices.
        //
        // The rule is "at least one", not "exactly one". Two voltage ports
        // on one conductor is a voltage-driven resistor: a well-posed
        // Dirichlet problem whose current falls out as a result. An earlier
        // draft of this code rejected that, and the test suite defended the
        // mistake; the cylinder passed either way, which is why it went
        // unnoticed.
        if (conductors_only) {
            for (ConductionPath& path : out.conduction_paths) {
                if (path.reference_count > 0) continue;

                if (!path.ports.empty()) {
                    fail(0, "a conductor carries port(s) " + join_names(path.ports) +
                                " but nothing fixes its potential, so Phi there is determined only "
                                "up to a constant. Add a voltage port on it");
                }

                // Floating: no port, no reference. Legitimate (an
                // unconnected shield or pad) but also what a forgotten
                // connection looks like, so it is said once. Pinning the
                // lowest-numbered node is deterministic on purpose -- "any
                // node" is fine, "a different node between runs" is not.
                path.is_floating = true;
                int lowest = -1;
                for (int t : path.tets) {
                    for (int v : mesh.tets[static_cast<std::size_t>(t)]) {
                        if (lowest < 0 || v < lowest) lowest = v;
                    }
                }
                path.pin_node = lowest;
                out.warnings.push_back(
                    "a conductor of " + std::to_string(path.tets.size()) +
                    " tets is floating: no port touches it and nothing fixes its potential. Node " +
                    std::to_string(lowest) +
                    " will be pinned to 0 V to keep the matrix non-singular. It carries no DC "
                    "current, so the value is arbitrary -- but check this is not a missing "
                    "connection.");
            }
        }
    }

    // ---- the tree-cotree gauge -------------------------------------------
    // Run here because the DOF map cannot classify an edge without it. The
    // mask is every edge on the outer boundary, which is what
    // `outer = flux_tangential` means and the only option supported today.
    out.dirichlet_edge = boundary_edge_mask(mesh);
    out.gauge = build_tree_cotree(mesh, out.dirichlet_edge);

    // ---- which side of each cut every touching tet is on ------------------
    for (BoundPort& bp : out.ports) {
        if (!bp.is_internal()) continue;

        // A point on the cut plane, and the normal pointing minus -> plus.
        const Vec3 origin = face_centroid(mesh, bp.surface.faces.front());

        std::vector<int> touching;
        for (int f : bp.surface.faces) {
            const FaceTets& ft = mesh.face_tets[static_cast<std::size_t>(f)];
            for (int k = 0; k < ft.count; ++k) touching.push_back(ft.tets[static_cast<std::size_t>(k)]);
        }
        // Every tet sharing a NODE with the cut, not merely a face: a tet
        // touching at one vertex still reads a value there.
        std::set<int> cut_nodes(bp.surface.nodes.begin(), bp.surface.nodes.end());
        for (int t = 0; t < num_tets; ++t) {
            for (int v : mesh.tets[static_cast<std::size_t>(t)]) {
                if (cut_nodes.count(v)) {
                    touching.push_back(t);
                    break;
                }
            }
        }
        bp.touching_tets = sorted_unique(std::move(touching));

        bp.tet_side.reserve(bp.touching_tets.size());
        for (int t : bp.touching_tets) {
            Vec3 centroid(0.0, 0.0, 0.0);
            for (int v : mesh.tets[static_cast<std::size_t>(t)]) {
                centroid = centroid + mesh.nodes[static_cast<std::size_t>(v)];
            }
            centroid = centroid * 0.25;
            const double signed_distance = (centroid - origin).dot(bp.direction);
            bp.tet_side.push_back(signed_distance > 0.0 ? +1 : -1);
        }
    }

    return out;
}

}  // namespace aphi_solver
