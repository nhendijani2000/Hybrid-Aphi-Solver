#include "aphi_solver/postprocess.hpp"

#include <array>
#include <charconv>
#include <cmath>
#include <chrono>
#include <cstring>
#include <stdexcept>
#include <string>
#include <utility>

#include "aphi_solver/basis_functions.hpp"
#include "aphi_solver/constants.hpp"
#include "aphi_solver/version.hpp"

namespace aphi_solver {

using Complex = std::complex<double>;

const char* field_set_name(FieldSet which) {
    switch (which) {
        case FieldSet::Potential: return "potential";
        case FieldSet::A: return "A";
        case FieldSet::B: return "B";
        case FieldSet::H: return "H";
        case FieldSet::E: return "E";
        case FieldSet::All: break;
    }
    return "fields";
}

std::complex<double> Solution::phi_scale() const {
    // Only ScaledPhi substitutes. Natural and RowScaled leave the Phi unknown as
    // Phi itself -- RowScaled scales the Phi ROWS, which changes the equations
    // and not the unknown.
    return conditioning == Conditioning::ScaledPhi ? Complex(0.0, omega) : Complex(1.0, 0.0);
}

NodalPotential potential_at_nodes(const Mesh& mesh, const BoundProblem& bound, const DofMap& dofs,
                                  const Solution& solution) {
    if (static_cast<int>(solution.x.size()) != dofs.num_total) {
        throw std::invalid_argument(
            "potential_at_nodes: the solution has " + std::to_string(solution.x.size()) +
            " entries for a DOF map of " + std::to_string(dofs.num_total) +
            ". They must describe the same problem.");
    }

    NodalPotential out;
    out.num_vertices = mesh.num_nodes();
    out.num_edges = mesh.num_edges();
    const int n = dofs.num_p2_nodes;
    out.position.resize(static_cast<std::size_t>(n));
    out.value.assign(static_cast<std::size_t>(n), Complex(0.0, 0.0));
    out.status.assign(static_cast<std::size_t>(n), PhiStatus::Absent);

    // Positions: the vertices as they are, then each edge's midpoint. The node
    // numbering is DofMap's -- vertex_p2(v) = v and edge_p2(e) = num_nodes + e --
    // so this stays in step with `phi_state` without a second convention.
    for (int v = 0; v < mesh.num_nodes(); ++v) {
        out.position[static_cast<std::size_t>(dofs.vertex_p2(v))] =
            mesh.nodes[static_cast<std::size_t>(v)];
    }
    for (int e = 0; e < mesh.num_edges(); ++e) {
        const std::pair<int, int>& ends = mesh.edges[static_cast<std::size_t>(e)];
        const Vec3& p = mesh.nodes[static_cast<std::size_t>(ends.first)];
        const Vec3& q = mesh.nodes[static_cast<std::size_t>(ends.second)];
        out.position[static_cast<std::size_t>(dofs.edge_p2(e))] =
            Vec3{0.5 * (p.x + q.x), 0.5 * (p.y + q.y), 0.5 * (p.z + q.z)};
    }

    const Complex scale = solution.phi_scale();

    for (int i = 0; i < n; ++i) {
        const std::size_t u = static_cast<std::size_t>(i);
        switch (dofs.phi_state[u]) {
            case PhiDof::Free:
                out.value[u] = scale * solution.x[static_cast<std::size_t>(dofs.phi_index[u])];
                out.status[u] = PhiStatus::Free;
                break;

            case PhiDof::Port: {
                const std::size_t k = static_cast<std::size_t>(dofs.phi_port[u]);
                out.value[u] = dofs.port_is_fixed[k]
                                   ? dofs.port_value[k]
                                   : scale * solution.x[static_cast<std::size_t>(
                                                 dofs.port_index[k])];
                out.status[u] = PhiStatus::Terminal;
                break;
            }

            case PhiDof::Cut: {
                // Two values, not one: the plus side carries the port's unknown
                // and the grounded minus side carries zero. The plus side is
                // reported and the status says so, rather than picking one
                // silently -- Phi really is discontinuous across a cut, and that
                // discontinuity is what drives the loop.
                const std::size_t k = static_cast<std::size_t>(dofs.phi_port[u]);
                out.value[u] = dofs.port_is_fixed[k]
                                   ? dofs.port_value[k]
                                   : scale * solution.x[static_cast<std::size_t>(
                                                 dofs.port_index[k])];
                out.status[u] = PhiStatus::Cut;
                break;
            }

            case PhiDof::Fixed:
                // A prescribed value is already physical: it was never scaled by
                // the substitution, so it must not be scaled back.
                out.value[u] = dofs.phi_fixed_value[u];
                out.status[u] = PhiStatus::Fixed;
                break;

            case PhiDof::Absent:
                out.value[u] = Complex(0.0, 0.0);
                out.status[u] = PhiStatus::Absent;
                break;
        }
    }
    (void)bound;
    return out;
}

// ------------------------------------------------------------------ TextBuffer

TextBuffer::TextBuffer(std::FILE* out, std::size_t capacity) : out_(out) {
    if (out_ == nullptr) throw std::invalid_argument("TextBuffer: null file");
    // Honour a genuinely small capacity rather than clamping to something
    // comfortable: the flush boundary is only exercised when `used_` actually
    // reaches the end, and a 4 kB floor meant a test could never get there often
    // enough to catch an off-by-one in `make_room`.
    buffer_.resize(capacity < 16 ? 16 : capacity);
}

TextBuffer::~TextBuffer() {
    // Best effort: a throw from a destructor would be worse than a lost flush,
    // and callers that care call flush() themselves and can see the error.
    if (used_ > 0 && out_ != nullptr) {
        std::fwrite(buffer_.data(), 1, used_, out_);
        written_ += used_;
        used_ = 0;
    }
}

void TextBuffer::make_room(std::size_t n) {
    if (used_ + n > buffer_.size()) {
        flush();
        if (n > buffer_.size()) buffer_.resize(n);
    }
    // The postcondition, checked UNCONDITIONALLY rather than assumed. An
    // off-by-one in the test above lets a write run one byte past the end,
    // which is undefined behaviour that a vector does not catch and that
    // usually does nothing visible -- so a test cannot find it, but this can.
    // One comparison per call, against silent memory corruption.
    if (used_ + n > buffer_.size()) {
        throw std::logic_error("TextBuffer::make_room: asked for " + std::to_string(n) +
                               " bytes with " + std::to_string(used_) + " used of " +
                               std::to_string(buffer_.size()) +
                               ", and could not make room. The capacity check and the "
                               "writes disagree.");
    }
}

void TextBuffer::flush() {
    if (used_ == 0) return;
    const std::size_t done = std::fwrite(buffer_.data(), 1, used_, out_);
    if (done != used_) throw std::runtime_error("TextBuffer: short write");
    written_ += used_;
    used_ = 0;
}

void TextBuffer::put(char c) {
    make_room(1);
    buffer_[used_++] = c;
}

void TextBuffer::put(const char* text) {
    const std::size_t n = std::strlen(text);
    make_room(n);
    std::memcpy(buffer_.data() + used_, text, n);
    used_ += n;
}

void TextBuffer::put(int value) { put(static_cast<long long>(value)); }

void TextBuffer::put(long long value) {
    make_room(32);
    const auto r = std::to_chars(buffer_.data() + used_, buffer_.data() + buffer_.size(), value);
    if (r.ec != std::errc()) throw std::runtime_error("TextBuffer: integer formatting failed");
    used_ = static_cast<std::size_t>(r.ptr - buffer_.data());
}

void TextBuffer::put(double value) {
    make_room(64);
    // No format and no precision: that asks for the SHORTEST representation that
    // round-trips exactly. Shorter than a fixed %.17g and lossless, where a
    // fixed precision would be longer and still lossy.
    const auto r = std::to_chars(buffer_.data() + used_, buffer_.data() + buffer_.size(), value);
    if (r.ec != std::errc()) throw std::runtime_error("TextBuffer: double formatting failed");
    used_ = static_cast<std::size_t>(r.ptr - buffer_.data());
}

// -------------------------------------------------------------- write_potential

namespace {

char status_letter(PhiStatus s) {
    switch (s) {
        case PhiStatus::Free: return 'F';
        case PhiStatus::Terminal: return 'T';
        case PhiStatus::Cut: return 'C';
        case PhiStatus::Fixed: return 'P';
        case PhiStatus::Absent: return 'A';
    }
    return '?';
}

}  // namespace


FieldOutput compute_fields(const Mesh& mesh, const BoundProblem& bound, const DofMap& dofs,
                           const Solution& solution) {
    if (static_cast<int>(solution.x.size()) != dofs.num_total) {
        throw std::invalid_argument(
            "compute_fields: the solution has " + std::to_string(solution.x.size()) +
            " entries for a DOF map of " + std::to_string(dofs.num_total) +
            ". They must describe the same problem.");
    }

    const int nt = mesh.num_tets();
    const int nv = mesh.num_nodes();

    const int ne = mesh.num_edges();
    const int np = dofs.num_p2_nodes;

    FieldOutput out;
    out.num_vertices = nv;
    out.num_edges = ne;
    out.a_tet.assign(static_cast<std::size_t>(nt), Vec3C{});
    out.b_tet.assign(static_cast<std::size_t>(nt), Vec3C{});
    out.e_tet.assign(static_cast<std::size_t>(nt), Vec3C{});
    out.sigma_tet.assign(static_cast<std::size_t>(nt), 0.0);
    // Sized for every P2 node. The vertex part is filled by accumulation below;
    // the mid-edge part afterwards, from the two endpoints.
    out.a_node.assign(static_cast<std::size_t>(np), Vec3C{});
    out.b_node.assign(static_cast<std::size_t>(np), Vec3C{});
    out.h_node.assign(static_cast<std::size_t>(np), Vec3C{});
    out.e_node.assign(static_cast<std::size_t>(np), Vec3C{});
    out.vertex_weight.assign(static_cast<std::size_t>(nv), 0.0);
    out.on_material_interface.assign(static_cast<std::size_t>(np), 0);
    // Which body first claimed each vertex; a second, different one marks an
    // interface.
    std::vector<int> first_body(static_cast<std::size_t>(nv), -1);

    // Phi is P2 and already exact at its nodes, so it is copied, not rebuilt.
    const NodalPotential phi = potential_at_nodes(mesh, bound, dofs, solution);
    out.phi_node = phi.value;

    const Complex scale = solution.phi_scale();
    const Complex jw(0.0, solution.omega);

    // The four vertices of a tet in barycentric coordinates: evaluating A and E
    // there rather than at the centroid keeps the linear variation that a
    // centroid value would average away.
    static const std::array<std::array<double, 4>, 4> kCorner = {{{{1.0, 0.0, 0.0, 0.0}},
                                                                  {{0.0, 1.0, 0.0, 0.0}},
                                                                  {{0.0, 0.0, 1.0, 0.0}},
                                                                  {{0.0, 0.0, 0.0, 1.0}}}};
    static const std::array<double, 4> kCentroid = {0.25, 0.25, 0.25, 0.25};

    for (int t = 0; t < nt; ++t) {
        const std::size_t ut = static_cast<std::size_t>(t);
        const TetGeometry g = compute_tet_geometry(mesh, t);
        const TetDofs td = dofs.local_dofs(t, mesh, bound);

        // A's coefficients. An eliminated edge is zero whichever reason
        // eliminated it: n x A = 0 on a Dirichlet surface, or the tree-cotree
        // gauge. `whitney_edge_value_global` carries the orientation, so the
        // raw global unknown is used here and the sign is NOT applied twice --
        // assembly instead uses the local basis and multiplies by the same
        // sign at scatter time. A test cross-checks the two against each other.
        std::array<Complex, 6> a{};
        for (int L = 0; L < 6; ++L) {
            const DofEntry& d = td.edge[static_cast<std::size_t>(L)];
            if (d.index >= 0) a[static_cast<std::size_t>(L)] = solution.x[static_cast<std::size_t>(d.index)];
        }

        // Phi's coefficients, resolved for THIS tet: a node on an internal cut
        // reads the port's unknown or zero depending which side this tet is on,
        // which is why local_dofs is consulted per tet rather than per node.
        std::array<Complex, 10> ph{};
        for (int b = 0; b < 10; ++b) {
            const std::size_t ub = static_cast<std::size_t>(b);
            const DofEntry& d = td.phi[ub];
            if (d.index >= 0) {
                ph[ub] = scale * d.coeff * solution.x[static_cast<std::size_t>(d.index)];
            } else {
                // A prescribed value is already physical -- it never went
                // through the ScaledPhi substitution, so it must not come back
                // through it either.
                ph[ub] = td.phi_fixed[ub];
            }
        }

        // B = curl A is constant over the tet: curl W_e is.
        Vec3C b_const{};
        for (int L = 0; L < 6; ++L) {
            const Vec3 c = whitney_edge_curl_global(mesh, t, g, L);
            const Complex ae = a[static_cast<std::size_t>(L)];
            b_const[0] += ae * c.x;
            b_const[1] += ae * c.y;
            b_const[2] += ae * c.z;
        }
        out.b_tet[ut] = b_const;

        const int body = bound.body_of_tet[ut];
        const double mu = kMu0 * bound.bodies[static_cast<std::size_t>(body)].mu_r;
        const double nu = 1.0 / mu;
        out.sigma_tet[ut] = bound.bodies[static_cast<std::size_t>(body)].sigma;

        // A and E at a point of this tet.
        auto evaluate = [&](const std::array<double, 4>& L) {
            Vec3C av{}, ev{};
            for (int k = 0; k < 6; ++k) {
                const Vec3 w = whitney_edge_value_global(mesh, t, g, k, L);
                const Complex ae = a[static_cast<std::size_t>(k)];
                av[0] += ae * w.x;
                av[1] += ae * w.y;
                av[2] += ae * w.z;
            }
            // E = -j*w*A - grad(Phi)
            Vec3C gp{};
            for (int n = 0; n < 10; ++n) {
                const Vec3 gr = p2_nodal_gradient(g, n, L);
                const Complex pn = ph[static_cast<std::size_t>(n)];
                gp[0] += pn * gr.x;
                gp[1] += pn * gr.y;
                gp[2] += pn * gr.z;
            }
            for (int i = 0; i < 3; ++i) ev[static_cast<std::size_t>(i)] =
                -jw * av[static_cast<std::size_t>(i)] - gp[static_cast<std::size_t>(i)];
            return std::pair<Vec3C, Vec3C>(av, ev);
        };

        const std::pair<Vec3C, Vec3C> mid = evaluate(kCentroid);
        out.a_tet[ut] = mid.first;
        out.e_tet[ut] = mid.second;

        // Volume-weighted accumulation onto the vertices.
        const double vol = g.volume < 0.0 ? -g.volume : g.volume;
        for (int c = 0; c < 4; ++c) {
            const int v = mesh.tets[ut][static_cast<std::size_t>(c)];
            const std::size_t uv = static_cast<std::size_t>(v);
            const std::pair<Vec3C, Vec3C> at = evaluate(kCorner[static_cast<std::size_t>(c)]);
            out.vertex_weight[uv] += vol;
            for (int i = 0; i < 3; ++i) {
                const std::size_t ui = static_cast<std::size_t>(i);
                out.a_node[uv][ui] += vol * at.first[ui];
                out.e_node[uv][ui] += vol * at.second[ui];
                out.b_node[uv][ui] += vol * b_const[ui];
                out.h_node[uv][ui] += vol * nu * b_const[ui];
            }
            // A vertex whose incident tets span more than one body sits on a
            // material interface, where averaging across the two sides is
            // meaningless -- see FieldOutput::on_material_interface.
            if (first_body[uv] < 0) {
                first_body[uv] = body;
            } else if (first_body[uv] != body) {
                out.on_material_interface[uv] = 1;
            }
        }
    }

    for (int v = 0; v < nv; ++v) {
        const std::size_t uv = static_cast<std::size_t>(v);
        const double w = out.vertex_weight[uv];
        if (w <= 0.0) {
            ++out.num_orphan_vertices;
            continue;
        }
        for (int i = 0; i < 3; ++i) {
            const std::size_t ui = static_cast<std::size_t>(i);
            out.a_node[uv][ui] /= w;
            out.b_node[uv][ui] /= w;
            out.h_node[uv][ui] /= w;
            out.e_node[uv][ui] /= w;
        }
    }

    // The mid-edge nodes: the mean of the two endpoint values.
    //
    // `A`, `B`, `H` and `E` have no mid-edge degree of freedom, so there is
    // nothing to read there and the endpoint mean is the standard
    // reconstruction. This is exactly what must NOT be done to `Phi`, whose
    // mid-edge value is a genuine P2 unknown -- see `phi_node`, which is
    // copied rather than averaged. The asymmetry is the point.
    //
    // An edge with an orphan endpoint inherits that endpoint's zero, which is
    // the same silence a vertex with no incident tet already carries.
    for (int e = 0; e < ne; ++e) {
        const std::pair<int, int>& ends = mesh.edges[static_cast<std::size_t>(e)];
        const std::size_t a0 = static_cast<std::size_t>(ends.first);
        const std::size_t a1 = static_cast<std::size_t>(ends.second);
        const std::size_t m = static_cast<std::size_t>(dofs.edge_p2(e));
        for (int i = 0; i < 3; ++i) {
            const std::size_t ui = static_cast<std::size_t>(i);
            out.a_node[m][ui] = 0.5 * (out.a_node[a0][ui] + out.a_node[a1][ui]);
            out.b_node[m][ui] = 0.5 * (out.b_node[a0][ui] + out.b_node[a1][ui]);
            out.h_node[m][ui] = 0.5 * (out.h_node[a0][ui] + out.h_node[a1][ui]);
            out.e_node[m][ui] = 0.5 * (out.e_node[a0][ui] + out.e_node[a1][ui]);
        }
        // An edge whose endpoint is on an interface straddles it too.
        out.on_material_interface[m] =
            (out.on_material_interface[a0] != 0 || out.on_material_interface[a1] != 0) ? 1 : 0;
    }

    for (int i = 0; i < np; ++i) {
        if (out.on_material_interface[static_cast<std::size_t>(i)] != 0) ++out.num_interface_nodes;
    }

    return out;
}

WriteStats write_potential(const std::string& path, const NodalPotential& potential,
                           const Solution& solution) {
    const auto started = std::chrono::steady_clock::now();

    std::FILE* f = std::fopen(path.c_str(), "wb");
    if (f == nullptr) {
        throw std::runtime_error("write_potential: could not open '" + path + "' for writing");
    }

    WriteStats stats;
    stats.nodes = potential.size();
    {
        TextBuffer buf(f);

        buf.put("# potential.out -- nodal electric scalar potential Phi, in volts\n");
        buf.put("# A-Phi solver ");
        buf.put(kVersion);
        buf.put("\n# conditioning ");
        buf.put(conditioning_keyword(solution.conditioning));
        buf.put("   omega ");
        buf.put(solution.omega);
        buf.put(" rad/s   frequency ");
        buf.put(solution.omega / 6.283185307179586476925286766559);
        buf.put(" Hz\n");
        buf.put("# nodes ");
        buf.put(potential.size());
        buf.put("  (vertices ");
        buf.put(potential.num_vertices);
        buf.put(", edge midpoints ");
        buf.put(potential.num_edges);
        buf.put(")\n");
        buf.put("# positions in metres. Mid-edge values are EXACT P2 unknowns, not interpolated.\n");
        buf.put("# status  F free   T terminal   C cut (value is the PLUS side; Phi jumps here)"
                "   P prescribed   A absent\n");
        buf.put("# index x y z Re(Phi) Im(Phi) status\n");

        const int n = potential.size();
        for (int i = 0; i < n; ++i) {
            const std::size_t u = static_cast<std::size_t>(i);
            buf.put(i);
            buf.put(' ');
            buf.put(potential.position[u].x);
            buf.put(' ');
            buf.put(potential.position[u].y);
            buf.put(' ');
            buf.put(potential.position[u].z);
            buf.put(' ');
            buf.put(potential.value[u].real());
            buf.put(' ');
            buf.put(potential.value[u].imag());
            buf.put(' ');
            buf.put(status_letter(potential.status[u]));
            buf.put('\n');
        }
        buf.flush();
        stats.bytes = buf.bytes_written();
    }
    std::fclose(f);

    stats.milliseconds =
        std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - started)
            .count();
    return stats;
}



WriteStats write_solution(const std::string& path, const NodalPotential& potential,
                          const FieldOutput& fields, const Solution& solution,
                          FieldSet which) {
    const auto started = std::chrono::steady_clock::now();
    if (fields.num_nodes() != potential.size()) {
        throw std::invalid_argument(
            "write_solution: the fields cover " + std::to_string(fields.num_nodes()) +
            " nodes and the potential " + std::to_string(potential.size()) +
            ". They must come from the same solve.");
    }

    std::FILE* f = std::fopen(path.c_str(), "wb");
    if (f == nullptr) {
        throw std::runtime_error("write_solution: could not open '" + path + "' for writing");
    }

    WriteStats stats;
    stats.nodes = potential.size();
    {
        TextBuffer buf(f);
        const int n = potential.size();

        const bool all = which == FieldSet::All;
        const bool w_phi = all || which == FieldSet::Potential;
        const bool w_a = all || which == FieldSet::A;
        const bool w_b = all || which == FieldSet::B;
        const bool w_h = all || which == FieldSet::H;
        const bool w_e = all || which == FieldSet::E;

        buf.put("# ");
        buf.put(field_set_name(which));
        buf.put(all ? " -- Phi, A, B, H and E at every P2 node\n"
                    : " at every P2 node\n");
        buf.put("# A-Phi solver ");
        buf.put(kVersion);
        buf.put("\n# conditioning ");
        buf.put(conditioning_keyword(solution.conditioning));
        buf.put("   omega ");
        buf.put(solution.omega);
        buf.put(" rad/s   frequency ");
        buf.put(solution.omega / 6.283185307179586476925286766559);
        buf.put(" Hz\n");
        buf.put("# nodes ");
        buf.put(n);
        buf.put("  (vertices ");
        buf.put(fields.num_vertices);
        buf.put(", edge midpoints ");
        buf.put(fields.num_edges);
        buf.put(")\n");
        buf.put("#\n");
        buf.put("# Units: positions m, Phi V, A Wb/m, B T, H A/m, E V/m. Every field is a\n");
        buf.put("# complex phasor under exp(+j w t), written as its real then imaginary part.\n");
        buf.put("#\n");
        buf.put("# Phi is EXACT at every node, including the mid-edge ones: it lives in a P2\n");
        buf.put("# space and a midpoint carries a genuine solved unknown. A, B, H and E have\n");
        buf.put("# no mid-edge freedom, so at a midpoint they are the mean of the two\n");
        buf.put("# endpoints. docs/FIELD_POSTPROCESSING.md Sec. 8.\n");
        buf.put("#\n");
        buf.put("# WHAT MAY BE READ FROM THIS FILE. The discrete system is exactly gauge\n");
        buf.put("# invariant, and the tree-cotree gauge picks one representative, so Phi and A\n");
        buf.put("# are gauge DEPENDENT -- a different spanning tree shifts Phi by -j*w*psi,\n");
        buf.put("# almost purely imaginary. B, H and E are not. Sec. 9.\n");
        buf.put("#\n");
        buf.put("# iface = 1 means the node's incident tets span more than one body. The nodal\n");
        buf.put("# average straddles a material interface there and E and H are MEANINGLESS:\n");
        buf.put("# E's normal component genuinely jumps across it. Filter on this column.\n");
        buf.put("# Sec. 12. The per-cell arrays in the .vtk have no interface to straddle.\n");
        buf.put("#\n");
        // The header names exactly the columns the rows below carry. A test
        // counts the two against each other: in a file this wide, adding a
        // column to one and not the other still parses, silently mislabelled
        // from that point on.
        buf.put("# index x y z");
        if (w_phi) buf.put("  Re(Phi) Im(Phi)");
        if (w_a) buf.put("  Re(Ax) Im(Ax) Re(Ay) Im(Ay) Re(Az) Im(Az)");
        if (w_b) buf.put("  Re(Bx) Im(Bx) Re(By) Im(By) Re(Bz) Im(Bz)");
        if (w_h) buf.put("  Re(Hx) Im(Hx) Re(Hy) Im(Hy) Re(Hz) Im(Hz)");
        if (w_e) buf.put("  Re(Ex) Im(Ex) Re(Ey) Im(Ey) Re(Ez) Im(Ez)");
        buf.put("  iface\n");

        auto put_vec = [&buf](const Vec3C& v) {
            for (int k = 0; k < 3; ++k) {
                const std::complex<double>& c = v[static_cast<std::size_t>(k)];
                buf.put(' ');
                buf.put(c.real());
                buf.put(' ');
                buf.put(c.imag());
            }
        };

        for (int i = 0; i < n; ++i) {
            const std::size_t u = static_cast<std::size_t>(i);
            buf.put(i);
            const Vec3& p = potential.position[u];
            buf.put(' ');
            buf.put(p.x);
            buf.put(' ');
            buf.put(p.y);
            buf.put(' ');
            buf.put(p.z);
            if (w_phi) {
                buf.put(' ');
                buf.put(fields.phi_node[u].real());
                buf.put(' ');
                buf.put(fields.phi_node[u].imag());
            }
            if (w_a) put_vec(fields.a_node[u]);
            if (w_b) put_vec(fields.b_node[u]);
            if (w_h) put_vec(fields.h_node[u]);
            if (w_e) put_vec(fields.e_node[u]);
            buf.put(' ');
            buf.put(static_cast<int>(fields.on_material_interface[u]));
            buf.put('\n');
        }

        buf.flush();
        stats.bytes = buf.bytes_written();
    }
    std::fclose(f);
    stats.milliseconds =
        std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - started)
            .count();
    return stats;
}

namespace {

/// One VTK `VECTORS` array: the real or the imaginary part of a complex field.
void put_vectors(TextBuffer& buf, const char* name, const std::vector<Vec3C>& v, int n,
                 bool real_part) {
    buf.put("\nVECTORS ");
    buf.put(name);
    buf.put(" double\n");
    for (int i = 0; i < n; ++i) {
        const Vec3C& x = v[static_cast<std::size_t>(i)];
        for (int k = 0; k < 3; ++k) {
            if (k != 0) buf.put(' ');
            const std::complex<double>& c = x[static_cast<std::size_t>(k)];
            buf.put(real_part ? c.real() : c.imag());
        }
        buf.put('\n');
    }
}

void put_magnitude(TextBuffer& buf, const char* name, const std::vector<Vec3C>& v, int n) {
    buf.put("\nSCALARS ");
    buf.put(name);
    buf.put(" double 1\nLOOKUP_TABLE default\n");
    for (int i = 0; i < n; ++i) {
        const Vec3C& x = v[static_cast<std::size_t>(i)];
        double s = 0.0;
        for (int k = 0; k < 3; ++k) s += std::norm(x[static_cast<std::size_t>(k)]);
        buf.put(std::sqrt(s));
        buf.put('\n');
    }
}

}  // namespace

WriteStats write_vtk(const std::string& path, const Mesh& mesh, const NodalPotential& potential,
                     const Solution& solution, const FieldOutput* fields, FieldSet which) {
    const auto started = std::chrono::steady_clock::now();
    if (potential.size() != mesh.num_nodes() + mesh.num_edges()) {
        throw std::invalid_argument(
            "write_vtk: the potential has " + std::to_string(potential.size()) +
            " nodes for a mesh of " + std::to_string(mesh.num_nodes()) + " vertices and " +
            std::to_string(mesh.num_edges()) + " edges");
    }

    std::FILE* f = std::fopen(path.c_str(), "wb");
    if (f == nullptr) {
        throw std::runtime_error("write_vtk: could not open '" + path + "' for writing");
    }

    WriteStats stats;
    stats.nodes = potential.size();
    {
        TextBuffer buf(f);
        const int np = potential.size();
        const int nc = mesh.num_tets();

        const bool all = which == FieldSet::All;
        const bool w_phi = all || which == FieldSet::Potential;
        const bool w_a = all || which == FieldSet::A;
        const bool w_b = all || which == FieldSet::B;
        const bool w_h = all || which == FieldSet::H;
        const bool w_e = all || which == FieldSet::E;

        buf.put("# vtk DataFile Version 3.0\n");
        buf.put("A-Phi ");
        buf.put(field_set_name(which));
        buf.put(", conditioning ");
        buf.put(conditioning_keyword(solution.conditioning));
        buf.put(", omega ");
        buf.put(solution.omega);
        buf.put(" rad/s\n");
        buf.put("ASCII\nDATASET UNSTRUCTURED_GRID\n");

        buf.put("POINTS ");
        buf.put(np);
        buf.put(" double\n");
        for (int i = 0; i < np; ++i) {
            const Vec3& p = potential.position[static_cast<std::size_t>(i)];
            buf.put(p.x);
            buf.put(' ');
            buf.put(p.y);
            buf.put(' ');
            buf.put(p.z);
            buf.put('\n');
        }

        // Quadratic tets: four corners then six mid-edge nodes, in VTK's order
        // rather than ours -- see kVtkQuadraticTetEdgeOrder.
        buf.put("\nCELLS ");
        buf.put(nc);
        buf.put(' ');
        buf.put(nc * 11);
        buf.put('\n');
        for (int t = 0; t < nc; ++t) {
            const std::size_t ut = static_cast<std::size_t>(t);
            buf.put("10");
            for (int v = 0; v < 4; ++v) {
                buf.put(' ');
                buf.put(mesh.tets[ut][static_cast<std::size_t>(v)]);
            }
            for (int slot = 0; slot < 6; ++slot) {
                const int local_edge = kVtkQuadraticTetEdgeOrder[static_cast<std::size_t>(slot)];
                const int e = mesh.tet_edges[ut][static_cast<std::size_t>(local_edge)];
                buf.put(' ');
                buf.put(mesh.num_nodes() + e);
            }
            buf.put('\n');
        }

        buf.put("\nCELL_TYPES ");
        buf.put(nc);
        buf.put('\n');
        for (int t = 0; t < nc; ++t) buf.put("24\n");  // VTK_QUADRATIC_TETRA

        buf.put("\nCELL_DATA ");
        buf.put(nc);
        buf.put("\nSCALARS body_tag int 1\nLOOKUP_TABLE default\n");
        for (int t = 0; t < nc; ++t) {
            buf.put(t < static_cast<int>(mesh.tet_tags.size())
                        ? mesh.tet_tags[static_cast<std::size_t>(t)]
                        : -1);
            buf.put('\n');
        }

        if (fields != nullptr) {
            // Per-cell quantities, written as they are computed rather than
            // smoothed. `B` is genuinely constant per tet, so `B_cell` is the
            // exact value and the nodal `B` below is a post-processing average
            // of it -- carrying both is what lets the smoothing be seen.
            buf.put("\nSCALARS sigma double 1\nLOOKUP_TABLE default\n");
            for (int t = 0; t < nc; ++t) {
                buf.put(fields->sigma_tet[static_cast<std::size_t>(t)]);
                buf.put('\n');
            }
            if (w_b || w_h) {
                put_vectors(buf, "B_cell_real", fields->b_tet, nc, true);
                put_vectors(buf, "B_cell_imag", fields->b_tet, nc, false);
            }
            // E per cell as well as per node. In an insulator `J` is zero, so
            // without this there is no averaging-free E anywhere outside the
            // conductors -- which is exactly where the nodal one is least
            // trustworthy, since the interface runs along that boundary.
            if (w_e) {
                put_vectors(buf, "E_cell_real", fields->e_tet, nc, true);
                put_vectors(buf, "E_cell_imag", fields->e_tet, nc, false);

                // J = sigma E, formed per tet because that is where sigma is
                // single-valued: at a node on a conductor/insulator interface it
                // is not. Zero throughout an insulator, which is correct and not
                // a gap. It travels with E because it IS E, scaled.
                std::vector<Vec3C> j(static_cast<std::size_t>(nc), Vec3C{});
                for (int t = 0; t < nc; ++t) {
                    const std::size_t ut = static_cast<std::size_t>(t);
                    const double s = fields->sigma_tet[ut];
                    for (int k = 0; k < 3; ++k) {
                        j[ut][static_cast<std::size_t>(k)] =
                            s * fields->e_tet[ut][static_cast<std::size_t>(k)];
                    }
                }
                put_vectors(buf, "J_real", j, nc, true);
                put_vectors(buf, "J_imag", j, nc, false);
            }
        }

        buf.put("\nPOINT_DATA ");
        buf.put(np);
        buf.put('\n');

        if (w_phi) {
            buf.put("SCALARS phi_real double 1\nLOOKUP_TABLE default\n");
            for (int i = 0; i < np; ++i) {
                buf.put(potential.value[static_cast<std::size_t>(i)].real());
                buf.put('\n');
            }
            buf.put("\nSCALARS phi_imag double 1\nLOOKUP_TABLE default\n");
            for (int i = 0; i < np; ++i) {
                buf.put(potential.value[static_cast<std::size_t>(i)].imag());
                buf.put('\n');
            }
            buf.put("\nSCALARS phi_magnitude double 1\nLOOKUP_TABLE default\n");
            for (int i = 0; i < np; ++i) {
                buf.put(std::abs(potential.value[static_cast<std::size_t>(i)]));
                buf.put('\n');
            }
            // Where Phi does not live, the value is zero because there is
            // nothing to report -- not because the solve found zero volts
            // there. Without this field a viewer cannot tell those apart.
            buf.put("\nSCALARS phi_present int 1\nLOOKUP_TABLE default\n");
            for (int i = 0; i < np; ++i) {
                buf.put(potential.status[static_cast<std::size_t>(i)] == PhiStatus::Absent ? 0 : 1);
                buf.put('\n');
            }
        }

        if (fields != nullptr) {
            if (fields->num_nodes() != np) {
                throw std::invalid_argument(
                    "write_vtk: the fields cover " + std::to_string(fields->num_nodes()) +
                    " nodes and the potential " + std::to_string(np) +
                    ". They must come from the same solve.");
            }
            // A complex phasor at every P2 node, real and imaginary parts as
            // separate VECTORS because the legacy format has no complex type.
            // ParaView's own Calculator can form any combination of them.
            // The `*_magnitude` scalars are the phasor amplitude
            // sqrt(|Xx|^2 + |Xy|^2 + |Xz|^2), which is NOT the length of either
            // the real or the imaginary vector: where the field is elliptically
            // polarised neither of those is the physical peak, and this is.
            if (w_a) {
                put_vectors(buf, "A_real", fields->a_node, np, true);
                put_vectors(buf, "A_imag", fields->a_node, np, false);
                put_magnitude(buf, "A_magnitude", fields->a_node, np);
            }
            if (w_b) {
                put_vectors(buf, "B_real", fields->b_node, np, true);
                put_vectors(buf, "B_imag", fields->b_node, np, false);
                put_magnitude(buf, "B_magnitude", fields->b_node, np);
            }
            if (w_h) {
                put_vectors(buf, "H_real", fields->h_node, np, true);
                put_vectors(buf, "H_imag", fields->h_node, np, false);
                put_magnitude(buf, "H_magnitude", fields->h_node, np);
            }
            if (w_e) {
                put_vectors(buf, "E_real", fields->e_node, np, true);
                put_vectors(buf, "E_imag", fields->e_node, np, false);
                put_magnitude(buf, "E_magnitude", fields->e_node, np);
            }

            // 1 where the nodal average straddles a material interface and the
            // nodal fields are therefore meaningless. Threshold this to 0
            // before reading E or H, or use the per-cell arrays instead.
            buf.put("\nSCALARS material_interface int 1\nLOOKUP_TABLE default\n");
            for (int i = 0; i < np; ++i) {
                buf.put(static_cast<int>(
                    fields->on_material_interface[static_cast<std::size_t>(i)]));
                buf.put('\n');
            }
        }

        buf.flush();
        stats.bytes = buf.bytes_written();
    }
    std::fclose(f);
    stats.milliseconds =
        std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - started)
            .count();
    return stats;
}

}  // namespace aphi_solver
