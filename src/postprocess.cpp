#include "aphi_solver/postprocess.hpp"

#include <charconv>
#include <cmath>
#include <chrono>
#include <cstring>
#include <stdexcept>
#include <string>

#include "aphi_solver/version.hpp"

namespace aphi_solver {

using Complex = std::complex<double>;

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


WriteStats write_vtk(const std::string& path, const Mesh& mesh, const NodalPotential& potential,
                     const Solution& solution) {
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

        buf.put("# vtk DataFile Version 3.0\n");
        buf.put("A-Phi potential, conditioning ");
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

        buf.put("\nPOINT_DATA ");
        buf.put(np);
        buf.put('\n');

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
        // Where Phi does not live, the value is zero because there is nothing to
        // report -- not because the solve found zero volts there. Without this
        // field a viewer cannot tell those apart.
        buf.put("\nSCALARS phi_present int 1\nLOOKUP_TABLE default\n");
        for (int i = 0; i < np; ++i) {
            buf.put(potential.status[static_cast<std::size_t>(i)] == PhiStatus::Absent ? 0 : 1);
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

}  // namespace aphi_solver
