#include "aphi_solver/gauge_free_potential.hpp"

#include <algorithm>
#include <chrono>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>

#include "aphi_solver/basis_functions.hpp"
#include "aphi_solver/factorization.hpp"

namespace aphi_solver {
namespace {

using Complex = std::complex<double>;
using Clock = std::chrono::steady_clock;

double ms_since(Clock::time_point t0) {
    return std::chrono::duration<double, std::milli>(Clock::now() - t0).count();
}

constexpr double kEps0 = 8.8541878128e-12;

/// Where each mesh vertex's Phi' lands.
///
/// Deliberately the SAME rules the main solve uses for Phi, so the recovered
/// potential can be compared with the one written beside it without having to
/// reason about two different conventions:
///
///   - a vertex on a current-driven port shares that port's ONE unknown, which
///     is `PhiDof::Port` and is what makes a terminal a terminal;
///   - a vertex on a voltage port is prescribed and leaves the system;
///   - a vertex outside Phi's support is absent;
///   - everything else is free.
///
/// `DofMap` numbers P2 nodes with the mesh vertices first, so `phi_state[v]`
/// for `v < num_nodes()` is this vertex's state.
struct VertexMap {
    std::vector<int> column;              ///< -1 if eliminated
    std::vector<Complex> prescribed;      ///< value when eliminated, else 0
    std::vector<int> port_column;         ///< per port, -1 if it has no unknown
    int unknowns = 0;
};

VertexMap build_vertex_map(const Mesh& mesh, const BoundProblem& bound, const DofMap& dofs) {
    const int nv = mesh.num_nodes();
    VertexMap m;
    m.column.assign(static_cast<std::size_t>(nv), -1);
    m.prescribed.assign(static_cast<std::size_t>(nv), Complex{});
    m.port_column.assign(bound.ports.size(), -1);

    // Ports first, so every vertex of one terminal gets the same column.
    for (std::size_t p = 0; p < bound.ports.size(); ++p) {
        if (dofs.port_is_fixed[p]) continue;
        m.port_column[p] = m.unknowns++;
    }
    for (int v = 0; v < nv; ++v) {
        const std::size_t uv = static_cast<std::size_t>(v);
        switch (dofs.phi_state[uv]) {
            case PhiDof::Port: {
                const int p = dofs.phi_port[uv];
                if (p < 0) break;
                const std::size_t up = static_cast<std::size_t>(p);
                if (dofs.port_is_fixed[up]) {
                    m.prescribed[uv] = dofs.port_value[up];
                } else {
                    m.column[uv] = m.port_column[up];
                }
                break;
            }
            case PhiDof::Free:
                m.column[uv] = m.unknowns++;
                break;
            case PhiDof::Cut:
            case PhiDof::Fixed:
            case PhiDof::Absent:
                // A cut's Phi jumps, so a single-valued P1 recovery cannot
                // represent it and the vertex is left prescribed at zero
                // rather than quietly averaged across the jump. Absent and
                // Fixed are eliminated as in the main solve.
                break;
        }
    }
    return m;
}

}  // namespace

GaugeFreePotential recover_gauge_free_potential(const Mesh& mesh, const BoundProblem& bound,
                                                const DofMap& dofs, const FieldOutput& fields,
                                                double omega) {
    const Clock::time_point t0 = Clock::now();
    GaugeFreePotential out;
    const int nv = mesh.num_nodes();
    const int nt = mesh.num_tets();
    if (static_cast<int>(fields.e_tet.size()) != nt) {
        throw std::invalid_argument(
            "recover_gauge_free_potential: fields.e_tet must have one entry per tet");
    }

    const VertexMap vm = build_vertex_map(mesh, bound, dofs);
    out.unknowns = vm.unknowns;
    out.phi_vertex.assign(static_cast<std::size_t>(nv), Complex{});
    out.port_potential.assign(bound.ports.size(), Complex{});
    if (vm.unknowns == 0) {
        out.converged = true;
        out.assemble_ms = ms_since(t0);
        return out;
    }

    // --- pattern: one row per unknown, columns from tet adjacency ----------
    std::vector<std::map<int, Complex>> rows(static_cast<std::size_t>(vm.unknowns));
    std::vector<Complex> rhs(static_cast<std::size_t>(vm.unknowns), Complex{});

    for (int t = 0; t < nt; ++t) {
        const std::size_t ut = static_cast<std::size_t>(t);
        const int body = bound.body_of_tet[ut];
        if (body < 0) continue;
        const BoundBody& mat = bound.bodies[static_cast<std::size_t>(body)];
        // beta = sigma + j w eps, the same weight the main assembly uses for
        // the Phi block. Any constant multiple of beta would do -- it appears
        // on both sides -- but matching the solver keeps the two comparable.
        const Complex beta(mat.sigma, omega * kEps0 * mat.eps_r);
        if (beta == Complex{}) continue;

        const TetGeometry g = compute_tet_geometry(mesh, t);
        const std::array<int, 4>& verts = mesh.tets[ut];
        const Vec3C& E = fields.e_tet[ut];

        for (int i = 0; i < 4; ++i) {
            const int ri = vm.column[static_cast<std::size_t>(verts[i])];
            if (ri < 0) continue;
            const Vec3& gi = g.grad_L[static_cast<std::size_t>(i)];
            // RHS: -<beta E, grad lam>. E is complex, grad lam is real.
            const Complex edotg = E[0] * gi.x + E[1] * gi.y + E[2] * gi.z;
            rhs[static_cast<std::size_t>(ri)] -= beta * g.volume * edotg;

            for (int j = 0; j < 4; ++j) {
                const Vec3& gj = g.grad_L[static_cast<std::size_t>(j)];
                const Complex k = beta * g.volume * (gi.x * gj.x + gi.y * gj.y + gi.z * gj.z);
                const int cj = vm.column[static_cast<std::size_t>(verts[j])];
                if (cj < 0) {
                    // Eliminated: its known value moves to the right-hand side.
                    const Complex val = vm.prescribed[static_cast<std::size_t>(verts[j])];
                    if (val != Complex{}) rhs[static_cast<std::size_t>(ri)] -= k * val;
                    continue;
                }
                if (cj < ri) continue;  // upper triangle only
                rows[static_cast<std::size_t>(ri)][cj] += k;
            }
        }
    }

    std::vector<int> row_ptr(static_cast<std::size_t>(vm.unknowns) + 1, 0);
    std::vector<int> col_index;
    for (int r = 0; r < vm.unknowns; ++r) {
        row_ptr[static_cast<std::size_t>(r) + 1] =
            row_ptr[static_cast<std::size_t>(r)] + static_cast<int>(rows[static_cast<std::size_t>(r)].size());
        for (const auto& [c, v] : rows[static_cast<std::size_t>(r)]) col_index.push_back(c);
    }
    SparseSymmetricZ a = SparseSymmetricZ::from_pattern(vm.unknowns, row_ptr, col_index);
    {
        std::vector<Complex>& vals = a.mutable_values();
        std::size_t k = 0;
        for (int r = 0; r < vm.unknowns; ++r) {
            for (const auto& [c, v] : rows[static_cast<std::size_t>(r)]) vals[k++] = v;
        }
    }
    out.assemble_ms = ms_since(t0);

    const Clock::time_point t1 = Clock::now();
    std::vector<Complex> x;
    SolveReport rep;
    SolveOptions opts;
    opts.equilibration_iterations = 10;
    out.converged = solve_symmetric(a, rhs, x, rep, opts);
    out.solve_ms = ms_since(t1);
    out.residual = rep.residual;
    if (!out.converged) return out;

    for (int v = 0; v < nv; ++v) {
        const std::size_t uv = static_cast<std::size_t>(v);
        const int c = vm.column[uv];
        out.phi_vertex[uv] = (c >= 0) ? x[static_cast<std::size_t>(c)] : vm.prescribed[uv];
    }
    for (std::size_t p = 0; p < bound.ports.size(); ++p) {
        out.port_potential[p] = dofs.port_is_fixed[p]
                                    ? dofs.port_value[p]
                                    : x[static_cast<std::size_t>(vm.port_column[p])];
    }
    return out;
}

WriteStats write_gauge_free_potential(const std::string& path, const Mesh& mesh,
                                      const BoundProblem& bound, const GaugeFreePotential& gfp,
                                      const RunInfo* run) {
    const Clock::time_point t0 = Clock::now();
    std::ostringstream buf;
    buf << "# potential_gaugefree.out -- electric scalar potential Phi, in volts,\n"
        << "# RECOVERED FROM THE FIELDS rather than read out of the gauged system.\n"
        << "#\n"
        << "# Solved from its own boundary value problem,\n"
        << "#\n"
        << "#     <beta grad(Phi), grad(lam)> = -<beta E, grad(lam)>,  beta = sigma + j w eps\n"
        << "#\n"
        << "# which is Stysch Sec. 3.1/6.3 (GAUGE_CHOICE.md Sec. 11.8). Because its only\n"
        << "# input is E, and E does not depend on the gauge, neither does this Phi.\n"
        << "#\n"
        << "# WHY IT EXISTS. Phi in potential.out comes from the tree-cotree gauged\n"
        << "# system. Where every terminal sits on the outer boundary that is fine: n x A\n"
        << "# = 0 pins psi there, and 07_GaugeInvariance measures the terminal quantities\n"
        << "# invariant to 1e-12 across three spanning trees. Where a terminal is INSIDE\n"
        << "# the domain nothing pins it, and 08_MixedPort_Interior measures the gauged\n"
        << "# terminal potential moving 3.8 % with the tree against 1.1e-14 for this one.\n"
        << "#\n"
        << "# WHAT IT DOES NOT FIX, AND THIS MATTERS. E, B, H and J are INPUTS here and\n"
        << "# are returned unchanged. At a mixed-material interior port they move 12 %,\n"
        << "# 91 %, 91 % and 77 % with the tree and remain wrong, and an energy-based\n"
        << "# inductance moves 42 % (tools/gauge_spike/RESULTS_SPIKEC.md). The three\n"
        << "# solutions differ by a SOLENOIDAL field, so div() annihilates the\n"
        << "# difference: this recovers the irrotational part exactly, which is where R\n"
        << "# lives, and is blind to the solenoidal part, which carries the flux and the\n"
        << "# inductance.\n"
        << "#\n"
        << "# So: read R from this file. Do NOT read L, B or H from the same solve and\n"
        << "# assume they are equally sound.\n"
        << "#\n"
        << "# One row per mesh VERTEX (P1), not per P2 node: the source div(beta E) is\n"
        << "# piecewise constant on a tet, so P2 would double the system and buy nothing.\n"
        << "#\n"
        << "# columns: index  x  y  z  Re(Phi)  Im(Phi)      positions in metres\n";
    if (run != nullptr) {
        buf << "#\n"
            << "# --- run -------------------------------------------------------------\n"
            << "# input        " << run->input_file << "\n"
            << "# mesh         " << run->mesh_file << "   " << run->num_tets << " tets\n";
    }
    buf << "# unknowns     " << gfp.unknowns << "   residual " << gfp.residual
        << (gfp.converged ? "" : "   *** DID NOT CONVERGE") << "\n";
    for (std::size_t p = 0; p < bound.ports.size(); ++p) {
        buf << "# port         " << bound.ports[p].name << "   Phi = "
            << gfp.port_potential[p].real() << " "
            << (gfp.port_potential[p].imag() < 0 ? "-" : "+") << " "
            << std::abs(gfp.port_potential[p].imag()) << "j V\n";
    }
    buf << "#\n";
    buf.setf(std::ios::scientific);
    buf.precision(16);
    for (int v = 0; v < mesh.num_nodes(); ++v) {
        const std::size_t uv = static_cast<std::size_t>(v);
        const Vec3& x = mesh.nodes[uv];
        buf << v << " " << x.x << " " << x.y << " " << x.z << " "
            << gfp.phi_vertex[uv].real() << " " << gfp.phi_vertex[uv].imag() << "\n";
    }
    const std::string text = buf.str();
    std::ofstream file(path, std::ios::binary);
    file << text;
    WriteStats st;
    st.nodes = mesh.num_nodes();
    st.bytes = text.size();
    st.milliseconds = ms_since(t0);
    return st;
}

}  // namespace aphi_solver
