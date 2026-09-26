// Command-line driver.
//
// **Reads, binds and assembles.** It parses an input file, opens the mesh it
// names, resolves every name against it, checks everything each stage can
// check, builds the DOF map and the sparsity pattern, and assembles the
// system at every frequency the file asks for.
//
// It does **not** solve: there is no linear solver yet. The closing lines say
// so, so the output never implies more has been done than has.

#include <chrono>
#include <memory>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "aphi_solver/assembly.hpp"
#include "aphi_solver/dof_map.hpp"
#include "aphi_solver/gmsh_reader.hpp"
#include "aphi_solver/input_file.hpp"
#include "aphi_solver/problem_binding.hpp"
#include "aphi_solver/version.hpp"

namespace {

using namespace aphi_solver;

/// A node belonging to no tet becomes its own reference group, so the root
/// count exceeds the number of connected pieces and the gauge looks
/// incomplete when it is not. Such a node has no edges and no Phi, so it
/// contributes nothing to the matrix -- but saying nothing about it leaves a
/// reader to wonder, which is how an hour got spent on exactly this.
std::string orphan_note(const Mesh& mesh, const BoundProblem& b) {
    std::vector<char> used(static_cast<std::size_t>(mesh.num_nodes()), 0);
    for (const TetVerts& t : mesh.tets) {
        for (int v : t) used[static_cast<std::size_t>(v)] = 1;
    }
    int orphans = 0;
    for (int i = 0; i < mesh.num_nodes(); ++i) {
        if (!used[static_cast<std::size_t>(i)]) ++orphans;
    }
    if (orphans == 0) return {};
    std::ostringstream out;
    out << "              " << orphans << " mesh node" << (orphans == 1 ? " belongs" : "s belong")
        << " to no tet, so " << (orphans == 1 ? "it counts as an extra root" : "they count as extra roots")
        << ";\n"
        << "              " << (orphans == 1 ? "it has" : "they have")
        << " no edges and no Phi, so the gauge over the real mesh is complete\n";
    (void)b;
    return out.str();
}

const char* unit_name(LengthUnit u) {
    switch (u) {
        case LengthUnit::Metre: return "m";
        case LengthUnit::Millimetre: return "mm";
        case LengthUnit::Micrometre: return "um";
        case LengthUnit::Nanometre: return "nm";
    }
    return "?";
}

const char* port_type_name(PortType t) {
    switch (t) {
        case PortType::BoundaryCurrent: return "boundary_current";
        case PortType::BoundaryVoltage: return "boundary_voltage";
        case PortType::InternalCurrent: return "internal_current";
        case PortType::InternalVoltage: return "internal_voltage";
    }
    return "?";
}

// Turns -0.0 into 0.0. Negative zero is mathematically fine but reads as a
// typo in a direction vector.
double tidy(double v) { return v + 0.0 == 0.0 ? 0.0 : v; }

std::string join(const std::vector<std::string>& v) {
    std::string out;
    for (std::size_t i = 0; i < v.size(); ++i) {
        if (i) out += " ";
        out += v[i];
    }
    return out;
}

std::string format_amplitude(const Port& p) {
    std::ostringstream s;
    s << (is_current_driven(p.type) ? "I = " : "V = ") << p.amplitude
      << (is_current_driven(p.type) ? " A" : " V");
    if (p.phase_deg != 0.0) s << " at " << p.phase_deg << " deg";
    return s.str();
}

// A long sweep is elided rather than dumped: the point is to confirm the
// endpoints and the count, which is what a mistyped sweep gets wrong.
std::string format_frequencies(const std::vector<double>& f) {
    std::ostringstream s;
    s << f.size() << " point" << (f.size() == 1 ? "" : "s") << ":  ";
    if (f.size() <= 6) {
        for (std::size_t i = 0; i < f.size(); ++i) {
            if (i) s << "  ";
            s << f[i];
        }
    } else {
        s << f[0] << "  " << f[1] << "  " << f[2] << "  ...  " << f[f.size() - 2] << "  "
          << f.back();
    }
    s << "  Hz";
    return s.str();
}

void print_summary(const std::string& path, const ParseResult& r, const Mesh& mesh,
                   const BoundProblem& b, const DofMap& d) {
    const Problem& p = r.problem;

    std::cout << "input     " << path << "\n";
    std::cout << "mesh      " << p.mesh_file << "\n";
    std::cout << "          " << mesh.num_nodes() << " nodes, " << mesh.num_edges() << " edges, "
              << mesh.num_faces() << " faces, " << mesh.num_tets() << " tets\n";
    std::cout << "units     " << unit_name(p.length_unit) << "  ->  " << length_scale(p.length_unit)
              << " m per unit\n";

    std::cout << "analysis  ";
    if (p.type == AnalysisType::DC) {
        std::cout << "dc\n";
    } else {
        std::cout << "frequency, " << format_frequencies(p.frequencies) << "\n";
        std::cout << "          formulation "
                  << (p.formulation == Formulation::FullWave ? "full_wave" : "reduced") << "\n";
    }
    std::cout << "          Phi lives on "
              << (phi_on_conductors_only(p) ? "conductors only (sigma > 0)" : "the whole domain")
              << "\n";
    std::cout << "solver    conditioning " << conditioning_keyword(p.conditioning);
    if (p.conditioning == Conditioning::Natural) {
        std::cout << "      unsymmetric matrix; defined at every frequency\n";
    } else {
        std::cout << "   symmetric matrix"
                  << (p.conditioning == Conditioning::ScaledPhi ? "; the unknown is Phi/(j*omega)"
                                                                : "")
                  << "\n";
    }
    std::cout << "boundary  flux_tangential   n x A = 0 on the outer boundary, Phi free\n";

    std::cout << "\nbodies (" << b.bodies.size() << ")\n";
    for (const BoundBody& body : b.bodies) {
        std::ostringstream vol;
        vol << "'" << body.volume << "' (tag " << body.tag << ")";
        std::cout << "  " << std::left << std::setw(6) << body.name << std::setw(22) << vol.str()
                  << std::right << std::setw(7) << body.tets.size() << " tets   sigma "
                  << body.sigma << " S/m  eps_r " << body.eps_r << "  mu_r " << body.mu_r
                  << (body.is_conductor() ? "   [conductor]" : "   [insulator]") << "\n";
    }

    std::cout << "\nports (" << b.ports.size() << ")\n";
    for (std::size_t i = 0; i < b.ports.size(); ++i) {
        const BoundPort& bp = b.ports[i];
        const Port& src = p.ports[i];
        std::cout << "  " << std::left << std::setw(6) << bp.name << std::setw(18)
                  << port_type_name(bp.type) << "'" << join(src.surface) << "'\n";
        std::cout << "        " << bp.surface.faces.size() << " faces, " << bp.surface.nodes.size()
                  << " vertices, " << bp.surface.edges.size() << " edges\n";
        std::cout << "        " << format_amplitude(src);
        if (!is_current_driven(bp.type) && src.amplitude == 0.0) {
            std::cout << "   [potential reference]";
        }
        std::cout << "\n";

        std::cout << "        d = (" << tidy(bp.direction.x) << ", " << tidy(bp.direction.y)
                  << ", " << tidy(bp.direction.z) << ")  ";
        if (!bp.is_internal()) {
            std::cout << "[inward normal -- positive current enters here]\n";
        } else if (bp.direction_from_hint) {
            std::cout << "[from current_direction]\n";
            std::cout << "        plus side: " << bp.plus_side_tet.size()
                      << " faces resolved;  rim: " << bp.rim_edges.size() << " edges\n";
        } else {
            std::cout << "[derived -- the sign of this port's I and V is arbitrary;\n"
                      << "             set current_direction to fix it]\n";
            std::cout << "        plus side: " << bp.plus_side_tet.size()
                      << " faces resolved;  rim: " << bp.rim_edges.size() << " edges\n";
        }
    }

    std::cout << "\nconduction paths (" << b.conduction_paths.size() << ")\n";
    for (std::size_t i = 0; i < b.conduction_paths.size(); ++i) {
        const ConductionPath& path = b.conduction_paths[i];
        std::cout << "  #" << i << "  " << std::setw(7) << path.tets.size() << " tets   bodies";
        for (int bi : path.bodies) std::cout << " " << b.bodies[static_cast<std::size_t>(bi)].name;
        std::cout << "   ports";
        if (path.ports.empty()) {
            std::cout << " (none)";
        } else {
            for (const std::string& n : path.ports) std::cout << " " << n;
        }
        std::cout << "   references " << path.reference_count;
        if (path.is_floating) std::cout << "   [floating -- node " << path.pin_node << " pinned]";
        std::cout << "\n";
    }

    int phi_tets = 0;
    for (bool v : b.phi_tet) {
        if (v) ++phi_tets;
    }
    std::cout << "\nPhi support   " << phi_tets << " of " << mesh.num_tets() << " tets\n";

    // The gauge: every A edge is free, a tree edge (a = 0 by the gauge) or a
    // Dirichlet edge (a = 0 from n x A = 0). The three must account for all
    // of them, which is worth showing rather than asserting quietly.
    int dirichlet = 0;
    for (bool v : b.dirichlet_edge) {
        if (v) ++dirichlet;
    }
    const int free_a = mesh.num_edges() - dirichlet - b.gauge.interior_tree_edge_count;
    std::cout << "gauge         tree-cotree, boundary-first\n"
              << "              " << b.gauge.tree_edge_count << " tree edges ("
              << b.gauge.surface_tree_edge_count << " on n x A = 0 surfaces, "
              << b.gauge.interior_tree_edge_count << " interior)\n"
              << "              " << dirichlet << " Dirichlet edges, "
              << b.gauge.num_dirichlet_components << " surface(s), "
              << b.gauge.num_reference_groups << " root(s)\n"
              << orphan_note(mesh, b)
              << "              " << free_a << " free A unknowns of " << mesh.num_edges()
              << " edges\n";

    std::vector<std::string> warnings = r.warnings;
    warnings.insert(warnings.end(), b.warnings.begin(), b.warnings.end());
    std::cout << "warnings      ";
    if (warnings.empty()) {
        std::cout << "none\n";
    } else {
        std::cout << warnings.size() << "\n";
        for (const std::string& w : warnings) std::cout << "  - " << w << "\n";
    }

    // The DOF map: every edge is exactly one of free / Dirichlet / tree, and
    // every P2 node exactly one of five states, so these account for the
    // whole mesh rather than sampling it.
    int e_free = 0, e_dir = 0, e_tree = 0;
    for (EdgeDof s : d.edge_state) {
        if (s == EdgeDof::Free) ++e_free;
        if (s == EdgeDof::Dirichlet) ++e_dir;
        if (s == EdgeDof::Tree) ++e_tree;
    }
    int n_free = 0, n_absent = 0, n_fixed = 0, n_port = 0, n_cut = 0;
    for (PhiDof s : d.phi_state) {
        if (s == PhiDof::Free) ++n_free;
        if (s == PhiDof::Absent) ++n_absent;
        if (s == PhiDof::Fixed) ++n_fixed;
        if (s == PhiDof::Port) ++n_port;
        if (s == PhiDof::Cut) ++n_cut;
    }

    // std::left is sticky from the tables above, so setw would pad the
    // wrong side of these counts.
    std::cout << std::right;
    std::cout << "\ndegrees of freedom\n"
              << "  A    " << std::setw(7) << d.num_a << " free   of " << mesh.num_edges()
              << " edges  (" << e_dir << " Dirichlet, " << e_tree << " tree)\n"
              << "  Phi  " << std::setw(7) << d.num_phi << " free   of " << d.num_p2_nodes
              << " P2 nodes  (" << n_absent << " absent, " << n_port << " on terminals, " << n_cut
              << " on cuts, " << n_fixed << " pinned)\n"
              << "  V    " << std::setw(7) << d.num_ports << "        one per port\n"
              << "  ----------------------------------------\n"
              << "  total" << std::setw(8) << d.num_total << " unknowns\n";

}

/// Assembles the system at every frequency the file asks for, and reports
/// what each pass produced.
///
/// Uses the sweep path from `docs/ASSEMBLY_PLAN.md` Sec. 11-13 rather than
/// calling `assemble` per frequency: the structure is frequency-independent,
/// so the pattern, the matrix and the scatter positions are built once and
/// only the values are recomputed.
///
/// Storage follows the conditioning. The two symmetric ones get a
/// `SparseSymmetric`, which holds the upper triangle for about half the
/// memory; `Natural` cannot, since its (Phi,A) block is `j*omega` times its
/// (A,Phi) block.
void assemble_all(const Problem& p, const Mesh& mesh, const BoundProblem& bound,
                  const DofMap& dofs) {
    using Clock = std::chrono::steady_clock;
    const auto ms = [](Clock::time_point a, Clock::time_point b) {
        return std::chrono::duration<double, std::milli>(b - a).count();
    };
    const double mb = 1048576.0;

    const bool symmetric = conditioning_needs_ac(p.conditioning);
    const SparsityStorage storage =
        symmetric ? SparsityStorage::UpperTriangle : SparsityStorage::Full;

    auto t0 = Clock::now();
    const SparsityPattern pattern = build_sparsity(dofs, bound, mesh, storage);
    auto t1 = Clock::now();

    std::cout << "\nassembly\n";
    std::cout << "  conditioning  " << conditioning_keyword(p.conditioning) << "  ->  "
              << (symmetric ? "symmetric, upper triangle stored"
                            : "unsymmetric, both triangles stored")
              << "\n";
    std::cout << "  pattern       " << pattern.nnz() << " nonzeros, "
              << static_cast<double>(pattern.nnz()) / dofs.num_total << " per row, "
              << pattern.nnz() * (sizeof(std::complex<double>) + sizeof(int)) / mb
              << " MB of values and columns   [" << ms(t0, t1) << " ms]\n";

    // The list of frequencies to assemble at. DC is one pass at omega = 0.
    std::vector<double> omegas;
    if (p.type == AnalysisType::DC) {
        omegas.push_back(0.0);
    } else {
        for (double f : p.frequencies) omegas.push_back(2.0 * M_PI * f);
    }

    // A ScatterMap is built once and reused, and pays for itself after the
    // first refill -- so it is worth it from two frequencies up, and roughly
    // a wash for one. ASSEMBLY_PLAN Sec. 13 has the numbers.
    std::unique_ptr<ScatterMap> map;
    double map_ms = 0.0;
    if (omegas.size() > 1) {
        auto m0 = Clock::now();
        map.reset(new ScatterMap(ScatterMap::build(pattern, dofs, mesh, bound)));
        auto m1 = Clock::now();
        map_ms = ms(m0, m1);
        std::cout << "  scatter map   " << map->bytes() / mb << " MB, reused at every frequency   ["
                  << map_ms << " ms]\n";
    } else {
        std::cout << "  scatter map   not built: one pass only, so it would cost more than it "
                     "saves\n";
    }

    auto a0 = Clock::now();
    AssembledSystem full;
    SymmetricSystem half;
    if (symmetric) {
        half = make_symmetric_system(pattern, dofs.num_total);
    } else {
        full = make_system(pattern, dofs.num_total);
    }
    auto a1 = Clock::now();
    std::cout << "  allocated     "
              << (symmetric ? half.matrix.nnz() : full.matrix.nnz())
              << " complex values once for the whole sweep   [" << ms(a0, a1) << " ms]\n\n";

    // `max |M_ij|` and `max |b_i|`, NOT anything about the A field -- an
    // earlier version of this table called the first column "largest |a|",
    // which reads as the vector potential, since `a` is what the A unknowns
    // are called everywhere else in this project.
    //
    // Neither is a correctness check. They are the cheapest signals that the
    // assembly produced something and that its scaling is what the
    // conditioning predicts: under row_scaled the port row's entry comes out
    // as 1/omega, and under scaled_phi the largest coefficient grows with
    // omega. A column of zeros would mean nothing was scattered at all.
    std::cout << std::right;
    std::cout << "  " << std::setw(12) << "frequency" << std::setw(14) << "max |M_ij|"
              << std::setw(14) << "max |b_i|" << std::setw(10) << "ms\n";
    double total_ms = 0.0;
    for (double omega : omegas) {
        auto r0 = Clock::now();
        if (symmetric) {
            refill(half, bound, mesh, dofs, pattern, omega, p.conditioning, map.get());
        } else {
            refill(full, bound, mesh, dofs, pattern, omega, p.conditioning, map.get());
        }
        auto r1 = Clock::now();
        total_ms += ms(r0, r1);

        double worst = 0.0, worst_rhs = 0.0;
        const std::vector<std::complex<double>>& values =
            symmetric ? half.matrix.values() : full.matrix.values();
        const std::vector<std::complex<double>>& rhs = symmetric ? half.rhs : full.rhs;
        for (const std::complex<double>& v : values) worst = std::max(worst, std::abs(v));
        for (const std::complex<double>& v : rhs) worst_rhs = std::max(worst_rhs, std::abs(v));

        std::ostringstream label;
        if (omega == 0.0) {
            label << "DC";
        } else {
            label << std::setprecision(4) << omega / (2.0 * M_PI) << " Hz";
        }
        std::cout << "  " << std::setw(12) << label.str() << std::setw(14) << worst
                  << std::setw(14) << worst_rhs << std::setw(10) << ms(r0, r1) << "\n";
    }

    std::cout << "\n  " << omegas.size() << " assembl" << (omegas.size() == 1 ? "y" : "ies")
              << " in " << total_ms << " ms";
    if (omegas.size() > 1) {
        std::cout << ", " << total_ms / omegas.size() << " ms each";
    }
    std::cout << "\n";

    std::cout << "\nok -- assembled. Nothing was SOLVED: there is no linear solver yet, so no\n"
              << "currents, voltages, R or L come out of this. Next in the plan: the solver.\n";
}

}  // namespace

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "A-Phi solver " << aphi_solver::kVersion << " -- assembles, does not solve\n\n"
                  << "usage: " << (argc > 0 ? argv[0] : "aphi_solver") << " <input-file>\n\n"
                  << "Reads and validates an input file, opens its mesh, numbers the unknowns\n"
                  << "and assembles the system at every frequency. It does not solve.\n\n"
                  << "try: aphi_solver examples/cylinder_box.aphi\n";
        return 2;
    }

    const std::string path = argv[1];
    try {
        const aphi_solver::ParseResult result = aphi_solver::parse_input_file(path);

        aphi_solver::Mesh mesh = aphi_solver::read_gmsh_msh(result.problem.mesh_file);
        aphi_solver::scale_mesh_to_metres(mesh, result.problem.length_unit);
        const aphi_solver::BoundProblem bound = aphi_solver::bind_to_mesh(result.problem, mesh);
        const aphi_solver::DofMap dofs = aphi_solver::build_dof_map(bound, mesh);

        std::cout << "A-Phi solver " << aphi_solver::kVersion
                  << " -- reads, checks and assembles; does not solve it yet\n\n";
        print_summary(path, result, mesh, bound, dofs);
        assemble_all(result.problem, mesh, bound, dofs);
        return 0;
    } catch (const aphi_solver::GmshReadError& e) {
        std::cerr << "mesh: " << e.what() << "\n";
        return 1;
    } catch (const aphi_solver::InputError& e) {
        // The message already carries "file:line: ", which is the form an
        // editor can jump to.
        std::cerr << e.what() << "\n";
        return 1;
    } catch (const std::exception& e) {
        std::cerr << path << ": " << e.what() << "\n";
        return 1;
    }
}
