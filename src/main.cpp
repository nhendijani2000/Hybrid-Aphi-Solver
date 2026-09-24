// Command-line driver.
//
// **Parse-only, as of step 2 of `Claude outputs/input_file_plan.md`.** It
// reads an input file, validates everything answerable from the text alone,
// and reports what it found. It does not open the mesh and it does not
// solve: binding is step 4 and the DOF map comes after it. Everything below
// that would need the mesh is marked as such rather than guessed at, so the
// output never implies more has been checked than actually has.

#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "aphi_solver/input_file.hpp"
#include "aphi_solver/version.hpp"

namespace {

using namespace aphi_solver;

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

void print_summary(const std::string& path, const ParseResult& r) {
    const Problem& p = r.problem;

    std::cout << "input     " << path << "\n";
    std::cout << "mesh      " << p.mesh_file << "   (not opened at this stage)\n";
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
    std::cout << "boundary  flux_tangential   n x A = 0 on the outer boundary, Phi free\n";

    std::cout << "\nbodies (" << p.bodies.size() << ")\n";
    for (const Body& b : p.bodies) {
        std::cout << "  " << std::left << std::setw(6) << b.name << "volume '" << b.volume << "'"
                  << std::string(b.volume.size() < 12 ? 12 - b.volume.size() : 1, ' ')
                  << "sigma " << std::setw(10) << b.sigma << " S/m   eps_r " << b.eps_r << "   mu_r "
                  << b.mu_r << (b.sigma > 0.0 ? "   [conductor]" : "   [insulator]") << "\n";
    }

    std::cout << "\nports (" << p.ports.size() << ")\n";
    for (const Port& port : p.ports) {
        std::cout << "  " << std::left << std::setw(6) << port.name << std::setw(18)
                  << port_type_name(port.type) << "surface " << join(port.surface) << "\n";
        std::cout << "        " << format_amplitude(port);
        if (!is_current_driven(port.type) && port.amplitude == 0.0) {
            std::cout << "   [potential reference]";
        }
        std::cout << "\n";
        if (is_internal(port.type)) {
            if (port.current_direction.has_value()) {
                const Vec3& d = *port.current_direction;
                std::cout << "        current_direction hint (" << d.x << ", " << d.y << ", " << d.z
                          << ")  -- resolved against the cut's normal when the mesh is read\n";
            } else {
                std::cout << "        no current_direction -- will be derived from the mesh\n";
            }
        }
    }

    std::cout << "\nwarnings  ";
    if (r.warnings.empty()) {
        std::cout << "none\n";
    } else {
        std::cout << r.warnings.size() << "\n";
        for (const std::string& w : r.warnings) std::cout << "  - " << w << "\n";
    }

    std::cout << "\nok -- this file describes " << num_solves(p) << " solve"
              << (num_solves(p) == 1 ? "" : "s") << ".\n"
              << "Nothing was solved: the mesh is not read and no DOFs exist yet.\n"
              << "Next in the plan: bind to the mesh (step 4), then the DOF map.\n";
}

}  // namespace

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "A-Phi solver " << aphi_solver::kVersion << " -- parse-only driver\n\n"
                  << "usage: " << (argc > 0 ? argv[0] : "aphi_solver") << " <input-file>\n\n"
                  << "Reads and validates an input file, then prints what it found.\n"
                  << "It does not open the mesh and does not solve.\n\n"
                  << "try: aphi_solver examples/cylinder_box.aphi\n";
        return 2;
    }

    const std::string path = argv[1];
    try {
        const aphi_solver::ParseResult result = aphi_solver::parse_input_file(path);
        std::cout << "A-Phi solver " << aphi_solver::kVersion
                  << " -- parse-only (no mesh is read, nothing is solved)\n\n";
        print_summary(path, result);
        return 0;
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
