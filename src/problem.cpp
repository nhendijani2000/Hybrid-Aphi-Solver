#include "aphi_solver/problem.hpp"

#include <cmath>
#include <stdexcept>

namespace aphi_solver {

double length_scale(LengthUnit unit) {
    switch (unit) {
        case LengthUnit::Metre:      return 1.0;
        case LengthUnit::Millimetre: return 1e-3;
        case LengthUnit::Micrometre: return 1e-6;
        case LengthUnit::Nanometre:  return 1e-9;
    }
    throw std::invalid_argument("length_scale: unknown LengthUnit");
}

bool is_internal(PortType type) {
    return type == PortType::InternalCurrent || type == PortType::InternalVoltage;
}

bool is_current_driven(PortType type) {
    return type == PortType::BoundaryCurrent || type == PortType::InternalCurrent;
}

std::complex<double> port_amplitude(const Port& port) {
    const double radians = port.phase_deg * M_PI / 180.0;
    // Built component-wise rather than with std::polar. A negative amplitude
    // is legal here (-1 at 0 degrees == +1 at 180), and std::polar's
    // behaviour is *undefined* for a negative magnitude -- MSVC happens to
    // give the same answer, so tests/test_problem.cpp cannot catch the
    // difference on this compiler, but the standard does not promise it
    // anywhere else. Defined behaviour is the whole reason for the longer
    // spelling.
    return {port.amplitude * std::cos(radians), port.amplitude * std::sin(radians)};
}

namespace {

Vec3 normalized(const Vec3& v, const char* who) {
    const double n = v.norm();
    if (n == 0.0) throw std::invalid_argument(std::string(who) + ": zero vector");
    return v * (1.0 / n);
}

}  // namespace

Vec3 canonical_orientation(const Vec3& normal) {
    const Vec3 n = normalized(normal, "canonical_orientation");

    // Largest-magnitude component decides the sign, with ties broken in
    // x, y, z order -- so the strict > comparisons below are deliberate:
    // an exact tie leaves the earlier component in charge.
    double largest = n.x;
    if (std::abs(n.y) > std::abs(largest)) largest = n.y;
    if (std::abs(n.z) > std::abs(largest)) largest = n.z;

    return largest < 0.0 ? n * -1.0 : n;
}

Vec3 resolve_cut_direction(const Vec3& hint, const Vec3& normal, double min_alignment) {
    const Vec3 h = normalized(hint, "resolve_cut_direction: hint");
    const Vec3 n = normalized(normal, "resolve_cut_direction: normal");

    const double alignment = h.dot(n);
    if (std::abs(alignment) < min_alignment) {
        throw std::invalid_argument(
            "resolve_cut_direction: the current-direction hint is nearly parallel to the cut "
            "plane, so it does not say which way through the cut the current flows");
    }
    return alignment > 0.0 ? n : n * -1.0;
}

bool phi_on_conductors_only(const Problem& problem) {
    return problem.type == AnalysisType::DC || problem.formulation == Formulation::Reduced;
}

std::size_t num_solves(const Problem& problem) {
    return problem.type == AnalysisType::DC ? 1u : problem.frequencies.size();
}

}  // namespace aphi_solver
