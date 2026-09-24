// Tests for the Problem description and its conventions.
//
// Step 1 of `Claude outputs/input_file_plan.md` Sec. 6. No parser and no
// mesh here on purpose. A convention error here -- an excitation entering
// the system with the wrong sign or a phase applied backwards -- would
// surface much later as a current of the wrong sign or an impedance with a
// flipped imaginary part, at which point it is indistinguishable from an
// assembly bug. This is where it is cheapest to catch.
//
// Minimal hand-written runner, matching the rest of tests/ (see the comment
// in tests/CMakeLists.txt for why no external framework).

#include <cmath>
#include <complex>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>

#include "aphi_solver/problem.hpp"

using aphi_solver::AnalysisType;
using aphi_solver::Body;
using aphi_solver::Formulation;
using aphi_solver::LengthUnit;
using aphi_solver::Port;
using aphi_solver::PortType;
using aphi_solver::Problem;
using aphi_solver::Vec3;
using Complex = std::complex<double>;

namespace {

int g_checks = 0;
int g_failures = 0;

void check(bool condition, const std::string& what) {
    ++g_checks;
    if (!condition) {
        ++g_failures;
        std::cerr << "FAIL: " << what << "\n";
    }
}

bool near(double a, double b, double tol = 1e-12) { return std::abs(a - b) <= tol; }
bool near(Complex a, Complex b, double tol = 1e-12) { return std::abs(a - b) <= tol; }

Port make_internal(PortType type, double amplitude) {
    Port p;
    p.name = "P";
    p.type = type;
    p.surface = {"cut"};
    p.amplitude = amplitude;
    return p;
}

Port make_boundary(PortType type, double amplitude) {
    Port p;
    p.name = "P";
    p.type = type;
    p.surface = {"terminal"};
    p.amplitude = amplitude;
    return p;
}

// ---------------------------------------------------------------------------

void test_port_classification() {
    check(!aphi_solver::is_internal(PortType::BoundaryCurrent), "boundary_current is not internal");
    check(!aphi_solver::is_internal(PortType::BoundaryVoltage), "boundary_voltage is not internal");
    check(aphi_solver::is_internal(PortType::InternalCurrent), "internal_current is internal");
    check(aphi_solver::is_internal(PortType::InternalVoltage), "internal_voltage is internal");

    check(aphi_solver::is_current_driven(PortType::BoundaryCurrent), "boundary_current is current driven");
    check(aphi_solver::is_current_driven(PortType::InternalCurrent), "internal_current is current driven");
    check(!aphi_solver::is_current_driven(PortType::BoundaryVoltage), "boundary_voltage is not current driven");
    check(!aphi_solver::is_current_driven(PortType::InternalVoltage), "internal_voltage is not current driven");
}

// The orientation convention (24 Sept 2026): d points from the grounded
// minus side to the live plus side, so V = Phi(plus) - 0 = Phi(plus) and
// the excitation enters the system unmodified. There is no sign factor to
// test -- which is the property worth asserting, because an earlier design
// had one and it is exactly the kind of thing that grows back.
void test_excitation_carries_no_sign() {
    // Every port type: what the file said is what the assembly gets.
    check(near(port_amplitude(make_internal(PortType::InternalCurrent, 1.0)), Complex(1.0, 0.0)),
          "internal_current: excitation is +I unmodified");
    check(near(port_amplitude(make_internal(PortType::InternalVoltage, 1.0)), Complex(1.0, 0.0)),
          "internal_voltage: excitation is +V unmodified");
    check(near(port_amplitude(make_boundary(PortType::BoundaryCurrent, 1.0)), Complex(1.0, 0.0)),
          "boundary_current: excitation is +I unmodified");
    check(near(port_amplitude(make_boundary(PortType::BoundaryVoltage, 0.0)), Complex(0.0, 0.0)),
          "boundary_voltage: 0 V reference is 0");

    // The excitation must not depend on the port's *kind* -- only on its
    // amplitude and phase. If internal and boundary ports ever diverge here,
    // a sign factor has crept back in.
    check(near(port_amplitude(make_internal(PortType::InternalCurrent, -7.5)),
               port_amplitude(make_boundary(PortType::BoundaryCurrent, -7.5))),
          "excitation depends on amplitude and phase alone, not on port kind");

    // The convention closes: with the minus side grounded at 0 and the plus
    // side carrying the excitation, V = Phi(plus) - Phi(minus) is what was
    // asked for.
    const Complex phi_minus(0.0, 0.0);
    const Complex phi_plus = port_amplitude(make_internal(PortType::InternalVoltage, 1.0));
    check(near(phi_plus - phi_minus, Complex(1.0, 0.0)),
          "V = Phi(plus) - Phi(minus) recovers the requested 1 V");
}

void test_complex_amplitude() {
    Port p = make_boundary(PortType::BoundaryCurrent, 1.0);

    p.phase_deg = 0.0;
    check(near(port_amplitude(p), Complex(1.0, 0.0)), "1.0 at 0 deg is +1");

    p.phase_deg = 90.0;
    check(near(port_amplitude(p), Complex(0.0, 1.0)), "1.0 at 90 deg is +j");

    p.phase_deg = 180.0;
    check(near(port_amplitude(p), Complex(-1.0, 0.0)), "1.0 at 180 deg is -1");

    p.phase_deg = -90.0;
    check(near(port_amplitude(p), Complex(0.0, -1.0)), "1.0 at -90 deg is -j");

    // A negative magnitude is legal and must not be silently absorbed.
    //
    // Note on what this does NOT pin (checked by negative control, Sept
    // 2026): swapping the implementation to std::polar, whose behaviour is
    // *undefined* for a negative magnitude, leaves this test green on MSVC,
    // which happens to return the same value. So the component-wise
    // construction in problem.cpp is a portability choice this test cannot
    // enforce -- it is defined behaviour where std::polar is not, and that
    // is the reason to keep it, not because a check here would catch the
    // difference on this compiler.
    Port negative = make_boundary(PortType::BoundaryCurrent, -1.0);
    negative.phase_deg = 0.0;
    Port rotated = make_boundary(PortType::BoundaryCurrent, 1.0);
    rotated.phase_deg = 180.0;
    check(near(port_amplitude(negative), port_amplitude(rotated)),
          "-1.0 at 0 deg equals 1.0 at 180 deg");
    check(near(port_amplitude(negative), Complex(-1.0, 0.0)), "-1.0 at 0 deg is -1");

    // Magnitude and phase together, against a hand value.
    Port both = make_boundary(PortType::BoundaryVoltage, 3.3);
    both.phase_deg = 45.0;
    const double r = 3.3 / std::sqrt(2.0);
    check(near(port_amplitude(both), Complex(r, r), 1e-12), "3.3 at 45 deg is 3.3/sqrt(2) * (1+j)");

    // An internal port's phase is applied exactly as a boundary port's is:
    // nothing about the cut geometry reaches back into the amplitude.
    Port internal = make_internal(PortType::InternalCurrent, 1.0);
    internal.phase_deg = 90.0;
    check(near(port_amplitude(internal), Complex(0.0, 1.0)),
          "an internal port's 90 deg phase gives +j, same as a boundary port's");
}

bool near(const aphi_solver::Vec3& a, const aphi_solver::Vec3& b, double tol = 1e-12) {
    return std::abs(a.x - b.x) <= tol && std::abs(a.y - b.y) <= tol && std::abs(a.z - b.z) <= tol;
}

bool throws_invalid_argument(const std::function<void()>& f) {
    try {
        f();
    } catch (const std::invalid_argument&) {
        return true;
    } catch (...) {
        return false;
    }
    return false;
}

// The hint only ever contributes a sign, so any vector within 90 degrees of
// the intended direction must give the same answer. That is what lets a user
// say "roughly +y" instead of computing a mesh normal.
void test_resolve_cut_direction() {
    using aphi_solver::resolve_cut_direction;
    const Vec3 n(0.0, 1.0, 0.0);

    check(near(resolve_cut_direction(Vec3(0.0, 1.0, 0.0), n), Vec3(0.0, 1.0, 0.0)),
          "hint along the normal keeps it");
    check(near(resolve_cut_direction(Vec3(0.0, -1.0, 0.0), n), Vec3(0.0, -1.0, 0.0)),
          "hint against the normal flips it");

    // The point of a hint: sloppy input, same answer. 0.1 0.9 -0.2 and a
    // clean +y must be indistinguishable.
    check(near(resolve_cut_direction(Vec3(0.1, 0.9, -0.2), n), Vec3(0.0, 1.0, 0.0)),
          "a sloppy hint within 90 degrees gives the same direction as an exact one");
    check(near(resolve_cut_direction(Vec3(-0.3, -0.8, 0.5), n), Vec3(0.0, -1.0, 0.0)),
          "a sloppy hint the other way flips it");

    // Magnitude is irrelevant -- only the sign of the dot product is used.
    check(near(resolve_cut_direction(Vec3(0.0, 1e-9, 0.0), n), Vec3(0.0, 1.0, 0.0)),
          "a tiny hint still resolves; only its direction matters");
    check(near(resolve_cut_direction(Vec3(0.0, 1e9, 0.0), n), Vec3(0.0, 1.0, 0.0)),
          "a huge hint resolves identically");

    // The returned direction is the *normal*, normalised -- never the hint.
    check(near(resolve_cut_direction(Vec3(1.0, 2.0, 0.0), Vec3(0.0, 5.0, 0.0)), Vec3(0.0, 1.0, 0.0)),
          "the result is the normalised normal, not the hint");

    // A hint lying in the cut plane picks no side.
    check(throws_invalid_argument([&] { resolve_cut_direction(Vec3(1.0, 0.0, 0.0), n); }),
          "a hint perpendicular to the normal is rejected");
    check(throws_invalid_argument([&] { resolve_cut_direction(Vec3(1.0, 0.05, 0.0), n); }),
          "a hint only 3 degrees out of the cut plane is rejected");
    check(!throws_invalid_argument([&] { resolve_cut_direction(Vec3(1.0, 0.5, 0.0), n); }),
          "a hint comfortably off the cut plane is accepted");

    check(throws_invalid_argument([&] { resolve_cut_direction(Vec3(0.0, 0.0, 0.0), n); }),
          "a zero hint is rejected");
    check(throws_invalid_argument([&] { resolve_cut_direction(n, Vec3(0.0, 0.0, 0.0)); }),
          "a zero normal is rejected");
}

// Used only when no hint is given. The requirement is reproducibility, not
// correctness: the same normal must always yield the same direction, so a
// test's reported sign cannot drift between builds.
void test_canonical_orientation() {
    using aphi_solver::canonical_orientation;

    check(near(canonical_orientation(Vec3(0.0, 1.0, 0.0)), Vec3(0.0, 1.0, 0.0)), "+y stays +y");
    check(near(canonical_orientation(Vec3(0.0, -1.0, 0.0)), Vec3(0.0, 1.0, 0.0)), "-y becomes +y");
    check(near(canonical_orientation(Vec3(0.0, 0.0, -3.0)), Vec3(0.0, 0.0, 1.0)), "-3z becomes +z");

    // Opposite normals of one plane must collapse to the same direction --
    // that is what makes the choice independent of node numbering.
    check(near(canonical_orientation(Vec3(1.0, -2.0, 0.5)), canonical_orientation(Vec3(-1.0, 2.0, -0.5))),
          "a normal and its negation give the same canonical direction");

    // The largest-magnitude component ends up positive.
    const Vec3 d = canonical_orientation(Vec3(0.3, -0.9, 0.2));
    check(d.y > 0.0, "the largest component is made positive");
    check(near(d.norm(), 1.0), "the result is a unit vector");

    // Idempotent, and ties break in x, y, z order.
    check(near(canonical_orientation(d), d), "canonical_orientation is idempotent");
    check(near(canonical_orientation(Vec3(-1.0, 1.0, 0.0)), Vec3(1.0, -1.0, 0.0) * (1.0 / std::sqrt(2.0))),
          "an exact tie leaves the earlier component (x) in charge");

    check(throws_invalid_argument([&] { canonical_orientation(Vec3(0.0, 0.0, 0.0)); }),
          "a zero normal is rejected");
}

// The struct must be able to carry the hint, and absence must be
// distinguishable from a zero vector -- the first means "derive it", the
// second is bad input.
void test_current_direction_is_optional() {
    Port without = make_internal(PortType::InternalCurrent, 1.0);
    check(!without.current_direction.has_value(), "a hint is absent by default");

    Port with = make_internal(PortType::InternalCurrent, 1.0);
    with.current_direction = Vec3(0.0, 1.0, 0.0);
    check(with.current_direction.has_value(), "a hint can be given");
    check(near(*with.current_direction, Vec3(0.0, 1.0, 0.0)), "the hint round-trips");

    Port zero = make_internal(PortType::InternalCurrent, 1.0);
    zero.current_direction = Vec3(0.0, 0.0, 0.0);
    check(zero.current_direction.has_value(),
          "an explicit zero hint is present, not absent -- so it can be rejected as bad input "
          "rather than silently derived");
}

// The Sec. 1.2 table: where Phi lives is forced at DC, chosen only between
// the two frequency-domain formulations.
void test_phi_support() {
    Problem dc;
    dc.type = AnalysisType::DC;
    check(phi_on_conductors_only(dc), "DC: Phi on conductors only");

    // Even if a formulation were somehow set, DC still forces it -- the Phi
    // rows are structurally empty where sigma = 0 regardless.
    dc.formulation = Formulation::FullWave;
    check(phi_on_conductors_only(dc), "DC forces conductors-only whatever the formulation says");

    Problem full;
    full.type = AnalysisType::Frequency;
    full.formulation = Formulation::FullWave;
    check(!phi_on_conductors_only(full), "full wave: Phi everywhere");

    Problem reduced;
    reduced.type = AnalysisType::Frequency;
    reduced.formulation = Formulation::Reduced;
    check(phi_on_conductors_only(reduced), "reduced: Phi on conductors only");
}

void test_length_scale() {
    check(near(length_scale(LengthUnit::Metre), 1.0), "m -> 1");
    check(near(length_scale(LengthUnit::Millimetre), 1e-3), "mm -> 1e-3");
    check(near(length_scale(LengthUnit::Micrometre), 1e-6), "um -> 1e-6");
    check(near(length_scale(LengthUnit::Nanometre), 1e-9), "nm -> 1e-9");
}

void test_num_solves() {
    Problem dc;
    dc.type = AnalysisType::DC;
    check(num_solves(dc) == 1u, "DC is one solve");

    Problem sweep;
    sweep.type = AnalysisType::Frequency;
    sweep.frequencies = {1e3, 1e4, 1e5, 1e6, 1e7};
    check(num_solves(sweep) == 5u, "a five-point sweep is five solves");
}

// The cylinder test case of Sec. 5, assembled by hand -- a check that the
// struct can actually express the first problem the solver will run, before
// any parser exists to build it.
void test_cylinder_problem_is_expressible() {
    Problem p;
    p.mesh_file = "../meshes/cylinder_box.msh";
    p.length_unit = LengthUnit::Millimetre;
    p.type = AnalysisType::DC;

    Body wire;
    wire.name = "B1";
    wire.volume = "wire";
    wire.sigma = 5.8e7;
    Body air;
    air.name = "B2";
    air.volume = "air";
    air.sigma = 0.0;
    p.bodies = {wire, air};

    Port feed = make_boundary(PortType::BoundaryCurrent, 1.0);
    feed.name = "P1";
    feed.surface = {"wire_bottom"};
    Port ground = make_boundary(PortType::BoundaryVoltage, 0.0);
    ground.name = "P2";
    ground.surface = {"wire_top"};
    p.ports = {feed, ground};

    check(p.bodies.size() == 2u && p.ports.size() == 2u, "cylinder problem has 2 bodies and 2 ports");
    check(phi_on_conductors_only(p), "cylinder at DC puts Phi in the wire only");
    check(near(length_scale(p.length_unit), 1e-3), "cylinder mesh is in mm");
    check(near(port_amplitude(p.ports[0]), Complex(1.0, 0.0)), "cylinder P1 drives +1 A");
    check(near(port_amplitude(p.ports[1]), Complex(0.0, 0.0)), "cylinder P2 is the 0 V reference");
    check(p.bodies[1].eps_r == 1.0 && p.bodies[1].mu_r == 1.0, "air defaults to eps_r = mu_r = 1");
}

}  // namespace

int main() {
    test_port_classification();
    test_excitation_carries_no_sign();
    test_complex_amplitude();
    test_resolve_cut_direction();
    test_canonical_orientation();
    test_current_direction_is_optional();
    test_phi_support();
    test_length_scale();
    test_num_solves();
    test_cylinder_problem_is_expressible();

    std::cout << (g_checks - g_failures) << "/" << g_checks << " checks passed\n";
    return g_failures == 0 ? 0 : 1;
}
