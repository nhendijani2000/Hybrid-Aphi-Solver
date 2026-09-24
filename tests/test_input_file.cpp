// Tests for the input-file parser.
//
// Step 2 of `Claude outputs/input_file_plan.md` Sec. 6. Every fixture is a
// string literal beside its assertions rather than a temp file, so a case
// reads as one piece.
//
// The pattern used throughout for error cases: assert the error fires, AND
// assert it reports the right LINE. An error pointing at the wrong line is
// nearly as unhelpful as no error. Each group is also paired with a
// negative control -- the same file with only the offending detail fixed
// must parse cleanly -- because a parser that rejected everything would
// otherwise pass every rejection test here.

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "aphi_solver/input_file.hpp"

using aphi_solver::AnalysisType;
using aphi_solver::Formulation;
using aphi_solver::InputError;
using aphi_solver::LengthUnit;
using aphi_solver::ParseResult;
using aphi_solver::PortType;
using aphi_solver::Problem;
using aphi_solver::Vec3;

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

// Asserts that `text` is rejected, on line `expected_line`. Returns the
// message so a caller can check its wording.
std::string expect_error(const std::string& text, int expected_line, const std::string& what) {
    try {
        aphi_solver::parse_input_string(text);
    } catch (const InputError& e) {
        check(e.line() == expected_line,
              what + " -- reported on line " + std::to_string(e.line()) + ", expected " +
                  std::to_string(expected_line) + " (" + e.message() + ")");
        return e.message();
    } catch (const std::exception& e) {
        check(false, what + " -- threw the wrong exception type: " + e.what());
        return {};
    }
    check(false, what + " -- was accepted, but should have been rejected");
    return {};
}

ParseResult expect_ok(const std::string& text, const std::string& what) {
    try {
        return aphi_solver::parse_input_string(text);
    } catch (const std::exception& e) {
        check(false, what + " -- was rejected: " + e.what());
        return {};
    }
}

// The Sec. 1 cylinder example, verbatim. Line 1 is the comment.
const char* kCylinderDC = R"(# Cylinder in a square box -- DC resistance and inductance.

[mesh]
file        = ../meshes/cylinder_box.msh
length_unit = mm

[analysis]
type        = dc

[boundary]
outer       = flux_tangential

[Body B1]
volume      = wire
sigma       = 5.8e7
eps_r       = 1
mu_r        = 1

[Body B2]
volume      = air
sigma       = 0

[port P1]
type        = boundary_current
surface     = wire_bottom
current     = 1.0

[port P2]
type        = boundary_voltage
surface     = wire_top
voltage     = 0.0
)";

// ---------------------------------------------------------------------------

void test_cylinder_example_parses() {
    const ParseResult r = expect_ok(kCylinderDC, "the Sec. 1 DC cylinder example");
    const Problem& p = r.problem;

    check(p.mesh_file == "../meshes/cylinder_box.msh", "mesh file is read verbatim");
    check(p.length_unit == LengthUnit::Millimetre, "length_unit = mm");
    check(p.type == AnalysisType::DC, "type = dc");
    check(p.frequencies.empty(), "DC carries no frequencies");
    check(aphi_solver::num_solves(p) == 1u, "DC is one solve");
    check(aphi_solver::phi_on_conductors_only(p), "DC puts Phi on conductors only");

    check(p.bodies.size() == 2u, "two bodies");
    check(p.bodies[0].name == "B1" && p.bodies[0].volume == "wire", "B1 is the wire");
    check(near(p.bodies[0].sigma, 5.8e7), "copper sigma");
    check(near(p.bodies[1].sigma, 0.0), "air sigma");
    check(near(p.bodies[1].eps_r, 1.0) && near(p.bodies[1].mu_r, 1.0),
          "eps_r and mu_r default to 1 when omitted");

    check(p.ports.size() == 2u, "two ports");
    check(p.ports[0].name == "P1" && p.ports[0].type == PortType::BoundaryCurrent, "P1 drives current");
    check(p.ports[0].surface.size() == 1u && p.ports[0].surface[0] == "wire_bottom", "P1 surface");
    check(near(p.ports[0].amplitude, 1.0), "P1 is 1 A");
    check(near(p.ports[0].phase_deg, 0.0), "phase defaults to 0");
    check(p.ports[1].type == PortType::BoundaryVoltage && near(p.ports[1].amplitude, 0.0),
          "P2 is the 0 V reference");
    check(!p.ports[0].current_direction.has_value(), "boundary ports carry no direction");

    check(r.warnings.empty(), "the cylinder example produces no warnings");
}

void test_frequency_example_parses() {
    const ParseResult r = expect_ok(R"(
[mesh]
file = m.msh
length_unit = mm

[analysis]
type        = frequency
frequencies = 1e3 1e4 1e5 1e6 1e7
formulation = full_wave

[Body B1]
volume = wire
sigma  = 5.8e7

[port P1]
type    = boundary_current
current = 1.0

surface = wire_bottom
[port P2]
type    = boundary_voltage
surface = wire_top
voltage = 0.0
)",
                                   "a frequency sweep");
    check(r.problem.type == AnalysisType::Frequency, "type = frequency");
    check(r.problem.frequencies.size() == 5u, "five frequencies");
    check(near(r.problem.frequencies[0], 1e3) && near(r.problem.frequencies[4], 1e7),
          "the listed frequencies are kept exactly");
    check(r.problem.formulation == Formulation::FullWave, "formulation = full_wave");
    check(!aphi_solver::phi_on_conductors_only(r.problem), "full wave puts Phi everywhere");
    check(aphi_solver::num_solves(r.problem) == 5u, "five solves");
}

void test_lexing() {
    // Comments anywhere, blank lines, odd spacing, CRLF, case folding.
    const ParseResult r = expect_ok(
        "# leading comment\r\n"
        "\r\n"
        "[MESH]\r\n"
        "   FILE   =   m.msh   # trailing comment\r\n"
        "length_unit=MM\r\n"
        "\r\n"
        "[Analysis]\r\n"
        "TYPE = DC\r\n"
        "\r\n"
        "[body B1]\r\n"
        "volume = Wire\r\n"
        "sigma  = 1\r\n"
        "\r\n"
        "[PORT P1]\r\n"
        "type    = BOUNDARY_VOLTAGE\r\n"
        "surface = Wire_Top\r\n"
        "voltage = 0\r\n",
        "comments, blank lines, CRLF and case folding");

    check(r.problem.mesh_file == "m.msh", "keys and section kinds fold case");
    check(r.problem.length_unit == LengthUnit::Millimetre, "keyword values fold case");
    check(r.problem.type == AnalysisType::DC, "type folds case");
    check(r.problem.bodies[0].volume == "Wire", "names keep their case -- Gmsh does not fold");
    check(r.problem.ports[0].surface[0] == "Wire_Top", "surface names keep their case");
    check(r.problem.ports[0].type == PortType::BoundaryVoltage, "port type folds case");
}

void test_syntax_errors() {
    expect_error("[mesh\nfile = m.msh\n", 1, "unclosed section header");
    expect_error("[]\n", 1, "empty section header");
    expect_error("[port P1 extra]\n", 1, "three-word section header");
    expect_error("file = m.msh\n", 1, "entry before any section");
    expect_error("[mesh]\nnonsense\n", 2, "line with no '='");
    expect_error("[mesh]\n= m.msh\n", 2, "missing key");
    expect_error("[mesh]\nfile =\n", 2, "missing value");
    expect_error("[mesh]\nfile = a\nfile = b\n", 3, "duplicate key reports the second one");
    expect_error("[wrong]\nx = 1\n", 1, "unknown section");
    expect_error("[solver]\nx = 1\n", 1, "reserved [solver] section");
    expect_error("[output]\nx = 1\n", 1, "reserved [output] section");
}

// Every rejection below is paired with the same file made valid, so the
// checks cannot be passed by a parser that simply rejects everything.
void test_value_errors_and_controls() {
    const std::string head =
        "[mesh]\nfile = m.msh\nlength_unit = mm\n[analysis]\ntype = dc\n";
    const std::string tail =
        "[port P1]\ntype = boundary_voltage\nsurface = s\nvoltage = 0\n";

    // sigma
    expect_error(head + "[Body B1]\nvolume = w\nsigma = -1\n" + tail, 8, "negative sigma");
    expect_ok(head + "[Body B1]\nvolume = w\nsigma = 0\n" + tail, "control: sigma = 0 is fine");

    expect_error(head + "[Body B1]\nvolume = w\nsigma = copper\n" + tail, 8, "non-numeric sigma");
    expect_error(head + "[Body B1]\nvolume = w\nsigma = 1 2\n" + tail, 8, "two values for sigma");
    expect_error(head + "[Body B1]\nvolume = w\nsigma = 1\neps_r = 0\n" + tail, 9, "eps_r = 0");
    expect_error(head + "[Body B1]\nvolume = w\nsigma = 1\nmu_r = -2\n" + tail, 9, "negative mu_r");
    expect_ok(head + "[Body B1]\nvolume = w\nsigma = 1\neps_r = 4.3\nmu_r = 1\n" + tail,
              "control: positive eps_r and mu_r");

    // unknown key, and required keys
    expect_error(head + "[Body B1]\nvolume = w\nsigma = 1\nsigmaa = 2\n" + tail, 9, "misspelled key");
    expect_error(head + "[Body B1]\nsigma = 1\n" + tail, 6, "body with no volume");
    expect_error(head + "[Body B1]\nvolume = w\n" + tail, 6, "body with no sigma");

    // one volume per body
    expect_error(head + "[Body B1]\nvolume = w1 w2\nsigma = 1\n" + tail, 7, "two volumes in one body");

    // duplicate names and volumes
    const std::string b1 = "[Body B1]\nvolume = w\nsigma = 1\n";
    expect_error(head + b1 + "[Body B1]\nvolume = x\nsigma = 1\n" + tail, 9, "duplicate body name");
    expect_error(head + b1 + "[Body B2]\nvolume = w\nsigma = 1\n" + tail, 9, "volume claimed twice");
    expect_ok(head + b1 + "[Body B2]\nvolume = x\nsigma = 1\n" + tail, "control: distinct bodies");

    // length_unit
    expect_error("[mesh]\nfile = m.msh\nlength_unit = furlong\n[analysis]\ntype = dc\n" + b1 + tail, 3,
                 "unknown length unit");
}

void test_port_errors_and_controls() {
    const std::string head =
        "[mesh]\nfile = m.msh\nlength_unit = mm\n[analysis]\ntype = dc\n"
        "[Body B1]\nvolume = w\nsigma = 1\n";

    expect_error(head + "[port P1]\ntype = boundary_current\nsurface = s\nvoltage = 1\n", 12,
                 "voltage given to a current port");
    expect_error(head + "[port P1]\ntype = boundary_voltage\nsurface = s\ncurrent = 1\n", 12,
                 "current given to a voltage port");
    expect_error(head + "[port P1]\ntype = boundary_current\nsurface = s\n", 9,
                 "current port with no amplitude");
    expect_error(head + "[port P1]\ntype = wave_port\nsurface = s\ncurrent = 1\n", 10,
                 "unknown port type");
    expect_error(head + "[port]\ntype = boundary_voltage\nsurface = s\nvoltage = 0\n", 9,
                 "port with no name");

    // current_direction belongs to internal ports only.
    expect_error(head +
                     "[port P1]\ntype = boundary_voltage\nsurface = s\nvoltage = 0\n"
                     "current_direction = +y\n",
                 13, "current_direction on a boundary port");
    expect_ok(head + "[port P1]\ntype = internal_voltage\nsurface = s\nvoltage = 0\n"
                     "current_direction = +y\n",
              "control: current_direction on an internal port");

    // phase at DC
    expect_error(head + "[port P1]\ntype = boundary_voltage\nsurface = s\nvoltage = 0\nphase_deg = 90\n",
                 13, "non-zero phase at DC");
    expect_ok(head + "[port P1]\ntype = boundary_voltage\nsurface = s\nvoltage = 0\nphase_deg = 0\n",
              "control: phase_deg = 0 at DC is harmless");

    // shared surfaces and duplicate names
    const std::string p1 = "[port P1]\ntype = boundary_voltage\nsurface = s\nvoltage = 0\n";
    expect_error(head + p1 + "[port P1]\ntype = boundary_current\nsurface = t\ncurrent = 1\n", 13,
                 "duplicate port name");
    expect_error(head + p1 + "[port P2]\ntype = boundary_current\nsurface = s\ncurrent = 1\n", 13,
                 "two ports on one surface");
    expect_ok(head + p1 + "[port P2]\ntype = boundary_current\nsurface = t\ncurrent = 1\n",
              "control: distinct surfaces");
}

void test_whole_problem_errors() {
    const std::string mesh = "[mesh]\nfile = m.msh\nlength_unit = mm\n";
    const std::string dc = "[analysis]\ntype = dc\n";
    const std::string body = "[Body B1]\nvolume = w\nsigma = 1\n";

    expect_error(dc + body + "[port P1]\ntype = boundary_voltage\nsurface = s\nvoltage = 0\n", 0,
                 "no [mesh] section");
    expect_error(mesh + body + "[port P1]\ntype = boundary_voltage\nsurface = s\nvoltage = 0\n", 0,
                 "no [analysis] section");
    expect_error(mesh + dc + "[port P1]\ntype = boundary_voltage\nsurface = s\nvoltage = 0\n", 0,
                 "no bodies");
    expect_error(mesh + dc + body, 0, "no ports");

    // The singular-matrix case: nothing anywhere fixes Phi's constant.
    const std::string message = expect_error(
        mesh + dc + body + "[port P1]\ntype = boundary_current\nsurface = s\ncurrent = 1\n", 0,
        "every port is a boundary current source, so Phi has no reference");
    check(message.find("reference") != std::string::npos,
          "the no-reference message says what is missing");

    expect_ok(mesh + dc + body +
                  "[port P1]\ntype = boundary_current\nsurface = s\ncurrent = 1\n"
                  "[port P2]\ntype = boundary_voltage\nsurface = t\nvoltage = 0\n",
              "control: a 0 V port supplies the reference");
    expect_ok(mesh + dc + body +
                  "[port P1]\ntype = internal_current\nsurface = s\ncurrent = 1\n"
                  "current_direction = +y\n",
              "control: an internal port supplies the reference by grounding its minus side");
}

void test_current_direction_escalation() {
    const std::string head =
        "[mesh]\nfile = m.msh\nlength_unit = mm\n[analysis]\ntype = dc\n"
        "[Body B1]\nvolume = w\nsigma = 1\n";

    // One internal port without a hint: legal, but warned about.
    const ParseResult one = expect_ok(
        head + "[port P1]\ntype = internal_current\nsurface = s\ncurrent = 1\n",
        "one internal port without a hint is legal");
    check(one.warnings.size() == 1u, "one warning is produced");
    check(one.warnings.empty() || one.warnings[0].find("arbitrary") != std::string::npos,
          "the warning says the sign is arbitrary");
    check(one.warnings.empty() || one.warnings[0].find("current_direction") != std::string::npos,
          "the warning names the key that fixes it");

    // Two, and the relative sign becomes meaningless -- refused, not warned.
    expect_error(head +
                     "[port P1]\ntype = internal_current\nsurface = s\ncurrent = 1\n"
                     "[port P2]\ntype = internal_current\nsurface = t\ncurrent = 1\n",
                 9, "two internal ports with no hints");

    // With hints, two internal ports are fine and silent.
    const ParseResult two = expect_ok(
        head +
            "[port P1]\ntype = internal_current\nsurface = s\ncurrent = 1\ncurrent_direction = +y\n"
            "[port P2]\ntype = internal_current\nsurface = t\ncurrent = 1\ncurrent_direction = -z\n",
        "two internal ports with hints");
    check(two.warnings.empty(), "hinted internal ports produce no warning");
    check(two.problem.ports[0].current_direction.has_value(), "P1 hint survives");
    check(two.problem.ports[1].current_direction.has_value(), "P2 hint survives");
}

void test_direction_values() {
    const std::string head =
        "[mesh]\nfile = m.msh\nlength_unit = mm\n[analysis]\ntype = dc\n"
        "[Body B1]\nvolume = w\nsigma = 1\n"
        "[port P1]\ntype = internal_current\nsurface = s\ncurrent = 1\ncurrent_direction = ";

    const auto direction_of = [&](const std::string& value) {
        return *expect_ok(head + value + "\n", "direction '" + value + "'").problem.ports[0].current_direction;
    };

    const Vec3 py = direction_of("+y");
    check(near(py.x, 0.0) && near(py.y, 1.0) && near(py.z, 0.0), "+y");
    const Vec3 mz = direction_of("-z");
    check(near(mz.x, 0.0) && near(mz.y, 0.0) && near(mz.z, -1.0), "-z");
    const Vec3 x = direction_of("x");
    check(near(x.x, 1.0), "a bare axis name means the positive direction");
    const Vec3 vec = direction_of("0 1 0");
    check(near(vec.y, 1.0), "an explicit 3-vector");

    // Stored normalised, so magnitude is irrelevant -- the hint only ever
    // contributes a sign.
    const Vec3 big = direction_of("0 5 0");
    check(near(big.y, 1.0), "a long vector is normalised on read");
    const Vec3 skew = direction_of("0.1 0.9 -0.2");
    check(near(skew.norm(), 1.0), "a skew hint is normalised too");

    expect_error(head + "0 0 0\n", 13, "zero direction vector");
    expect_error(head + "+w\n", 13, "unknown axis letter");
    expect_error(head + "1 2\n", 13, "two-component direction");
    expect_error(head + "1 2 3 4\n", 13, "four-component direction");
}

void test_frequency_rules() {
    const std::string head = "[mesh]\nfile = m.msh\nlength_unit = mm\n";
    const std::string tail =
        "[Body B1]\nvolume = w\nsigma = 1\n"
        "[port P1]\ntype = boundary_voltage\nsurface = s\nvoltage = 0\n";

    expect_error(head + "[analysis]\ntype = dc\nfrequencies = 1e3\n" + tail, 6,
                 "frequencies at DC");
    expect_error(head + "[analysis]\ntype = dc\nf_start = 1\nf_stop = 2\nf_points = 3\n" + tail, 6,
                 "sweep keys at DC");
    expect_error(head + "[analysis]\ntype = dc\nformulation = reduced\n" + tail, 6,
                 "formulation at DC");
    expect_error(head + "[analysis]\ntype = frequency\n" + tail, 4,
                 "frequency analysis with no frequencies");
    expect_error(head + "[analysis]\ntype = frequency\nfrequencies = 1e3\nf_start = 1\n" + tail, 6,
                 "both frequency forms at once");
    expect_error(head + "[analysis]\ntype = frequency\nfrequencies = 0\n" + tail, 6,
                 "zero frequency");
    expect_error(head + "[analysis]\ntype = frequency\nfrequencies = -1e3\n" + tail, 6,
                 "negative frequency");
    expect_error(head + "[analysis]\ntype = frequency\nf_start = 1e3\nf_stop = 1e2\nf_points = 5\n" + tail,
                 7, "f_stop below f_start");
    expect_error(head + "[analysis]\ntype = frequency\nf_start = 1e3\nf_stop = 1e7\nf_points = 1\n" + tail,
                 8, "f_points of 1");
    expect_error(head + "[analysis]\ntype = frequency\nf_start = 1e3\nf_stop = 1e7\nf_points = 2.5\n" + tail,
                 8, "fractional f_points");
    expect_error(head + "[analysis]\ntype = frequency\nf_start = 1e3\nf_stop = 1e7\n" + tail, 4,
                 "sweep missing f_points");

    const ParseResult ok = expect_ok(
        head + "[analysis]\ntype = frequency\nsweep = log\nf_start = 1e3\nf_stop = 1e7\nf_points = 41\n" +
            tail,
        "control: a complete log sweep");
    check(ok.problem.frequencies.size() == 41u, "41 points");
}

void test_sweep_expansion() {
    using aphi_solver::expand_frequency_sweep;

    const std::vector<double> log_sweep = expand_frequency_sweep(1e3, 1e7, 41, true);
    check(log_sweep.size() == 41u, "41 points requested, 41 returned");

    // Endpoints must be exact, not merely close: a swept result has to line
    // up against a reference table without a tolerance on the frequency
    // itself.
    check(log_sweep.front() == 1e3, "log sweep starts exactly at f_start");
    check(log_sweep.back() == 1e7, "log sweep ends exactly at f_stop");

    // 4 decades over 40 intervals is 10 points per decade, so index 20 is
    // exactly two decades up.
    check(near(log_sweep[20], 1e5, 1e-6), "the log sweep midpoint is 1e5");
    check(near(log_sweep[10], 1e4, 1e-7), "one decade up is 1e4");

    for (std::size_t i = 1; i < log_sweep.size(); ++i) {
        check(log_sweep[i] > log_sweep[i - 1], "log sweep is strictly increasing at " + std::to_string(i));
        if (i > 1) break;  // one representative check, not 40 identical ones
    }

    const std::vector<double> lin = expand_frequency_sweep(0.0, 100.0, 5, false);
    check(lin.front() == 0.0 && lin.back() == 100.0, "linear sweep endpoints are exact");
    check(near(lin[1], 25.0) && near(lin[2], 50.0) && near(lin[3], 75.0), "linear sweep spacing");

    bool threw = false;
    try {
        expand_frequency_sweep(1e3, 1e7, 1, true);
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    check(threw, "fewer than 2 points is rejected");

    threw = false;
    try {
        expand_frequency_sweep(0.0, 1e7, 5, true);
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    check(threw, "a log sweep from zero is rejected");
}

void test_boundary_section() {
    const std::string head = "[mesh]\nfile = m.msh\nlength_unit = mm\n[analysis]\ntype = dc\n";
    const std::string tail =
        "[Body B1]\nvolume = w\nsigma = 1\n"
        "[port P1]\ntype = boundary_voltage\nsurface = s\nvoltage = 0\n";

    expect_ok(head + "[boundary]\nouter = flux_tangential\n" + tail, "flux_tangential is supported");
    expect_ok(head + tail, "control: [boundary] may be omitted entirely");

    // Recognised but not implemented -- rejected by name, so the message can
    // say "not yet" rather than "unknown".
    const std::string pec = expect_error(head + "[boundary]\nouter = pec\n" + tail, 7, "outer = pec");
    check(pec.find("not supported yet") != std::string::npos,
          "pec is refused as unimplemented, not as unknown");
    expect_error(head + "[boundary]\nouter = abc\n" + tail, 7, "outer = abc");
    expect_error(head + "[boundary]\nouter = nonsense\n" + tail, 7, "unknown outer value");

    // Per-surface boundary conditions are reserved, not repurposed.
    expect_error(head + "[boundary shield]\nouter = pec\n" + tail, 6, "named [boundary] section");
}

}  // namespace

int main() {
    test_cylinder_example_parses();
    test_frequency_example_parses();
    test_lexing();
    test_syntax_errors();
    test_value_errors_and_controls();
    test_port_errors_and_controls();
    test_whole_problem_errors();
    test_current_direction_escalation();
    test_direction_values();
    test_frequency_rules();
    test_sweep_expansion();
    test_boundary_section();

    std::cout << (g_checks - g_failures) << "/" << g_checks << " checks passed\n";
    return g_failures == 0 ? 0 : 1;
}
