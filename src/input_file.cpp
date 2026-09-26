#include "aphi_solver/input_file.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <map>
#include <set>
#include <sstream>

namespace aphi_solver {

InputError::InputError(std::string path, int line, const std::string& message)
    : std::runtime_error(path + (line > 0 ? ":" + std::to_string(line) : "") + ": " + message),
      path_(std::move(path)),
      line_(line),
      message_(message) {}

std::vector<double> expand_frequency_sweep(double start, double stop, int points, bool logarithmic) {
    if (points < 2) throw std::invalid_argument("expand_frequency_sweep: points must be >= 2");
    if (stop <= start) throw std::invalid_argument("expand_frequency_sweep: stop must exceed start");
    if (logarithmic && start <= 0.0) {
        throw std::invalid_argument("expand_frequency_sweep: log spacing needs a positive start");
    }

    std::vector<double> out(static_cast<std::size_t>(points));

    // The endpoints are written in directly rather than evaluated, so they
    // come out exact. Evaluating them alongside the interior points would
    // leave 1e7 as 9.999999999999998e6 -- harmless numerically, but it
    // makes a swept result impossible to line up against a reference table
    // by eye or by exact comparison.
    out.front() = start;
    out.back() = stop;

    const double n = static_cast<double>(points - 1);
    if (logarithmic) {
        const double log_start = std::log(start);
        const double step = (std::log(stop) - log_start) / n;
        for (int i = 1; i < points - 1; ++i) {
            out[static_cast<std::size_t>(i)] = std::exp(log_start + step * i);
        }
    } else {
        const double step = (stop - start) / n;
        for (int i = 1; i < points - 1; ++i) {
            out[static_cast<std::size_t>(i)] = start + step * i;
        }
    }
    return out;
}

namespace {

// ---------------------------------------------------------------------------
// Lexing: text -> sections of key/value entries, each carrying its line.

struct Entry {
    std::string key;    // lowercased
    std::string value;  // trimmed, comment stripped, case preserved
    int line = 0;
};

struct Section {
    std::string kind;  // lowercased: mesh, analysis, boundary, body, port, ...
    std::string name;  // as written; empty for the unnamed sections
    int line = 0;
    std::vector<Entry> entries;
};

std::string lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return s;
}

std::string trim(const std::string& s) {
    const auto is_space = [](unsigned char c) { return std::isspace(c) != 0; };
    std::size_t b = 0;
    while (b < s.size() && is_space(static_cast<unsigned char>(s[b]))) ++b;
    std::size_t e = s.size();
    while (e > b && is_space(static_cast<unsigned char>(s[e - 1]))) --e;
    return s.substr(b, e - b);
}

std::vector<std::string> split_tokens(const std::string& s) {
    std::vector<std::string> out;
    std::istringstream in(s);
    std::string token;
    while (in >> token) out.push_back(token);
    return out;
}

std::string join(const std::vector<std::string>& items, const char* sep) {
    std::string out;
    for (std::size_t i = 0; i < items.size(); ++i) {
        if (i) out += sep;
        out += items[i];
    }
    return out;
}

class Parser {
public:
    Parser(const std::string& text, std::string path, std::string base_dir)
        : text_(text), path_(std::move(path)), base_dir_(std::move(base_dir)) {}

    ParseResult run();

private:
    [[noreturn]] void fail(int line, const std::string& message) const {
        throw InputError(path_, line, message);
    }

    // -- lexing --
    std::vector<Section> lex() const;

    // -- per-section handlers --
    void read_mesh(const Section& s);
    void read_analysis(const Section& s);
    void read_boundary(const Section& s);
    void read_solver(const Section& s);
    void read_body(const Section& s);
    void read_port(const Section& s);

    // -- cross-section checks --
    void check_whole_problem(const std::vector<Section>& sections);

    // -- value accessors, all reporting the entry's own line --
    const Entry* find(const Section& s, const std::string& key) const;
    const Entry& require(const Section& s, const std::string& key) const;
    void reject_unknown_keys(const Section& s, const std::vector<std::string>& allowed) const;

    double as_number(const Entry& e) const;
    int as_int(const Entry& e) const;
    std::vector<double> as_numbers(const Entry& e) const;
    std::string as_keyword(const Entry& e, const std::vector<std::string>& allowed) const;
    std::vector<std::string> as_names(const Entry& e) const;
    Vec3 as_direction(const Entry& e) const;

    const std::string& text_;
    std::string path_;
    std::string base_dir_;

    ParseResult result_;
    bool saw_mesh_ = false;
    bool saw_analysis_ = false;
    int analysis_line_ = 0;
    int conditioning_line_ = 0;  ///< where [solver] conditioning was set, for the DC check
};

// ---------------------------------------------------------------------------

std::vector<Section> Parser::lex() const {
    std::vector<Section> sections;
    std::istringstream in(text_);
    std::string raw;
    int line_number = 0;

    while (std::getline(in, raw)) {
        ++line_number;
        if (!raw.empty() && raw.back() == '\r') raw.pop_back();  // tolerate CRLF

        // Strip a UTF-8 byte-order mark. Notepad, Visual Studio and
        // PowerShell's Set-Content all write one by default on Windows, and
        // without this the first line fails to parse with a message about an
        // invisible character -- which is about as unhelpful as an error can
        // be. Found by running the driver on a file PowerShell had written.
        if (line_number == 1 && raw.size() >= 3 && static_cast<unsigned char>(raw[0]) == 0xEF &&
            static_cast<unsigned char>(raw[1]) == 0xBB && static_cast<unsigned char>(raw[2]) == 0xBF) {
            raw.erase(0, 3);
        }

        const std::size_t hash = raw.find('#');
        const std::string line = trim(hash == std::string::npos ? raw : raw.substr(0, hash));
        if (line.empty()) continue;

        if (line.front() == '[') {
            if (line.back() != ']') {
                fail(line_number, "section header is missing its closing ']'");
            }
            const std::vector<std::string> parts = split_tokens(line.substr(1, line.size() - 2));
            if (parts.empty()) fail(line_number, "empty section header '[]'");
            if (parts.size() > 2) {
                fail(line_number, "a section header is '[kind]' or '[kind name]', but this has " +
                                      std::to_string(parts.size()) + " words");
            }
            Section s;
            s.kind = lower(parts[0]);
            s.name = parts.size() == 2 ? parts[1] : std::string();
            s.line = line_number;
            sections.push_back(std::move(s));
            continue;
        }

        const std::size_t eq = line.find('=');
        if (eq == std::string::npos) {
            fail(line_number, "expected 'key = value' or a '[section]' header, but found '" + line + "'");
        }
        if (sections.empty()) {
            fail(line_number, "'" + trim(line.substr(0, eq)) + "' appears before any [section] header");
        }

        Entry e;
        e.key = lower(trim(line.substr(0, eq)));
        e.value = trim(line.substr(eq + 1));
        e.line = line_number;
        if (e.key.empty()) fail(line_number, "missing key before '='");
        if (e.value.empty()) fail(line_number, "'" + e.key + "' has no value");

        Section& current = sections.back();
        for (const Entry& existing : current.entries) {
            if (existing.key == e.key) {
                fail(line_number, "'" + e.key + "' is set twice in this section (first at line " +
                                      std::to_string(existing.line) + ")");
            }
        }
        current.entries.push_back(std::move(e));
    }
    return sections;
}

// ---------------------------------------------------------------------------
// Value accessors.

const Entry* Parser::find(const Section& s, const std::string& key) const {
    for (const Entry& e : s.entries) {
        if (e.key == key) return &e;
    }
    return nullptr;
}

const Entry& Parser::require(const Section& s, const std::string& key) const {
    if (const Entry* e = find(s, key)) return *e;
    const std::string what = s.name.empty() ? "[" + s.kind + "]" : "[" + s.kind + " " + s.name + "]";
    fail(s.line, what + " is missing the required key '" + key + "'");
}

void Parser::reject_unknown_keys(const Section& s, const std::vector<std::string>& allowed) const {
    for (const Entry& e : s.entries) {
        if (std::find(allowed.begin(), allowed.end(), e.key) == allowed.end()) {
            fail(e.line, "unknown key '" + e.key + "' in [" + s.kind + "]; valid keys are " +
                             join(allowed, ", "));
        }
    }
}

double Parser::as_number(const Entry& e) const {
    const std::vector<std::string> tokens = split_tokens(e.value);
    if (tokens.size() != 1) {
        fail(e.line, "'" + e.key + "' takes one number, but found " + std::to_string(tokens.size()) +
                         " values");
    }
    try {
        std::size_t used = 0;
        const double v = std::stod(tokens[0], &used);
        if (used != tokens[0].size()) throw std::invalid_argument("trailing");
        return v;
    } catch (const std::exception&) {
        fail(e.line, "'" + e.key + "' expects a number, but found '" + tokens[0] + "'");
    }
}

int Parser::as_int(const Entry& e) const {
    const double v = as_number(e);
    if (v != std::floor(v)) {
        fail(e.line, "'" + e.key + "' expects a whole number, but found '" + e.value + "'");
    }
    return static_cast<int>(v);
}

std::vector<double> Parser::as_numbers(const Entry& e) const {
    const std::vector<std::string> tokens = split_tokens(e.value);
    std::vector<double> out;
    out.reserve(tokens.size());
    for (const std::string& t : tokens) {
        try {
            std::size_t used = 0;
            const double v = std::stod(t, &used);
            if (used != t.size()) throw std::invalid_argument("trailing");
            out.push_back(v);
        } catch (const std::exception&) {
            fail(e.line, "'" + e.key + "' expects numbers, but found '" + t + "'");
        }
    }
    return out;
}

std::string Parser::as_keyword(const Entry& e, const std::vector<std::string>& allowed) const {
    const std::string v = lower(e.value);
    if (std::find(allowed.begin(), allowed.end(), v) == allowed.end()) {
        fail(e.line, "'" + e.key + "' must be one of " + join(allowed, ", ") + ", but found '" +
                         e.value + "'");
    }
    return v;
}

std::vector<std::string> Parser::as_names(const Entry& e) const {
    // Names stay case-sensitive: they have to match Gmsh physical names
    // exactly, and Gmsh does not fold case.
    return split_tokens(e.value);
}

Vec3 Parser::as_direction(const Entry& e) const {
    const std::vector<std::string> tokens = split_tokens(e.value);

    if (tokens.size() == 1) {
        // Axis shorthand: +x, -y, z. Sugar for the obvious vector; the hint
        // only contributes a sign, so this covers most real cuts.
        std::string t = lower(tokens[0]);
        double sign = 1.0;
        if (t.size() == 2 && (t[0] == '+' || t[0] == '-')) {
            sign = t[0] == '-' ? -1.0 : 1.0;
            t = t.substr(1);
        }
        if (t == "x") return Vec3(sign, 0.0, 0.0);
        if (t == "y") return Vec3(0.0, sign, 0.0);
        if (t == "z") return Vec3(0.0, 0.0, sign);
        fail(e.line, "'" + e.key + "' expects an axis (+x, -y, z, ...) or three numbers, but found '" +
                         e.value + "'");
    }

    if (tokens.size() != 3) {
        fail(e.line, "'" + e.key + "' expects an axis (+x, -y, z, ...) or three numbers, but found " +
                         std::to_string(tokens.size()) + " values");
    }
    const std::vector<double> v = as_numbers(e);
    const Vec3 d(v[0], v[1], v[2]);
    if (d.norm() == 0.0) {
        fail(e.line, "'" + e.key + "' is the zero vector, which points nowhere");
    }
    return d * (1.0 / d.norm());
}

// ---------------------------------------------------------------------------
// Section handlers.

void Parser::read_mesh(const Section& s) {
    reject_unknown_keys(s, {"file", "length_unit"});

    // The path is taken verbatim rather than tokenised: a directory name may
    // contain spaces, and splitting one would be a baffling failure.
    const Entry& file = require(s, "file");
    std::string path = file.value;
    if (!base_dir_.empty()) {
        const bool absolute = path.size() > 1 && (path[0] == '/' || path[0] == '\\' || path[1] == ':');
        if (!absolute) path = base_dir_ + "/" + path;
    }
    result_.problem.mesh_file = path;

    const Entry& unit = require(s, "length_unit");
    const std::string u = as_keyword(unit, {"m", "mm", "um", "nm"});
    if (u == "m") result_.problem.length_unit = LengthUnit::Metre;
    else if (u == "mm") result_.problem.length_unit = LengthUnit::Millimetre;
    else if (u == "um") result_.problem.length_unit = LengthUnit::Micrometre;
    else result_.problem.length_unit = LengthUnit::Nanometre;
}

void Parser::read_analysis(const Section& s) {
    reject_unknown_keys(s, {"type", "frequencies", "sweep", "f_start", "f_stop", "f_points", "formulation"});
    analysis_line_ = s.line;

    Problem& p = result_.problem;
    const std::string type = as_keyword(require(s, "type"), {"dc", "frequency"});
    p.type = type == "dc" ? AnalysisType::DC : AnalysisType::Frequency;

    const Entry* explicit_list = find(s, "frequencies");
    const Entry* sweep = find(s, "sweep");
    const Entry* f_start = find(s, "f_start");
    const Entry* f_stop = find(s, "f_stop");
    const Entry* f_points = find(s, "f_points");
    const bool any_sweep_key = sweep || f_start || f_stop || f_points;

    if (p.type == AnalysisType::DC) {
        // At DC these keys do not describe a choice the user has, so a file
        // carrying one says something its author cannot have meant.
        if (explicit_list) fail(explicit_list->line, "'frequencies' does not apply when type = dc");
        for (const Entry* e : {sweep, f_start, f_stop, f_points}) {
            if (e) fail(e->line, "'" + e->key + "' does not apply when type = dc");
        }
        if (const Entry* f = find(s, "formulation")) {
            fail(f->line,
                 "'formulation' does not apply when type = dc: at omega = 0 the Phi equation is "
                 "empty wherever sigma = 0, so Phi lives on the conductors whatever is written here");
        }
        return;
    }

    if (const Entry* f = find(s, "formulation")) {
        p.formulation =
            as_keyword(*f, {"full_wave", "reduced"}) == "reduced" ? Formulation::Reduced : Formulation::FullWave;
    }

    if (explicit_list && any_sweep_key) {
        fail(explicit_list->line,
             "give either 'frequencies' or a sweep (sweep / f_start / f_stop / f_points), not both");
    }
    if (!explicit_list && !any_sweep_key) {
        fail(s.line, "type = frequency needs either 'frequencies' or a sweep "
                     "(f_start, f_stop, f_points)");
    }

    if (explicit_list) {
        p.frequencies = as_numbers(*explicit_list);
        if (p.frequencies.empty()) fail(explicit_list->line, "'frequencies' is empty");
        for (double f : p.frequencies) {
            if (f <= 0.0) {
                fail(explicit_list->line,
                     "every frequency must be positive, but found " + std::to_string(f) +
                         " (use type = dc for zero)");
            }
        }
        return;
    }

    if (!f_start) fail(s.line, "a sweep needs 'f_start'");
    if (!f_stop) fail(s.line, "a sweep needs 'f_stop'");
    if (!f_points) fail(s.line, "a sweep needs 'f_points'");

    const bool logarithmic = !sweep || as_keyword(*sweep, {"log", "linear"}) == "log";
    const double start = as_number(*f_start);
    const double stop = as_number(*f_stop);
    const int points = as_int(*f_points);

    if (start <= 0.0) fail(f_start->line, "'f_start' must be positive");
    if (stop <= start) fail(f_stop->line, "'f_stop' must be greater than 'f_start'");
    if (points < 2) fail(f_points->line, "'f_points' must be at least 2 (it counts both endpoints)");

    p.frequencies = expand_frequency_sweep(start, stop, points, logarithmic);
}

void Parser::read_boundary(const Section& s) {
    reject_unknown_keys(s, {"outer"});
    if (const Entry* outer = find(s, "outer")) {
        const std::string v = as_keyword(*outer, {"flux_tangential", "pec", "abc"});
        if (v != "flux_tangential") {
            fail(outer->line, "'outer = " + v + "' is recognised but not supported yet; only "
                              "flux_tangential (n x A = 0 with Phi free) is implemented");
        }
        result_.problem.outer = OuterBoundary::FluxTangential;
    }
}

void Parser::read_solver(const Section& s) {
    reject_unknown_keys(s, {"conditioning"});
    if (const Entry* e = find(s, "conditioning")) {
        const std::string word = as_keyword(*e, {"natural", "row_scaled", "scaled_phi"});
        Conditioning picked = Conditioning::Natural;
        if (!conditioning_from_keyword(word, picked)) {
            fail(e->line, "'conditioning = " + word + "' is not one I know");
        }
        result_.problem.conditioning = picked;
        conditioning_line_ = e->line;
    }
}

void Parser::read_body(const Section& s) {
    if (s.name.empty()) fail(s.line, "a [Body] section needs a name, as in '[Body B1]'");
    reject_unknown_keys(s, {"volume", "sigma", "eps_r", "mu_r"});

    Body b;
    b.name = s.name;
    b.line = s.line;

    const std::vector<std::string> volumes = as_names(require(s, "volume"));
    if (volumes.size() != 1) {
        fail(require(s, "volume").line,
             "'volume' names exactly one Physical Volume; for several volumes of one material, "
             "write one [Body ...] section each");
    }
    b.volume = volumes[0];

    const Entry& sigma = require(s, "sigma");
    b.sigma = as_number(sigma);
    if (b.sigma < 0.0) fail(sigma.line, "'sigma' cannot be negative");

    if (const Entry* e = find(s, "eps_r")) {
        b.eps_r = as_number(*e);
        if (b.eps_r <= 0.0) fail(e->line, "'eps_r' must be positive");
    }
    if (const Entry* e = find(s, "mu_r")) {
        b.mu_r = as_number(*e);
        if (b.mu_r <= 0.0) fail(e->line, "'mu_r' must be positive");
    }

    for (const Body& existing : result_.problem.bodies) {
        if (existing.name == b.name) fail(s.line, "a body named '" + b.name + "' is already defined");
        if (existing.volume == b.volume) {
            fail(s.line, "volume '" + b.volume + "' is already claimed by body '" + existing.name + "'");
        }
    }
    result_.problem.bodies.push_back(std::move(b));
}

void Parser::read_port(const Section& s) {
    if (s.name.empty()) fail(s.line, "a [port] section needs a name, as in '[port P1]'");
    reject_unknown_keys(s, {"type", "surface", "current", "voltage", "phase_deg", "current_direction"});

    Port p;
    p.name = s.name;
    p.line = s.line;

    const std::string type = as_keyword(
        require(s, "type"),
        {"boundary_current", "boundary_voltage", "internal_current", "internal_voltage"});
    if (type == "boundary_current") p.type = PortType::BoundaryCurrent;
    else if (type == "boundary_voltage") p.type = PortType::BoundaryVoltage;
    else if (type == "internal_current") p.type = PortType::InternalCurrent;
    else p.type = PortType::InternalVoltage;

    p.surface = as_names(require(s, "surface"));
    if (p.surface.empty()) fail(require(s, "surface").line, "'surface' names no surface");

    // The amplitude key has to match the drive type. Supplying the other one
    // almost always means the port type itself is wrong, so say that rather
    // than silently ignoring it.
    const Entry* current = find(s, "current");
    const Entry* voltage = find(s, "voltage");
    const bool wants_current = is_current_driven(p.type);
    const Entry* wanted = wants_current ? current : voltage;
    const Entry* unwanted = wants_current ? voltage : current;
    const char* wanted_key = wants_current ? "current" : "voltage";
    const char* unwanted_key = wants_current ? "voltage" : "current";

    if (unwanted) {
        fail(unwanted->line, std::string("'") + unwanted_key + "' does not belong on a " + type +
                                 " port; it takes '" + wanted_key + "'");
    }
    if (!wanted) {
        fail(s.line, "[port " + s.name + "] is missing '" + wanted_key + "'");
    }
    p.amplitude = as_number(*wanted);

    if (const Entry* phase = find(s, "phase_deg")) {
        p.phase_deg = as_number(*phase);
        if (result_.problem.type == AnalysisType::DC && p.phase_deg != 0.0) {
            fail(phase->line, "a non-zero 'phase_deg' has no meaning when type = dc");
        }
    }

    if (const Entry* dir = find(s, "current_direction")) {
        if (!is_internal(p.type)) {
            fail(dir->line,
                 "'current_direction' applies only to internal ports; a boundary port's positive "
                 "direction is into the domain, which the mesh already fixes");
        }
        p.current_direction = as_direction(*dir);
    }

    for (const Port& existing : result_.problem.ports) {
        if (existing.name == p.name) fail(s.line, "a port named '" + p.name + "' is already defined");
        for (const std::string& a : existing.surface) {
            for (const std::string& b : p.surface) {
                if (a == b) {
                    fail(s.line, "surface '" + a + "' is already used by port '" + existing.name + "'");
                }
            }
        }
    }
    result_.problem.ports.push_back(std::move(p));
}

// ---------------------------------------------------------------------------

void Parser::check_whole_problem(const std::vector<Section>& sections) {
    const Problem& p = result_.problem;

    if (!saw_mesh_) fail(0, "no [mesh] section");
    if (!saw_analysis_) fail(0, "no [analysis] section");
    if (p.bodies.empty()) fail(0, "no [Body ...] section: the model has no material anywhere");
    if (p.ports.empty()) fail(0, "no [port ...] section: nothing would drive the problem");

    // Both symmetric forms divide by j*omega -- RowScaled in the Phi rows,
    // ScaledPhi in any prescribed potential -- so neither exists at DC. That
    // is caught here, with the line the user wrote, rather than left to throw
    // from inside assembly with no idea where it came from.
    if (p.type == AnalysisType::DC && conditioning_needs_ac(p.conditioning)) {
        fail(conditioning_line_,
             std::string("conditioning = ") + conditioning_keyword(p.conditioning) +
                 " needs a non-zero frequency: it scales the Phi equations by 1/(j*omega), which "
                 "is undefined at DC. At DC the system decouples and there is no coupling block "
                 "to symmetrise, so use conditioning = natural");
    }

    // A potential reference is needed or Phi is fixed only up to a constant
    // and the matrix is singular. A voltage port supplies one directly; an
    // internal port supplies one too, since it holds its minus side at 0 V.
    bool has_reference = false;
    for (const Port& port : p.ports) {
        if (!is_current_driven(port.type) || is_internal(port.type)) has_reference = true;
    }
    if (!has_reference) {
        fail(0, "no potential reference: every port is a boundary current source, so Phi is "
                "determined only up to a constant. Add a boundary_voltage port with voltage = 0");
    }

    // An internal port with no hint gets a direction derived from the mesh,
    // which fixes the sign of its reported I and V arbitrarily. Tolerable
    // for one port -- Z, R, L and P are unaffected -- but between two ports
    // the *relative* sign is what a coupling term is made of, so there it is
    // refused rather than warned about.
    std::vector<std::string> unhinted;
    int internal_ports = 0;
    for (const Port& port : p.ports) {
        if (!is_internal(port.type)) continue;
        ++internal_ports;
        if (!port.current_direction.has_value()) unhinted.push_back(port.name);
    }
    if (!unhinted.empty()) {
        if (internal_ports >= 2) {
            // Point at the first offending section so the message has a line.
            int line = 0;
            for (const Section& s : sections) {
                if (s.kind == "port" && s.name == unhinted.front()) line = s.line;
            }
            fail(line, "with more than one internal port, every one needs 'current_direction': "
                       "without it the signs of their currents are unrelated, and any coupling "
                       "between them would be meaningless. Missing on: " + join(unhinted, ", "));
        }
        result_.warnings.push_back(
            "port '" + unhinted.front() +
            "' has no current_direction, so the sign of its reported I and V is arbitrary "
            "(Z, R, L and P are not). Set current_direction to fix it.");
    }
}

ParseResult Parser::run() {
    const std::vector<Section> sections = lex();

    // [analysis] is read before the ports, which need to know whether this
    // is a DC problem to judge 'phase_deg'.
    for (const Section& s : sections) {
        if (s.kind == "analysis") {
            if (saw_analysis_) fail(s.line, "a second [analysis] section");
            read_analysis(s);
            saw_analysis_ = true;
        }
    }

    for (const Section& s : sections) {
        if (s.kind == "mesh") {
            if (saw_mesh_) fail(s.line, "a second [mesh] section");
            if (!s.name.empty()) fail(s.line, "[mesh] takes no name");
            read_mesh(s);
            saw_mesh_ = true;
        } else if (s.kind == "analysis") {
            if (!s.name.empty()) fail(s.line, "[analysis] takes no name");
        } else if (s.kind == "boundary") {
            if (!s.name.empty()) {
                fail(s.line, "per-surface boundary conditions are not supported yet; [boundary] "
                             "takes no name and applies to the whole outer boundary");
            }
            read_boundary(s);
        } else if (s.kind == "body") {
            read_body(s);
        } else if (s.kind == "port") {
            read_port(s);
        } else if (s.kind == "solver") {
            read_solver(s);
        } else if (s.kind == "output") {
            fail(s.line, "[" + s.kind + "] is reserved but not supported yet");
        } else {
            fail(s.line, "unknown section '[" + s.kind + "]'; valid sections are mesh, analysis, "
                         "boundary, solver, Body, port");
        }
    }

    check_whole_problem(sections);
    return std::move(result_);
}

}  // namespace

ParseResult parse_input_string(const std::string& text, const std::string& path_label,
                               const std::string& base_dir) {
    Parser parser(text, path_label, base_dir);
    return parser.run();
}

ParseResult parse_input_file(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw InputError(path, 0, "cannot open file");

    std::ostringstream buffer;
    buffer << in.rdbuf();

    std::string dir;
    const std::size_t slash = path.find_last_of("/\\");
    if (slash != std::string::npos) dir = path.substr(0, slash);

    return parse_input_string(buffer.str(), path, dir);
}

}  // namespace aphi_solver
