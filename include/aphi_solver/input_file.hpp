#pragma once

#include <stdexcept>
#include <string>
#include <vector>

#include "aphi_solver/problem.hpp"

namespace aphi_solver {

/// Thrown for any malformed input file. Carries the file and line so the
/// message can point at the offending text rather than merely describing
/// it -- an error that names the wrong line is nearly as unhelpful as no
/// error, so `line` is threaded through every check that has one.
///
/// `line` is 0 for a problem with the file as a whole (it does not exist,
/// it has no ports at all) rather than with one line of it.
class InputError : public std::runtime_error {
public:
    InputError(std::string path, int line, const std::string& message);

    const std::string& path() const { return path_; }
    int line() const { return line_; }

    /// Just the explanation, without the "file:line: " prefix that what()
    /// carries. Useful to tests that check the wording.
    const std::string& message() const { return message_; }

private:
    std::string path_;
    int line_ = 0;
    std::string message_;
};

/// A parsed problem plus anything worth saying about it that is not fatal.
///
/// Warnings exist because some things are legal, deliberate, and still
/// worth stating once -- chiefly an internal port with no
/// `current_direction`, whose reported I and V then carry an arbitrary
/// sign (see problem.hpp's orientation note). Silence there would be
/// misleading; an error would be wrong, since the choice is a real one.
struct ParseResult {
    Problem problem;
    std::vector<std::string> warnings;
};

/// Reads and validates an input file. See
/// `docs/INPUT_FILE_PLAN.md` Sec. 1 for the format and Sec. 3.1
/// for the complete list of checks.
///
/// **This stage never touches the mesh.** Every check here is answerable
/// from the text alone: syntax, required keys, value ranges, and the
/// cross-section rules (a potential reference exists; no two ports share a
/// surface). Whether a name like `wire_bottom` exists in the mesh, or
/// whether its faces really are on the boundary, is the binding stage's
/// job -- those are a different class of mistake and reporting them
/// separately is what lets each name the right thing.
///
/// A relative `[mesh] file` is resolved against the directory containing
/// `path`, so an input file can sit beside its mesh and be run from
/// anywhere.
///
/// Throws InputError on any failure, including a file that cannot be read.
ParseResult parse_input_file(const std::string& path);

/// Same parser, reading from memory. This is the form the tests use --
/// every fixture is then a string literal beside its assertions rather
/// than a temp file, which keeps a test readable as one piece.
///
/// `path_label` appears in error messages only. `base_dir` resolves a
/// relative `[mesh] file`; empty (the default) leaves the path exactly as
/// written, which is what a test wants.
ParseResult parse_input_string(const std::string& text, const std::string& path_label = "<input>",
                               const std::string& base_dir = "");

/// Expands a frequency sweep into the explicit list the `Problem` carries.
/// Exposed for its own tests: it is the only arithmetic in this stage, and
/// the endpoints must come out exactly (a 41-point log sweep from 1e3 to
/// 1e7 has to end at 1e7, not 0.99999e7), which is easy to get subtly
/// wrong and invisible afterwards.
///
/// Throws std::invalid_argument if `points < 2`, if either endpoint is
/// non-positive under `log` spacing, or if `stop <= start`.
std::vector<double> expand_frequency_sweep(double start, double stop, int points, bool logarithmic);

}  // namespace aphi_solver
