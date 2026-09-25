#pragma once

#include <complex>
#include <optional>
#include <string>
#include <vector>

#include "aphi_solver/mesh.hpp"  // Vec3

namespace aphi_solver {

/// The in-memory description of one problem: everything the input file
/// says, parsed and range-checked, but **not yet bound to a mesh**. Names
/// here are still strings -- resolving them to physical groups, faces,
/// nodes and edges is the binding step's job (see
/// `docs/INPUT_FILE_PLAN.md` Sec. 3.3), and it is deliberately
/// separate: a name that does not exist in the mesh is a different class of
/// error from a malformed file, and reporting them at different stages is
/// what lets each one name the right thing.
///
/// This struct, not the file format, is the interface the rest of the
/// solver is written against. The parser fills it; binding consumes it; the
/// DOF map consumes what binding produces. Sequencing it first means the
/// later steps are written against something settled.

/// Mesh files carry no units -- Gmsh writes bare numbers -- so the input
/// file must say what they mean.
enum class LengthUnit { Metre, Millimetre, Micrometre, Nanometre };

/// Metres per unit, the factor mesh coordinates are multiplied by at bind
/// time so everything downstream is SI.
double length_scale(LengthUnit unit);

enum class AnalysisType {
    DC,        ///< omega = 0. The system decouples (FORMULATION.md Sec. 2).
    Frequency  ///< one solve per entry of Problem::frequencies.
};

/// Which form of the A-Phi system to assemble. Only meaningful at
/// AnalysisType::Frequency -- see phi_on_conductors_only below for the one
/// place this actually changes.
enum class Formulation {
    FullWave,  ///< eps retained; the all-frequency system, Phi everywhere.
    Reduced    ///< eps -> 0; the eddy-current variant, Phi on conductors.
};

/// Treatment of the outer domain boundary. Only one value today; the enum
/// exists so adding PEC and ABC later is a new enumerator rather than a
/// change of representation.
enum class OuterBoundary {
    FluxTangential  ///< n x A = 0 with Phi left free. See ROADMAP Phase 04.
};

enum class PortType {
    BoundaryCurrent,   ///< terminal on the domain boundary, current driven
    BoundaryVoltage,   ///< terminal on the domain boundary, voltage driven
    InternalCurrent,   ///< cut through a conductor, current driven
    InternalVoltage    ///< cut through a conductor, voltage driven
};

/// **Port orientation convention** (decided 24 Sept 2026). Every internal
/// cut has a direction `d`, and it fixes everything:
///
///  - the **minus** side is the one `d` emits from; it is **held at 0 V**;
///  - the **plus** side is the one `d` points into; it carries the port's
///    single unknown potential;
///  - `I` is the current crossing the cut along `d`, and
///    `V = Phi(plus) - Phi(minus) = Phi(plus)`.
///
/// So `P = 0.5*Re(V*conj(I))` is the power the source delivers: current runs
/// minus -> plus inside the gap, then round the conductor from plus back to
/// minus, which is the battery convention.
///
/// There is deliberately **no sign factor anywhere**. An earlier design let
/// the input file ground either side, which needed a per-port `s = +/-1`
/// multiplying the prescribed potential and the port row. Grounding the
/// minus side always makes `s` identically +1, so it was removed rather
/// than left as an untestable branch: the only thing the choice ever
/// controlled was Phi's arbitrary offset, never V, I, Z or P.
///
/// `d` itself is **not in the input file** and is not stored here -- it is
/// derived from the port's own faces at bind time, because the file has no
/// good way to let a user point at a direction. The consequence, which must
/// be reported with any internal-port result: the *sign* of that port's I
/// and V is arbitrary (mesh-determined, though deterministic for a given
/// mesh). Z, R, L, |I| and P are unaffected, since I and V flip together.
/// Letting the user choose `d` is the open question for a future interface.

/// One body: exactly one mesh volume plus its material. One `[Body ...]`
/// section of the input file.
struct Body {
    std::string name;    ///< the section's own handle, e.g. "B1"
    int line = 0;        ///< its line in the input file; 0 if built in code
    std::string volume;  ///< Physical Volume name, or a tag written as digits
    double sigma = 0.0;  ///< S/m. Required in the file -- never defaulted, so
                         ///< a conductor cannot silently become an insulator.
    double eps_r = 1.0;
    double mu_r = 1.0;
};

/// One port. `amplitude` holds whichever of `current` / `voltage` the type
/// calls for -- a port has exactly one excitation, so storing one field
/// keeps "which key applies" a property of `type` alone rather than a pair
/// of optional fields that could both be set.
///
/// Note what is *not* here: no ground side. Which side of a cut sits at 0 V
/// follows from the direction (see the orientation note above), so there is
/// nothing to store.
struct Port {
    std::string name;
    int line = 0;        ///< its line in the input file; 0 if built in code
    PortType type = PortType::BoundaryCurrent;
    std::vector<std::string> surface;  ///< Physical Surface names or tags

    double amplitude = 0.0;  ///< magnitude: amperes or volts, peak
    double phase_deg = 0.0;  ///< degrees; positive leads, under e^{+j*omega*t}

    /// Internal ports only, optional. **A hint, not the normal.**
    ///
    /// The cut's normal is computed from its own faces; all this has to do
    /// is resolve the remaining +/- sign, via
    /// `d = sign(hint . n) * n` (resolve_cut_direction below). So any vector
    /// within 90 degrees of the intended current direction gives the same
    /// answer, and a user who knows only that current should flow "roughly
    /// +y" can say exactly that. They never have to compute a normal --
    /// which is the whole point, since a mesh normal is not something
    /// anyone can write down by hand.
    ///
    /// Absent means "derive it": `canonical_orientation` picks a
    /// deterministic sign and the port's reported I and V then carry an
    /// arbitrary (though mesh-stable) sign. That is a warning for one
    /// internal port and an error for two or more, where *relative* signs
    /// between ports would otherwise be meaningless.
    std::optional<Vec3> current_direction;
};

struct Problem {
    std::string mesh_file;  ///< as written, relative to the input file
    LengthUnit length_unit = LengthUnit::Metre;

    AnalysisType type = AnalysisType::DC;

    /// One entry per solve; empty at DC. A `sweep` in the input file is
    /// expanded into this list by the parser, so nothing downstream sees a
    /// second way of saying the same thing.
    std::vector<double> frequencies;

    Formulation formulation = Formulation::FullWave;
    OuterBoundary outer = OuterBoundary::FluxTangential;

    std::vector<Body> bodies;
    std::vector<Port> ports;
};

/// True for the two cut-based port types.
bool is_internal(PortType type);

/// True for the two current-driven types. A current port prescribes I and
/// solves for V; a voltage port does the reverse.
bool is_current_driven(PortType type);

/// The port's complex excitation, `amplitude * exp(j * phase_deg * pi/180)`,
/// under the `e^{+j*omega*t}` convention -- so a positive phase leads.
///
/// This is the quantity the assembly uses directly, with no sign applied:
/// for a current port it is the right-hand side of the port row, and for a
/// voltage port it is the prescribed potential on the live side. The
/// orientation note above is why there is nothing to multiply it by.
///
/// A negative `amplitude` is legal and means what it says: -1 at 0 degrees
/// is the same excitation as +1 at 180.
std::complex<double> port_amplitude(const Port& port);

/// The deterministic sign convention used when a cut port gives no
/// `current_direction`: returns `normal` normalised and oriented so that its
/// largest-magnitude component is positive (ties broken in x, y, z order).
///
/// Which of the two directions this picks is arbitrary with respect to what
/// the user meant -- but it must be *reproducible*, or a regression test's
/// reported sign could flip between builds on one unchanged mesh. That is
/// the whole requirement on it.
///
/// Throws std::invalid_argument if `normal` is zero.
Vec3 canonical_orientation(const Vec3& normal);

/// Resolves a cut's direction **d** from its face normal and the user's
/// hint: `d = sign(hint . n) * n`, with both normalised.
///
/// The hint is only ever used for its sign, so it need not be the normal --
/// any vector within 90 degrees of the intended direction gives the same
/// **d**. `min_alignment` guards the case where it is nearly *in* the cut
/// plane and so picks no side at all.
///
/// Throws std::invalid_argument if either vector is zero, or if
/// `|hint_hat . n_hat| < min_alignment` -- the hint does not say which way
/// through the cut the current goes.
Vec3 resolve_cut_direction(const Vec3& hint, const Vec3& normal, double min_alignment = 0.1);

/// Where Phi lives, which is *not* a free choice: at omega = 0 the Phi
/// equation `-div(sigma*grad(Phi)) = 0` is identically empty wherever
/// sigma = 0, so Phi DOFs in an insulator would be zero rows and the matrix
/// singular. True therefore at DC, and in the reduced (eps -> 0) variant;
/// false only for the full-wave system, where displacement current gives
/// Phi something to do in the dielectric.
bool phi_on_conductors_only(const Problem& problem);

/// Number of linear systems this problem implies: one per frequency, or one
/// at DC.
std::size_t num_solves(const Problem& problem);

}  // namespace aphi_solver
