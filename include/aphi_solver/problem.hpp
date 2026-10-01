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

/// How the assembled system is scaled, from `docs/CONDITIONING.md`. This is
/// **not** `Formulation` above, which decides where Phi lives; this decides
/// only how the Phi rows and columns are scaled, and therefore whether the
/// matrix comes out symmetric. The two were both called "formulation" in
/// earlier drafts, which is exactly the confusion this comment exists to
/// prevent.
///
/// All three describe the same problem and have the same solution. They
/// differ in a row scale `r` on every Phi row and a column scale `c` from
/// the substitution `Phi = c * Phi'`; the matrix is symmetric exactly when
/// `c == r * j*omega`.
enum class Conditioning {
    Natural,    ///< r = 1, c = 1. Unsymmetric; needs a general (LU) solver.
    RowScaled,  ///< r = 1/(j*omega), c = 1. Symmetric. CONDITIONING.md F1.
    ScaledPhi   ///< r = 1, c = j*omega. Symmetric, and the unknown is Phi'.
};

/// The input-file keyword for one conditioning choice, and back. `name_to`
/// returns false for anything else rather than guessing.
const char* conditioning_keyword(Conditioning c);
bool conditioning_from_keyword(const std::string& word, Conditioning& out);

/// True for the two that need `j*omega` and so cannot be used at DC.
bool conditioning_needs_ac(Conditioning c);

/// Which factorization does the work.
///
/// `Internal` is this project's own LDL^T and is the DEFAULT, deliberately:
/// a clean checkout builds and runs the whole suite with nothing installed.
/// See SOLVER_PLAN.md Sec. 12 -- having a library available as an option is a
/// different thing from depending on one.
///
/// `Mumps` is only selectable in a build configured with `-DAPHI_WITH_MUMPS=ON`.
/// Asking for it in a build without it is a clear error at load time, not a
/// silent fall back to the internal solver: a run that quietly used a different
/// factorization than the one asked for would invalidate any timing comparison
/// made with it, which is the main reason to have the option at all.
enum class SolverBackend {
    Internal,  ///< this project's LDL^T. Always available.
    Mumps      ///< MUMPS complex symmetric indefinite. Needs APHI_WITH_MUMPS.
};

const char* solver_backend_keyword(SolverBackend b);
bool solver_backend_from_keyword(const std::string& word, SolverBackend& out);

/// True when this build can actually run that backend.
bool solver_backend_available(SolverBackend b);

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

// ---------------------------------------------------------------------------
// Post-processing requests.
//
// One [postprocess] section is one picture (or one line plot). The solver does
// not render anything: it validates the request and records it, and the
// rendering tool reads it back. That split keeps ParaView out of this build and
// still catches a typo in a second rather than after a solve.

enum class PostGeometry {
    Body,   ///< the cells of one named Physical Volume
    Plane,  ///< a cut plane, normal to x, y or z
    Points  ///< an explicit list of probe coordinates
};

enum class PostPlane { XY, YZ, ZX };

enum class PostField { Phi, A, E, B, H, J };

/// What to draw. The names are deliberately unlike each other, because the
/// first three are easy to confuse and the difference is not cosmetic:
/// a vector phasor's tip traces an ELLIPSE, so it has two characteristic
/// lengths, and three different numbers are all reasonably called "the
/// magnitude". See docs/ComplexVectorPhasorConcept.md for the derivation.
///
/// Writing Ehat = P + jQ and E(t) = P cos(wt) - Q sin(wt), with semi-axes
/// a (major) and b (minor):
///
///     ComplexMagnitude  N = sqrt(|P|^2+|Q|^2) = sqrt(a^2+b^2) = sqrt(2)*RMS
///     MagnitudeAtPhase  |P cos(th) - Q sin(th)|          -- depends on th
///     Peak              a = max over th                  -- what it reaches
///
/// with b <= MagnitudeAtPhase <= Peak <= ComplexMagnitude <= sqrt(2)*Peak.
enum class PostDisplay {
    ComplexMagnitude,  ///< DEFAULT. Phase-independent; what every figure in the
                       ///< regression suite is already coloured by.
    MagnitudeAtPhase,  ///< the instantaneous magnitude at `phase_deg`
    Peak,              ///< the semi-major axis: the largest value ever reached
    AxialRatio,        ///< b/a in [0,1]; 0 linear, 1 circular. A diagnostic that
                       ///< says whether ComplexMagnitude may be read as a peak.
    Phase,             ///< arg of a scalar, or of one component of a vector
    Vector,            ///< the instantaneous vector at `phase_deg`
    Real,              ///< P alone
    Imag               ///< Q alone
};

enum class PostComponent { None, X, Y, Z };

/// The colour map. The spelling here is this project's; the renderer maps each
/// to a ParaView preset name, which is the thing that varies between versions.
///
/// RAINBOW IS THE DEFAULT, and that is a reversal: the renderer used Viridis
/// first. Two reasons to prefer the rainbow. It is what every existing figure
/// in the suite already uses -- cases 01-03 with `Jet`, case 05 with
/// `Rainbow Uniform` -- so a new picture sits beside the published ones without
/// reading as a different tool. And it is what the field is read with: blue
/// cold, red hot, no legend needed to know which end is which.
///
/// `Rainbow` rather than `Jet` because Jet is not perceptually uniform: it
/// compresses the greens and stretches the cyans, which invents contour-like
/// bands that are not in the data. Rainbow Uniform looks the same at a glance
/// and does not do that. `Jet` remains available for matching an older figure.
enum class PostColormap {
    Rainbow,     ///< blue -> red, perceptually uniform. The default.
    Jet,         ///< blue -> red, the classic. Bands where the data does not.
    Turbo,       ///< a modern Jet, uniform
    CoolToWarm,  ///< blue -> white -> red; for signed data about zero
    Viridis,     ///< dark blue -> yellow; uniform, colour-blind safe
    BlueToRed,   ///< a wider rainbow
    BlackBody,   ///< black -> red -> yellow -> white
    Grayscale,
    XRay         ///< white -> black
};

/// One `[postprocess]` section.
///
/// WHICH FIELDS MATTER DEPENDS ON THE OTHERS, which is why they are optional
/// here and why the parser rejects any that the choices make meaningless rather
/// than ignoring them. A `phase_deg` that is silently dropped because the
/// display is phase-independent is exactly the kind of thing that makes someone
/// distrust the solver.
///
///     geometry = body    -> `body` is required
///     geometry = plane   -> `plane` is required, `offset` optional (0)
///     geometry = points  -> `points` is required
///     display needs a phase (MagnitudeAtPhase, Vector) -> `phase_deg`, else 0
///     display = Phase on a VECTOR field -> `component` is required
///
/// `body` is not checked against the mesh here. Like `Body::volume`, that is
/// the binding stage's job -- see the contract at the top of input_file.hpp.
struct PostprocessRequest {
    std::string name;  ///< the section's own handle, e.g. "PP1"
    int line = 0;      ///< its line in the input file; 0 if built in code

    PostGeometry geometry = PostGeometry::Body;
    PostField field = PostField::E;
    PostDisplay display = PostDisplay::ComplexMagnitude;

    std::string body;                  ///< geometry = Body
    PostPlane plane = PostPlane::XY;   ///< geometry = Plane
    double offset = 0.0;               ///< geometry = Plane, in `length_unit`
    std::vector<Vec3> points;          ///< geometry = Points, in `length_unit`

    /// The instant wt = phase_deg, under the e^{+jwt} convention: a field
    /// stored as P + jQ is drawn as P cos(phase) - Q sin(phase). Zero unless
    /// the display uses it, and the parser refuses it when it does not.
    double phase_deg = 0.0;

    PostComponent component = PostComponent::None;  ///< display = Phase, vector

    /// Applies to every display that produces a colour bar, which is all of
    /// them except `vector` -- and even there the arrows are coloured by it.
    PostColormap colormap = PostColormap::Rainbow;
};

/// True for the fields that are vectors. `Phi` is the only scalar, and the
/// distinction decides which displays are legal.
inline bool post_field_is_vector(PostField f) { return f != PostField::Phi; }

/// True for the displays that read `phase_deg`.
inline bool post_display_uses_phase(PostDisplay d) {
    return d == PostDisplay::MagnitudeAtPhase || d == PostDisplay::Vector;
}

struct Problem {
    std::string mesh_file;  ///< as written, relative to the input file

    /// Where solve_mesh writes potential.out, potential.vtk and anything else
    /// it produces. Resolved relative to the INPUT FILE, not the working
    /// directory, so a case runs identically from anywhere -- which is what
    /// lets a regression case own its own outputs. Empty means the working
    /// directory, which is the old behaviour.
    std::string output_dir;
    LengthUnit length_unit = LengthUnit::Metre;

    AnalysisType type = AnalysisType::DC;

    /// One entry per solve; empty at DC. A `sweep` in the input file is
    /// expanded into this list by the parser, so nothing downstream sees a
    /// second way of saying the same thing.
    std::vector<double> frequencies;

    Formulation formulation = Formulation::FullWave;
    OuterBoundary outer = OuterBoundary::FluxTangential;

    /// Defaults to Natural: it is the only one defined at every frequency,
    /// and choosing between the three is a measurement nobody can make until
    /// a solver reports a condition number (ASSEMBLY_PLAN Sec. 9 item 5).
    Conditioning conditioning = Conditioning::Natural;

    /// `[solver] backend`. Internal unless the file says otherwise.
    SolverBackend backend = SolverBackend::Internal;

    std::vector<Body> bodies;
    std::vector<Port> ports;

    /// Empty unless the file asks for pictures. Order is the order written.
    std::vector<PostprocessRequest> postprocess;
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
