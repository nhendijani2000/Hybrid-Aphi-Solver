#pragma once

#include <array>
#include <complex>
#include <cstdio>
#include <string>
#include <vector>

#include "aphi_solver/dof_map.hpp"
#include "aphi_solver/mesh.hpp"
#include "aphi_solver/problem.hpp"
#include "aphi_solver/problem_binding.hpp"

namespace aphi_solver {

/// Post-processing: turning a solution vector back into fields.
/// `docs/POSTPROCESSING_PLAN.md`.

/// What one solve produced, with everything needed to interpret it.
///
/// The conditioning is carried here for one reason: under
/// `Conditioning::ScaledPhi` the unknown is `Φ' = Φ/(jω)`, not `Φ`. Bundling it
/// with the vector means the `jω` is applied **once**, in `phi_scale()`, rather
/// than at every call site where one of them will eventually forget. That
/// omission is listed as a negative control in `docs/SOLVER_PLAN.md` §8 and is
/// the likeliest way to produce a wrong field that still looks smooth.
struct Solution {
    std::vector<std::complex<double>> x;  ///< the raw vector, laid out [a | phi | V]
    Conditioning conditioning = Conditioning::Natural;
    double omega = 0.0;

    /// The factor turning a stored Φ unknown into the physical Φ: `jω` under
    /// `ScaledPhi`, 1 otherwise.
    std::complex<double> phi_scale() const;
};

/// Why a P2 node has the potential it has -- because not every node has one, and
/// one kind has two.
enum class PhiStatus : char {
    Free,      ///< solved for
    Terminal,  ///< on a port terminal: the port's single unknown
    Cut,       ///< on an internal cut: Φ is DISCONTINUOUS here, see below
    Fixed,     ///< prescribed
    Absent     ///< Φ does not live here at all (an insulator at DC)
};

/// Φ at every P2 node: the mesh's vertices first, then one per edge midpoint.
///
/// **Mid-edge values are exact, not interpolated.** Φ lives in a P2 space, so an
/// edge's midpoint carries a genuine unknown that the solve determined. Averaging
/// the two endpoints -- which is the right thing for `B`, `H` and `E`, and what a
/// P1 code has to do -- would replace the quadratic with its linear interpolant
/// and throw away the term that makes the space second order.
/// `docs/POSTPROCESSING_PLAN.md` §2.
///
/// **A node on an internal cut has two values, not one.** The plus side reads the
/// port's unknown and the grounded minus side reads zero; that is the whole
/// mechanism by which a single-valued Φ can drive a closed loop. `value` holds
/// the **plus side** and `status` is `Cut`, so the discontinuity is visible
/// rather than silently resolved.
struct NodalPotential {
    std::vector<Vec3> position;
    std::vector<std::complex<double>> value;
    std::vector<PhiStatus> status;

    int num_vertices = 0;  ///< entries [0, num_vertices) are mesh vertices
    int num_edges = 0;     ///< the rest are edge midpoints, in mesh edge order

    int size() const { return static_cast<int>(value.size()); }
};

NodalPotential potential_at_nodes(const Mesh& mesh, const BoundProblem& bound, const DofMap& dofs,
                                  const Solution& solution);

/// An append-only text buffer that formats with `std::to_chars` and writes in
/// large blocks.
///
/// **The formatting is about 90 % of the win, not the buffering** -- measured,
/// because the reverse is the natural assumption. On 200000 values:
///
///     to_chars + 1 MB buffer   7.5 ms     <- this
///     to_chars, fwrite/line   11.6 ms     buffering is worth 1.6x
///     ostream, one big write 108.1 ms     to_chars is worth 14.4x
///     ostream, value by value 115.9 ms    15.5x overall
///     ostream + std::endl     464.2 ms    61.9x
///
/// `ostream <<` on a double runs through locale-aware machinery; `std::to_chars`
/// does not, and its shortest representation **round-trips exactly**, so the file
/// is smaller *and* lossless.
///
/// Buffering earns only 1.6x for a reason worth knowing: a `FILE*` is **already**
/// buffered by the C runtime, so `fwrite` per line was never a system call per
/// line. This buffer is a second layer that removes `fwrite`'s per-call overhead,
/// not syscalls. What genuinely costs is `std::endl`, which forces a real flush
/// every line and is 62x slower -- the one thing this code never does.
class TextBuffer {
public:
    explicit TextBuffer(std::FILE* out, std::size_t capacity = 1u << 20);
    ~TextBuffer();

    TextBuffer(const TextBuffer&) = delete;
    TextBuffer& operator=(const TextBuffer&) = delete;

    void put(char c);
    void put(const char* text);
    void put(const std::string& text) { put(text.c_str()); }
    void put(int value);
    void put(long long value);

    /// Shortest round-tripping representation. No precision argument, because a
    /// fixed precision would be both longer and lossy.
    void put(double value);

    void flush();

    std::size_t bytes_written() const { return written_; }

private:
    void make_room(std::size_t n);

    std::FILE* out_;
    std::vector<char> buffer_;
    std::size_t used_ = 0;
    std::size_t written_ = 0;
};

struct WriteStats {
    std::size_t bytes = 0;
    int nodes = 0;
    double milliseconds = 0.0;
};

/// What produced a result file.
///
/// Written into the header of every `.out` so a file says how it was made:
/// which mesh and input, how many unknowns, how well the solve converged and
/// what each stage cost. Without it a result file is undated and unattributed,
/// and two runs on different meshes are indistinguishable once the console has
/// scrolled away -- which is exactly the position a mesh comparison leaves you
/// in.
///
/// Optional everywhere. Passing `nullptr` writes the old header unchanged, so
/// the writers stay usable from tests that have no solve behind them.
struct RunInfo {
    std::string input_file;   ///< the .aphi that was run
    std::string mesh_file;    ///< the .msh it named
    std::string ordering;     ///< natural | rcm | amd
    int num_tets = 0;
    int unknowns = 0;         ///< the size of the linear system
    int free_edges = 0;       ///< A degrees of freedom actually solved
    int tree_edges = 0;       ///< zeroed by the tree-cotree gauge
    int dirichlet_edges = 0;  ///< pinned by n x A = 0
    std::size_t stored_nonzeros = 0;
    std::size_t factor_nnz = 0;  ///< nnz(L); the fill-in, and the memory driver
    double backward_error = 0.0;
    double residual = 0.0;
    double assemble_ms = 0.0;
    double analyze_ms = 0.0;
    double factorize_ms = 0.0;
    double solve_ms = 0.0;
};

/// Writes the nodal potential to `path` -- `potential.out` by convention.
///
/// One line per P2 node: index, position, `Re(Φ)`, `Im(Φ)`, and a one-character
/// status so a reader can tell a solved value from a prescribed one, and can see
/// which nodes are on a cut and therefore discontinuous.

WriteStats write_potential(const std::string& path, const NodalPotential& potential,
                           const Solution& solution, const RunInfo* run = nullptr);


/// VTK's node order for a quadratic tetrahedron (`VTK_QUADRATIC_TETRA`, type
/// 24), expressed as our local edge index for each of its six mid-edge slots.
///
/// VTK wants the mid-edge nodes in the order (0,1), (1,2), (0,2), (0,3), (1,3),
/// (2,3). Ours are numbered (0,1), (0,2), (0,3), (1,2), (2,3), (1,3)
/// (`kTetLocalEdgeVerts`). So the two disagree in four of six places, and a
/// wrong mapping produces a mesh that still renders -- smoothly, plausibly, and
/// with the quadratic nodes on the wrong edges. One test checks each slot's
/// POSITION against the midpoint VTK expects there, which is the only way to
/// catch that.
inline constexpr std::array<int, 6> kVtkQuadraticTetEdgeOrder = {0, 3, 1, 2, 5, 4};

/// A complex vector field value. Plain `std::array` rather than a `Vec3` of
/// complexes, because nothing here does vector algebra on it -- it is carried,
/// averaged componentwise, and written.
using Vec3C = std::array<std::complex<double>, 3>;

/// The fields of `docs/FORMULATION.md` Sec. 1, reconstructed from a solution.
///
/// **Each quantity is stored where it actually lives, not forced to a common
/// grid.** `docs/POSTPROCESSING_PLAN.md` Sec. 2:
///
/// - `Φ` is P2, so `phi_node` holds **exact solved values**.
///   Nothing is interpolated or averaged, and averaging the mid-edge values
///   would throw away the term that makes the space second order.
/// - `B = curl A` is genuinely **constant per tet**, because `curl W_e` is.
///   `b_tet` is that constant; `b_node` is a volume-weighted average and is
///   *smoothing*, not refinement.
/// - `A` and `E` are linear within a tet. They are evaluated **at each vertex
///   from each incident tet** and then averaged -- not evaluated once at the
///   centroid and smeared, which would lose the linear variation.
///
/// Both the per-tet and the per-vertex forms are kept so the smoothing can be
/// seen rather than assumed: where they disagree, the mesh is too coarse.
///
/// `Φ` here is gauge-dependent and `A` is too; `B`, `H` and `E` are not.
/// See Sec. 9 of the plan.
struct FieldOutput {
    /// Exact P2 values of `Φ`, one per P2 node: `DofMap`'s numbering, so
    /// entries `[0, num_vertices)` are the mesh vertices and the rest are edge
    /// midpoints in mesh edge order.
    std::vector<std::complex<double>> phi_node;

    // As computed: one per tet. `b_tet` is exact; `a_tet` and `e_tet` are the
    // value at the tet's centroid, which for a linear field is its mean.
    std::vector<Vec3C> a_tet, b_tet, e_tet;

    /// The conductivity of each tet's body, so `J = sigma E` can be formed
    /// per tet -- which is where it is well defined. At a node on a
    /// conductor/insulator interface it is not.
    std::vector<double> sigma_tet;

    /// `A`, `B`, `H` and `E` at **every P2 node**, indexed exactly like
    /// `phi_node`, so a field and the potential can be read at the same node
    /// index without a second convention.
    ///
    /// The two halves are produced differently, and the difference is the whole
    /// point of `docs/POSTPROCESSING_PLAN.md` Sec. 2:
    ///
    /// - **At a vertex**: evaluated there from each incident tet, then averaged
    ///   with the tet volumes as weights.
    /// - **At an edge midpoint**: the mean of the two endpoint values. These
    ///   quantities have **no mid-edge degree of freedom** -- unlike `Φ`, whose
    ///   mid-edge value is a genuine unknown the solve determined and which is
    ///   therefore never averaged. Averaging is right here and wrong there, and
    ///   that asymmetry is deliberate.
    std::vector<Vec3C> a_node, b_node, h_node, e_node;

    /// `J = sigma E` at every P2 node.
    ///
    /// **`sigma` is a property of the tet, not of the node**, so this is only
    /// defined where a node's incident tets all share one body -- exactly where
    /// `on_material_interface` is 0. At an interface node `J` is genuinely
    /// two-valued (a finite current density on the conductor side, zero on the
    /// insulator side), and no single number is right, so this is left at zero
    /// there and the flag says why.
    ///
    /// `j_tet` has no such trouble and is the one to integrate.
    std::vector<Vec3C> j_node;

    /// `sigma E` per tet, which is where `J` is unambiguous. Zero throughout an
    /// insulator, which is correct and not a gap.
    std::vector<Vec3C> j_tet;

    /// Total volume of the tets that contributed to each **P2 node** -- size
    /// `num_p2_nodes`, mid-edge nodes included, since those are accumulated
    /// directly rather than averaged from their endpoints. Zero means no tet
    /// reached that node at all, which leaves its fields at zero.
    ///
    /// This is `A`'s weight. `E` and `B` have their own, because at a material
    /// interface each averages over only one side and the three differ.
    std::vector<double> vertex_weight;

    /// 1 at a P2 node whose incident tets do **not** all belong to one body.
    ///
    /// **The nodal fields are meaningless at these nodes**, and quietly so.
    /// `E`'s normal component genuinely jumps across a conductor/insulator
    /// interface -- `J_n = 0` at a free conductor surface requires it -- and so
    /// does `H`'s tangential component across a change of `mu`. Averaging over
    /// tets on both sides produces a value that is neither.
    ///
    /// Measured on the 50 Hz cylinder: inside the wire the nodal `E_z` is
    /// -999.999 V/m against an exact -1000, with a transverse component of
    /// 9e-07 V/m; **on the surface** the same average gives a transverse
    /// component of 2705 V/m, larger than the axial field itself.
    ///
    /// Threshold this to 0 before reading a nodal field, or use the per-tet
    /// arrays, which are exact and have no interface to straddle.
    std::vector<unsigned char> on_material_interface;

    int num_interface_nodes = 0;

    int num_vertices = 0;  ///< entries [0, num_vertices) of the node arrays
    int num_edges = 0;     ///< the rest, in mesh edge order
    int num_orphan_vertices = 0;

    int num_nodes() const { return static_cast<int>(phi_node.size()); }
};

/// Reconstructs `A`, `B`, `H`, `E` and `Φ` over the whole mesh.
///
/// `H = B / mu` uses the tet's own body, so it is discontinuous across a
/// material interface exactly where `B` is continuous. The per-vertex `H` at
/// such an interface averages across it and is therefore meaningless there;
/// `b_node` is the one to look at on a boundary.
FieldOutput compute_fields(const Mesh& mesh, const BoundProblem& bound, const DofMap& dofs,
                           const Solution& solution);

/// Which quantity a writer should emit.
///
/// One file per field is the default because it is what a reader actually
/// wants: opening `B_field.vtk` in ParaView gives a source with `B` on it and
/// nothing else to pick through, and `E_field.out` is a table of `E` rather
/// than thirty-one columns of which six are `E`. `All` writes the combined
/// file, which is still the one to open when comparing two fields in a single
/// ParaView session.
enum class FieldSet { All, Potential, A, B, H, E, J };

/// The short name used in filenames and headers: "potential", "A", "B", ...
const char* field_set_name(FieldSet which);

/// Writes `Φ`, `A`, `B`, `H` and `E` at every P2 node as text -- the plan's
/// `WriteSolution` (`docs/POSTPROCESSING_PLAN.md` §3), in this project's
/// naming. `which` selects one field, or `All` for every one of them.
///
/// One line per node: index, position, then each field as a real/imaginary
/// pair per component, and finally the material-interface flag. Thirty-one
/// columns, so the header names them.
///
/// This is the companion to `write_vtk`, not a replacement: the VTK carries the
/// same nodal values **plus** the per-cell ones, which are the exact,
/// unaveraged quantities. Use this when the values are wanted as numbers rather
/// than as a picture.

/// Terminal current, from the volume integral of `J` along an axis.
///
///     I = -(1/length) * integral_conductor  J . axis  dV
///
/// This is the standard volume form of the terminal current, and it is used in
/// preference to a surface integral over the terminal face for a concrete
/// reason: `J = sigma E` is exact per tetrahedron and has no material interface
/// to straddle, whereas a face integral would need `J` at nodes lying exactly
/// on the conductor boundary, where its normal component jumps. See Sec. 12 of
/// docs/FIELD_POSTPROCESSING.md.
///
/// The conductor is every tet with `sigma > 0`; current cannot flow elsewhere,
/// so no body tag is needed.
///
/// The SIGN is current flowing OUT of the high-potential terminal, which is the
/// convention that makes `Z = V/I` positive for a passive load. `axis` points
/// from the low-potential terminal to the high one and need not be normalised
/// beyond being a unit vector; `length` is the conductor's extent along it.
///
/// `axis` and `length` are geometry and cannot be inferred from the fields, so
/// they are arguments rather than guesses.
std::complex<double> terminal_current(const Mesh& mesh, const FieldOutput& fields,
                                      const Vec3& axis, double length);

/// `Z = V / I`, for a case driven by a known terminal voltage.
///
/// `Re(Z)` is the resistance and `Im(Z)/omega` the inductance. Trivial, but
/// written once here so every caller divides the same way round -- inverting it
/// is a mistake that produces a plausible-looking number.
inline std::complex<double> impedance(double voltage, std::complex<double> current) {
    return std::complex<double>(voltage, 0.0) / current;
}

WriteStats write_solution(const std::string& path, const NodalPotential& potential,
                          const FieldOutput& fields, const Solution& solution,
                          FieldSet which = FieldSet::All, const RunInfo* run = nullptr);

/// Writes a legacy VTK unstructured grid for ParaView.
///
/// The mesh is written as **quadratic** tetrahedra, so the P2 mid-edge values go
/// into the file as the genuine unknowns they are rather than being discarded or
/// averaged away. ParaView interpolates them correctly.
///
/// Point data: `phi_real`, `phi_imag`, `phi_magnitude`, and `phi_present`.
/// That last one is 0 where Φ does not live -- which under
/// `formulation = reduced`, or at DC, is the whole insulating region. Writing
/// those nodes as 0 without saying so would make the air look like it was
/// solved and found to be at zero volts. Threshold on `phi_present` in ParaView
/// to see only where the potential means something.
///
/// **`phi_imag` is gauge-dependent -- do not plot it as a result.** The discrete
/// system is exactly gauge-invariant (`alpha = j*omega*beta`), its freedom is
/// `psi` in P1, and the tree-cotree constraint picks one representative. A
/// different tree shifts `Phi` by `-j*omega*psi`, which for real `psi` is almost
/// entirely imaginary. Measured on the 50 Hz cylinder: `Re(Phi)` matches the
/// exact `z/l` to 4.0e-07, while `Im(Phi)` violates the problem's own mirror
/// antisymmetry by twice its own peak. `E = -j*omega*A - grad(Phi)` is the
/// gauge-invariant object. `docs/POSTPROCESSING_PLAN.md` section 9.
///
/// Cell data: `body_tag`, the mesh's physical-volume tag, so conductor and
/// insulator can be separated without consulting the input file.
WriteStats write_vtk(const std::string& path, const Mesh& mesh, const NodalPotential& potential,
                     const Solution& solution, const FieldOutput* fields = nullptr,
                     FieldSet which = FieldSet::All, const RunInfo* run = nullptr);

/// Writes the `[postprocess]` requests as JSON, for `tools/postprocess.py`.
///
/// THE SOLVER RENDERS NOTHING. It owns the grammar, the validation and the
/// naming; the Python tool owns ParaView. This file is the whole interface
/// between them, which is why it carries more than the requests themselves:
///
///   * **the resolved `body_tag`** for each `geometry = body` request. The
///     input file names a Physical Volume; the VTK carries integer tags. Only
///     binding knows the mapping, so resolving it here means the renderer never
///     has to open a mesh.
///   * **metres.** `offset` and `points` are written in the file's
///     `length_unit`, but `scale_mesh_to_metres` has already run by the time
///     any VTK is written, so the manifest converts them. A renderer slicing at
///     a millimetre offset in a metre mesh would miss by a factor of 1000 and
///     produce an empty, plausible-looking picture.
///   * **which file holds which field, and what its arrays are called.**
///     Hard-coding `E_cell_real` in Python would be a second place to edit
///     whenever a writer changes its mind.
///
/// `frequency` is the solve this manifest belongs to, and `suffix` is the
/// filename suffix that distinguishes a sweep's solves (empty for a single
/// one), so the renderer reads the files that go with it.
///
/// Writes nothing and reports zero bytes when there are no requests: a file
/// that asked for no pictures should not grow an empty manifest.
WriteStats write_postprocess_manifest(const std::string& path, const Problem& problem,
                                      const BoundProblem& bound, double frequency,
                                      const std::string& suffix, const RunInfo* run = nullptr);

}  // namespace aphi_solver
