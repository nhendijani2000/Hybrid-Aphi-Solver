#pragma once

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
/// Both halves matter, and the formatting is the larger one: `ostream <<` on a
/// double runs through locale-aware machinery and is typically an order of
/// magnitude slower than `std::to_chars`, which also produces the **shortest
/// representation that round-trips exactly** -- so the file is smaller *and*
/// loses nothing. The block writes then remove one `fwrite` per value.
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

/// Writes the nodal potential to `path` -- `potential.out` by convention.
///
/// One line per P2 node: index, position, `Re(Φ)`, `Im(Φ)`, and a one-character
/// status so a reader can tell a solved value from a prescribed one, and can see
/// which nodes are on a cut and therefore discontinuous.
WriteStats write_potential(const std::string& path, const NodalPotential& potential,
                           const Solution& solution);

}  // namespace aphi_solver
