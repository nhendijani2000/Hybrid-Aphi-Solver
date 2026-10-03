#pragma once

/// A gauge-free electric scalar potential, recovered from the fields.
///
/// WHY THIS EXISTS
///
/// `08_MixedPort_Interior` measures a single-potential terminal on an interior
/// conductor/air face losing 3.1 % of its inductance and 77 % of the nearby
/// current density to the choice of spanning tree, against 1e-12 and 1e-11 for
/// the same wire with both caps on the outer boundary. `GAUGE_CHOICE.md` Sec.
/// 15.9 then derives that **no tree can fix it**: every route to pinning `psi`
/// on an interior face either leaves it path dependent, hence unprotected, or
/// imposes a spurious flux constraint, hence wrong.
///
/// Phi read out of the gauged A-Phi system is therefore a tree artefact
/// wherever the terminals are not on the boundary. This recovers a Phi that is
/// not, by solving its OWN boundary value problem driven by the field:
///
///     < beta grad(Phi'), grad(lam) >  =  - < beta E, grad(lam) >
///
/// with `beta = sigma + j w eps`. That is Stysch's Sec. 3.1 / 6.3 construction
/// (`GAUGE_CHOICE.md` Sec. 11.8), and the third independent route in the
/// literature to a unique Phi after Chew's `chi` and Ostrowski & Hiptmair's
/// electro-quasistatic gauge.
///
/// WHAT IT FIXES, AND WHAT IT DOES NOT. Measured in `tools/gauge_spike`,
/// milestone 4, on a stub whose exact DC resistance is known:
///
/// | quantity                | from the gauged system | from this  |
/// |-------------------------|------------------------|------------|
/// | R                       | 11.5007 uOhm           | **11.494253, exact to 8 digits** |
/// | tree dependence of R    | 3.8e-02                | **1.1e-14** |
/// | reactance               | L ~ 3.3 nH             | **absent**  |
///
/// It fixes **Phi**, and through Phi the resistance and the resistive part of a
/// terminal voltage. It does **not** touch `E`, `B`, `H` or `J`, which are
/// inputs here and come out of the gauged solve unchanged -- at a mixed port
/// they move 12 %, 91 %, 91 % and 77 % with the tree and remain wrong.
///
/// The reason it works at all is the reason it cannot do more: the three
/// solutions differ by a **solenoidal** field, `div(beta dE) = 0`, so `div`
/// annihilates the difference. This recovers the irrotational (galvanic) part
/// of the solution exactly and is blind to the solenoidal part -- which is
/// where the flux, the inductance and the local field errors live. An
/// energy-based inductance does not rescue them either: measured at 42 %
/// spread across three trees (`tools/gauge_spike/RESULTS_SPIKEC.md`).
///
/// So: ship this for `R`, and do not read it as "mixed ports are fixed".

#include <complex>
#include <vector>

#include "aphi_solver/dof_map.hpp"
#include "aphi_solver/mesh.hpp"
#include "aphi_solver/postprocess.hpp"
#include "aphi_solver/problem_binding.hpp"

namespace aphi_solver {

/// The recovered potential, and enough about the solve to judge it.
struct GaugeFreePotential {
    /// One value per mesh VERTEX. P1, not P2: the source `div(beta E)` is
    /// piecewise constant on a tet, so the extra P2 freedom buys nothing here
    /// while doubling the system.
    std::vector<std::complex<double>> phi_vertex;

    /// One per port, in `BoundProblem::ports` order: the terminal's single
    /// recovered potential. A voltage port's comes back as prescribed, which
    /// is a free check that the assembly is right.
    std::vector<std::complex<double>> port_potential;

    int unknowns = 0;
    double assemble_ms = 0.0;
    double solve_ms = 0.0;
    bool converged = false;

    /// Relative residual of the recovered system, so a silent failure is not
    /// mistaken for a gauge-free answer.
    double residual = 0.0;
};

/// Recover `Phi` from `fields.e_tet`. `omega` is needed for `beta`; at
/// `omega == 0` the weight degenerates to `sigma` alone, which is the DC
/// conduction problem and is correct.
///
/// Terminals are handled exactly as the main solve handles them: every vertex
/// of a current-driven port shares ONE unknown (`PhiDof::Port`), and a voltage
/// port's vertices are prescribed. The reference comes from the same place it
/// does in the main solve, so the recovered `Phi` is directly comparable with
/// the one written next to it.
GaugeFreePotential recover_gauge_free_potential(const Mesh& mesh, const BoundProblem& bound,
                                                const DofMap& dofs, const FieldOutput& fields,
                                                double omega);

/// Write the recovered potential beside the solved one, in the same plain
/// column format every other .out file uses: position, then the value.
///
/// One row per mesh VERTEX, not per P2 node, because that is what was
/// recovered. The header states what the file is and what it does not fix, so
/// a reader who opens only this file is not left to infer it.
WriteStats write_gauge_free_potential(const std::string& path, const Mesh& mesh,
                                      const BoundProblem& bound, const GaugeFreePotential& gfp,
                                      const RunInfo* run = nullptr);

}  // namespace aphi_solver
