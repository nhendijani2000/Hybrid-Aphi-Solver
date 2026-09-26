#include "aphi_solver/assembly.hpp"

#include <array>
#include <stdexcept>
#include <string>

namespace aphi_solver {

using Complex = std::complex<double>;

FormulationScales formulation_scales(Formulation3 f, double omega) {
    const Complex jw(0.0, omega);
    FormulationScales s;
    switch (f) {
        case Formulation3::Natural:
            s.row = Complex(1.0, 0.0);
            s.column = Complex(1.0, 0.0);
            s.symmetric = false;
            return s;
        case Formulation3::RowScaled:
            if (omega == 0.0) {
                throw std::invalid_argument(
                    "formulation_scales: the row-scaled formulation divides the Phi rows by "
                    "j*omega, which is undefined at DC. At omega = 0 the system decouples and "
                    "there is no coupling block to symmetrise -- use the natural form.");
            }
            s.row = Complex(1.0, 0.0) / jw;
            s.column = Complex(1.0, 0.0);
            s.symmetric = true;
            return s;
        case Formulation3::ScaledPhi:
            if (omega == 0.0) {
                throw std::invalid_argument(
                    "formulation_scales: the scaled-potential formulation substitutes "
                    "Phi = j*omega*Phi', so a prescribed potential would need dividing by "
                    "j*omega, which is undefined at DC. Use the natural form.");
            }
            s.row = Complex(1.0, 0.0);
            s.column = jw;
            s.symmetric = true;
            return s;
    }
    throw std::invalid_argument("formulation_scales: unknown formulation");
}

namespace {

/// One `matrix[i][j] += value`, refusing to lose the term if the symbolic
/// pass did not provide for it.
inline void add_at(const SparsityPattern& pattern, std::vector<Complex>& values, int row, int col,
                   const Complex& value) {
    const int slot = pattern.find_slot(row, col);
    if (slot < 0) {
        throw std::logic_error("assemble: no slot for entry (" + std::to_string(row) + ", " +
                               std::to_string(col) +
                               "). The sparsity pattern and the scatter disagree about the "
                               "matrix's shape; dropping the term would give a quietly wrong "
                               "matrix.");
    }
    values[static_cast<std::size_t>(slot)] += value;
}

/// The prescribed value of local Phi node `q` expressed in the SCALED
/// unknown, `Phi_given / c`. Eliminating a column multiplies this by that
/// column's matrix coefficient, which already carries `c`.
inline Complex prescribed(const TetDofs& d, int q, const Complex& column_scale) {
    const Complex v = d.phi_fixed[static_cast<std::size_t>(q)];
    return v == Complex(0.0, 0.0) ? v : v / column_scale;
}

}  // namespace

AssembledSystem assemble(const BoundProblem& bound, const Mesh& mesh, const DofMap& dofs,
                         const SparsityPattern& pattern, double omega, Formulation3 formulation) {
    const FormulationScales scales = formulation_scales(formulation, omega);
    const Complex jw(0.0, omega);

    AssembledSystem out;
    out.matrix = SparseMatrixZ::from_pattern(pattern.rows, pattern.cols, pattern.row_ptr,
                                             pattern.col_index);
    out.rhs.assign(static_cast<std::size_t>(dofs.num_total), Complex(0.0, 0.0));
    std::vector<Complex>& values = out.matrix.mutable_values();

    // Reused for every tet. ~2.6 KB on the stack, never heap.
    std::array<Complex, 36> aa{};
    std::array<Complex, 60> ap{};
    std::array<Complex, 100> pp{};

    // Bodies then tets, so the coefficients are a per-body cost rather than
    // a per-tet one. The bodies partition the tets, so each is visited once.
    for (const BoundBody& body : bound.bodies) {
        const ElementCoefficients c = element_coefficients(body, omega);

        for (int tet : body.tets) {
            const TetGeometry g = compute_tet_geometry(mesh, tet);
            const TetDofs d = dofs.local_dofs(tet, mesh, bound);

            kernel_AA(g, c, aa.data());
            kernel_APhi(g, c, ap.data());
            kernel_PhiPhi(g, c, pp.data());

            // --- (A, A) -------------------------------------------------
            for (int p = 0; p < 6; ++p) {
                const DofEntry& rp = d.edge[static_cast<std::size_t>(p)];
                if (rp.index < 0) continue;
                for (int q = 0; q < 6; ++q) {
                    const DofEntry& cq = d.edge[static_cast<std::size_t>(q)];
                    // A prescribed edge has a = 0, so it contributes
                    // nothing to the right-hand side either -- unlike a
                    // prescribed Phi, which generally does.
                    if (cq.index < 0) continue;
                    add_at(pattern, values, rp.index, cq.index,
                           (rp.coeff * cq.coeff) * aa[static_cast<std::size_t>(p * 6 + q)]);
                }
            }

            // --- (A, Phi) and (Phi, A) ----------------------------------
            //
            // One kernel serves both: the assembled (Phi,A) block is
            // r * alpha * C^T = r * j*omega * (beta*C)^T, and `ap` already
            // holds beta*C.
            for (int p = 0; p < 6; ++p) {
                const DofEntry& e = d.edge[static_cast<std::size_t>(p)];
                if (e.index < 0) continue;
                for (int q = 0; q < 10; ++q) {
                    const DofEntry& n = d.phi[static_cast<std::size_t>(q)];
                    const Complex block = ap[static_cast<std::size_t>(p * 10 + q)];

                    // (A, Phi): the Phi DOF is the column, so a prescribed
                    // one moves to the A row's right-hand side.
                    if (n.index >= 0) {
                        add_at(pattern, values, e.index, n.index,
                               scales.column * (e.coeff * n.coeff) * block);
                    } else {
                        // The eliminated column carries its own matrix
                        // coefficient -- including the column scale --
                        // times the prescribed value OF THE SCALED
                        // UNKNOWN, which is Phi_given / c. The two factors
                        // of c cancel, so a prescribed potential reaches
                        // the right-hand side unscaled; multiplying by c
                        // and using the physical value would be wrong by c
                        // under ScaledPhi.
                        const Complex fixed = prescribed(d, q, scales.column);
                        if (fixed != Complex(0.0, 0.0)) {
                            out.rhs[static_cast<std::size_t>(e.index)] -=
                                scales.column * e.coeff * block * fixed;
                        }
                    }

                    // (Phi, A): the Phi DOF is now the row. A prescribed
                    // row has no equation at all, so it is simply absent.
                    if (n.index >= 0) {
                        add_at(pattern, values, n.index, e.index,
                               scales.row * jw * (n.coeff * e.coeff) * block);
                    }
                }
            }

            // --- (Phi, Phi) ---------------------------------------------
            for (int p = 0; p < 10; ++p) {
                const DofEntry& rp = d.phi[static_cast<std::size_t>(p)];
                if (rp.index < 0) continue;
                for (int q = 0; q < 10; ++q) {
                    const DofEntry& cq = d.phi[static_cast<std::size_t>(q)];
                    const Complex block = pp[static_cast<std::size_t>(p * 10 + q)];
                    if (cq.index >= 0) {
                        add_at(pattern, values, rp.index, cq.index,
                               scales.row * scales.column * (rp.coeff * cq.coeff) * block);
                    } else {
                        const Complex fixed = prescribed(d, q, scales.column);
                        if (fixed != Complex(0.0, 0.0)) {
                            out.rhs[static_cast<std::size_t>(rp.index)] -=
                                scales.row * scales.column * rp.coeff * block * fixed;
                        }
                    }
                }
            }
        }
    }

    // --- port excitations ------------------------------------------------
    //
    // The two types are duals: a current port knows its row's equation and
    // solves for V, a voltage port knows V and loses that equation.
    for (std::size_t k = 0; k < bound.ports.size(); ++k) {
        const int row = dofs.port_index[k];
        if (is_current_driven(bound.ports[k].type)) {
            // The port row is the terminal's current balance, scattered
            // from the tets and therefore carrying the row scale, so the
            // driving current must carry it too.
            out.rhs[static_cast<std::size_t>(row)] += scales.row * bound.ports[k].amplitude;
        } else {
            // Nothing was scattered into this row: the terminal is a
            // prescribed Phi, so the port has no column. The row holds the
            // constraint itself, `V' = V_given / c`, which keeps the matrix
            // square and the row and column both otherwise empty -- so the
            // symmetry of the symmetric formulations survives.
            //
            // The diagonal is 1 rather than something comparable to the
            // Phi block's ~1e8, which leaves one badly scaled row. That is
            // equilibration's job (`equilibration.hpp`) and the right value
            // is a measurement once the solver exists, not a guess now.
            //
            // The current through a voltage port is therefore NOT a
            // residual of this row any more; extraction has to re-form the
            // terminal's current balance. See docs/ASSEMBLY_PLAN.md Sec. 10.
            add_at(pattern, values, row, row, Complex(1.0, 0.0));
            out.rhs[static_cast<std::size_t>(row)] += dofs.port_value[k] / scales.column;
        }
    }

    return out;
}

}  // namespace aphi_solver
