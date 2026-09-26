#include "aphi_solver/assembly.hpp"

#include <algorithm>
#include <array>
#include <stdexcept>
#include <string>

namespace aphi_solver {

using Complex = std::complex<double>;

FormulationScales formulation_scales(Conditioning f, double omega) {
    const Complex jw(0.0, omega);
    FormulationScales s;
    switch (f) {
        case Conditioning::Natural:
            s.row = Complex(1.0, 0.0);
            s.column = Complex(1.0, 0.0);
            s.symmetric = false;
            return s;
        case Conditioning::RowScaled:
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
        case Conditioning::ScaledPhi:
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

/// One `matrix[i][j] += value` by global index, for the few entries that do
/// not come from a tet. Refuses to lose the term if the symbolic pass did
/// not provide for it.
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

/// The live global DOFs of one tet, sorted and uniqued, and where each local
/// DOF sits among them.
///
/// Locating entries -- not computing them -- was 68 % of assembly time when
/// each was found by its own `find_slot`: 75 ms of 107 on the cylinder, over
/// 3.6 M searches (docs/ASSEMBLY_PLAN.md Sec. 11). A tet contributes at most
/// 16 distinct global DOFs, so all ~225 of its entries lie in 16 rows and 16
/// columns.
struct TetLive {
    int n = 0;
    std::array<int, 16> live{};
    std::array<int, 6> edge_at{};
    std::array<int, 10> phi_at{};
};

void gather_live(const TetDofs& d, TetLive& s) {
    s.n = 0;
    for (const DofEntry& e : d.edge) {
        if (e.index >= 0) s.live[static_cast<std::size_t>(s.n++)] = e.index;
    }
    for (const DofEntry& e : d.phi) {
        if (e.index >= 0) s.live[static_cast<std::size_t>(s.n++)] = e.index;
    }
    const auto first = s.live.begin();
    std::sort(first, first + s.n);
    // Several local DOFs can share one global index -- every node of a port
    // terminal reads that port's single unknown -- so the list is uniqued,
    // exactly as `build_sparsity` does before forming its pairs. Duplicates
    // then land on the same slot and accumulate, which is what is wanted.
    s.n = static_cast<int>(std::unique(first, first + s.n) - first);

    const auto position = [&](int global) {
        return global < 0 ? -1
                          : static_cast<int>(std::lower_bound(first, first + s.n, global) - first);
    };
    for (int i = 0; i < 6; ++i) {
        s.edge_at[static_cast<std::size_t>(i)] = position(d.edge[static_cast<std::size_t>(i)].index);
    }
    for (int i = 0; i < 10; ++i) {
        s.phi_at[static_cast<std::size_t>(i)] = position(d.phi[static_cast<std::size_t>(i)].index);
    }
}

/// Fills `out` with the value-array index of every (row, column) pair among
/// `s.live`, row-major with stride `s.n`.
///
/// Each row is swept in ASCENDING column order so every `lower_bound`
/// resumes where the previous one stopped instead of at the row's beginning.
/// That is the whole optimisation, and it is why the live list is sorted.
///
/// Every pair must exist, because `build_sparsity` formed the pattern from
/// exactly these lists; a miss means the symbolic and numeric passes
/// disagree, so it throws rather than dropping a term.
void locate(const SparsityPattern& pattern, const TetLive& s, int* out,
            std::size_t capacity) {
    // `out` holds an n x n block. Checking that here, once per tet, turns a
    // mismatch between this and whoever sized the buffer into an exception
    // instead of a write past the end -- which in a `ScatterMap` means heap
    // corruption during construction, before any test can look at it.
    if (static_cast<std::size_t>(s.n) * static_cast<std::size_t>(s.n) > capacity) {
        throw std::logic_error("locate: a tet with " + std::to_string(s.n) +
                               " live DOFs needs " + std::to_string(s.n * s.n) +
                               " slots but was given " + std::to_string(capacity) + ".");
    }
    for (int r = 0; r < s.n; ++r) {
        const int row = s.live[static_cast<std::size_t>(r)];
        auto it = pattern.col_index.begin() + pattern.row_ptr[static_cast<std::size_t>(row)];
        const auto end =
            pattern.col_index.begin() + pattern.row_ptr[static_cast<std::size_t>(row) + 1];
        for (int col = 0; col < s.n; ++col) {
            const int want = s.live[static_cast<std::size_t>(col)];
            it = std::lower_bound(it, end, want);
            if (it == end || *it != want) {
                throw std::logic_error(
                    "assemble: no slot for entry (" + std::to_string(row) + ", " +
                    std::to_string(want) +
                    "). The sparsity pattern and the scatter disagree about the matrix's "
                    "shape; dropping the term would give a quietly wrong matrix.");
            }
            out[r * s.n + col] = static_cast<int>(it - pattern.col_index.begin());
        }
    }
}

/// One tet's scatter positions, however they were obtained -- worked out on
/// the spot, or read from a `ScatterMap`. The scatter loop is written once,
/// against this.
struct TetScatter {
    int n = 0;
    const int* edge_at = nullptr;
    const int* phi_at = nullptr;
    const int* slot = nullptr;  ///< n x n, row-major

    int of(int row_at, int col_at) const { return slot[row_at * n + col_at]; }
};

/// The prescribed value of local Phi node `q` expressed in the SCALED
/// unknown, `Phi_given / c`. Eliminating a column multiplies this by that
/// column's matrix coefficient, which already carries `c`.
inline Complex prescribed(const TetDofs& d, int q, const Complex& column_scale) {
    const Complex v = d.phi_fixed[static_cast<std::size_t>(q)];
    return v == Complex(0.0, 0.0) ? v : v / column_scale;
}

}  // namespace

// ---------------------------------------------------------------- ScatterMap

ScatterMap ScatterMap::build(const SparsityPattern& pattern, const DofMap& dofs, const Mesh& mesh,
                             const BoundProblem& bound) {
    const int num_tets = mesh.num_tets();
    ScatterMap m;
    m.live_.assign(static_cast<std::size_t>(num_tets), 0);
    m.pos_.assign(static_cast<std::size_t>(num_tets) * 16, -1);
    m.offset_.assign(static_cast<std::size_t>(num_tets) + 1, 0);

    // Two passes, so `slot_` is allocated once at its exact size rather than
    // grown. The blocks are not all the same size: a tet on the Dirichlet
    // boundary, or outside Phi's support, has fewer than 16 live DOFs.
    std::vector<TetLive> live(static_cast<std::size_t>(num_tets));
    for (int t = 0; t < num_tets; ++t) {
        const std::size_t u = static_cast<std::size_t>(t);
        gather_live(dofs.local_dofs(t, mesh, bound), live[u]);
        m.live_[u] = live[u].n;
        m.offset_[u + 1] =
            m.offset_[u] + static_cast<std::size_t>(live[u].n) * static_cast<std::size_t>(live[u].n);
    }
    m.slot_.assign(m.offset_[static_cast<std::size_t>(num_tets)], -1);

    for (int t = 0; t < num_tets; ++t) {
        const std::size_t u = static_cast<std::size_t>(t);
        for (int i = 0; i < 6; ++i) {
            m.pos_[u * 16 + static_cast<std::size_t>(i)] =
                live[u].edge_at[static_cast<std::size_t>(i)];
        }
        for (int i = 0; i < 10; ++i) {
            m.pos_[u * 16 + 6 + static_cast<std::size_t>(i)] =
                live[u].phi_at[static_cast<std::size_t>(i)];
        }
        locate(pattern, live[u], m.slot_.data() + m.offset_[u],
               m.offset_[u + 1] - m.offset_[u]);
    }
    return m;
}

std::size_t ScatterMap::bytes() const {
    return live_.size() * sizeof(int) + pos_.size() * sizeof(int) +
           offset_.size() * sizeof(std::size_t) + slot_.size() * sizeof(int);
}

// ------------------------------------------------------------------ assembly

AssembledSystem make_system(const SparsityPattern& pattern, int num_unknowns) {
    AssembledSystem out;
    out.matrix = SparseMatrixZ::from_pattern(pattern.rows, pattern.cols, pattern.row_ptr,
                                             pattern.col_index);
    out.rhs.assign(static_cast<std::size_t>(num_unknowns), Complex(0.0, 0.0));
    return out;
}

void refill(AssembledSystem& system, const BoundProblem& bound, const Mesh& mesh,
            const DofMap& dofs, const SparsityPattern& pattern, double omega,
            Conditioning formulation, const ScatterMap* scatter) {
    const FormulationScales scales = formulation_scales(formulation, omega);
    const Complex jw(0.0, omega);

    if (system.matrix.rows() != pattern.rows || system.matrix.nnz() != pattern.nnz()) {
        throw std::invalid_argument(
            "refill: this system was not made from this pattern. Their shapes must agree, or "
            "the slots the scatter writes to mean something else.");
    }
    if (system.rhs.size() != static_cast<std::size_t>(dofs.num_total)) {
        throw std::invalid_argument("refill: the right-hand side is not one entry per unknown.");
    }
    if (scatter != nullptr && scatter->num_tets() != mesh.num_tets()) {
        throw std::invalid_argument("refill: the ScatterMap was built for a different mesh (" +
                                    std::to_string(scatter->num_tets()) + " tets against " +
                                    std::to_string(mesh.num_tets()) + ").");
    }

    std::vector<Complex>& values = system.matrix.mutable_values();

    // The scatter accumulates, so anything left from a previous frequency has
    // to go first. Zeroing an array already held is 0.6 ms against 11.4 to
    // allocate a fresh one, which is why `make_system` is separate.
    std::fill(values.begin(), values.end(), Complex(0.0, 0.0));
    std::fill(system.rhs.begin(), system.rhs.end(), Complex(0.0, 0.0));

    // Reused for every tet. ~4.8 KB on the stack, never heap.
    std::array<Complex, 36> aa{};
    std::array<Complex, 60> ap{};
    std::array<Complex, 100> pp{};
    TetLive live;
    std::array<int, 256> scratch{};

    // Bodies then tets, so the coefficients are a per-body cost rather than
    // a per-tet one. The bodies partition the tets, so each is visited once.
    for (const BoundBody& body : bound.bodies) {
        const ElementCoefficients c = element_coefficients(body, omega);

        for (int tet : body.tets) {
            const TetGeometry g = compute_tet_geometry(mesh, tet);
            const TetDofs d = dofs.local_dofs(tet, mesh, bound);

            TetScatter where;
            if (scatter != nullptr) {
                where = {scatter->live_count(tet), scatter->edge_at(tet), scatter->phi_at(tet),
                         scatter->slots(tet)};
            } else {
                gather_live(d, live);
                locate(pattern, live, scratch.data(), scratch.size());
                where = {live.n, live.edge_at.data(), live.phi_at.data(), scratch.data()};
            }

            kernel_AA(g, c, aa.data());
            kernel_APhi(g, c, ap.data());
            kernel_PhiPhi(g, c, pp.data());

            // --- (A, A) -------------------------------------------------
            for (int p = 0; p < 6; ++p) {
                const DofEntry& rp = d.edge[static_cast<std::size_t>(p)];
                if (rp.index < 0) continue;
                const int row_at = where.edge_at[p];
                for (int q = 0; q < 6; ++q) {
                    const DofEntry& cq = d.edge[static_cast<std::size_t>(q)];
                    // A prescribed edge has a = 0, so it contributes
                    // nothing to the right-hand side either -- unlike a
                    // prescribed Phi, which generally does.
                    if (cq.index < 0) continue;
                    values[static_cast<std::size_t>(where.of(row_at, where.edge_at[q]))] +=
                        (rp.coeff * cq.coeff) * aa[static_cast<std::size_t>(p * 6 + q)];
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
                const int row_at = where.edge_at[p];
                for (int q = 0; q < 10; ++q) {
                    const DofEntry& n = d.phi[static_cast<std::size_t>(q)];
                    const Complex block = ap[static_cast<std::size_t>(p * 10 + q)];

                    // (A, Phi): the Phi DOF is the column, so a prescribed
                    // one moves to the A row's right-hand side.
                    if (n.index >= 0) {
                        values[static_cast<std::size_t>(where.of(row_at, where.phi_at[q]))] +=
                            scales.column * (e.coeff * n.coeff) * block;
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
                            system.rhs[static_cast<std::size_t>(e.index)] -=
                                scales.column * e.coeff * block * fixed;
                        }
                    }

                    // (Phi, A): the Phi DOF is now the row. A prescribed
                    // row has no equation at all, so it is simply absent.
                    if (n.index >= 0) {
                        values[static_cast<std::size_t>(where.of(where.phi_at[q], row_at))] +=
                            scales.row * jw * (n.coeff * e.coeff) * block;
                    }
                }
            }

            // --- (Phi, Phi) ---------------------------------------------
            for (int p = 0; p < 10; ++p) {
                const DofEntry& rp = d.phi[static_cast<std::size_t>(p)];
                if (rp.index < 0) continue;
                const int row_at = where.phi_at[p];
                for (int q = 0; q < 10; ++q) {
                    const DofEntry& cq = d.phi[static_cast<std::size_t>(q)];
                    const Complex block = pp[static_cast<std::size_t>(p * 10 + q)];
                    if (cq.index >= 0) {
                        values[static_cast<std::size_t>(where.of(row_at, where.phi_at[q]))] +=
                            scales.row * scales.column * (rp.coeff * cq.coeff) * block;
                    } else {
                        const Complex fixed = prescribed(d, q, scales.column);
                        if (fixed != Complex(0.0, 0.0)) {
                            system.rhs[static_cast<std::size_t>(rp.index)] -=
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
            system.rhs[static_cast<std::size_t>(row)] += scales.row * bound.ports[k].amplitude;
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
            system.rhs[static_cast<std::size_t>(row)] += dofs.port_value[k] / scales.column;
        }
    }
}

AssembledSystem assemble(const BoundProblem& bound, const Mesh& mesh, const DofMap& dofs,
                         const SparsityPattern& pattern, double omega, Conditioning formulation) {
    AssembledSystem out = make_system(pattern, dofs.num_total);
    refill(out, bound, mesh, dofs, pattern, omega, formulation);
    return out;
}

}  // namespace aphi_solver
