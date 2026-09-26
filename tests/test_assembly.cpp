// Tests for from_pattern and global assembly.
//
// Steps 5-7 of `docs/ASSEMBLY_PLAN.md` §8. This is the first point at which
// a whole matrix exists, so the checks are about the matrix as a whole --
// symmetry where it is claimed, the exact relation between the three
// formulations, and agreement with the independent triplet path.

#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "aphi_solver/assembly.hpp"
#include "aphi_solver/gmsh_reader.hpp"

using namespace aphi_solver;
using Complex = std::complex<double>;

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

bool near(Complex a, Complex b, double tol) { return std::abs(a - b) <= tol; }

Mesh make_cube() {
    Mesh m;
    m.nodes = {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
               {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};
    m.tets = {{0, 1, 2, 6}, {0, 2, 3, 6}, {0, 3, 7, 6},
              {0, 7, 4, 6}, {0, 4, 5, 6}, {0, 5, 1, 6}};
    m.build_topology();
    m.tet_tags = {1, 1, 1, 1, 1, 1};
    m.physical_names = {{3, 1, "metal"}, {2, 10, "bottom"}, {2, 11, "top"}};
    return m;
}

void tag_face(Mesh& m, int a, int b, int c, int tag) {
    TaggedFace f;
    f.nodes = {a, b, c};
    std::sort(f.nodes.begin(), f.nodes.end());
    f.tag = tag;
    m.tagged_boundary_faces.push_back(f);
}

/// A cube problem. `frequency` selects DC or an AC analysis.
BoundProblem bind_cube(Mesh& m, double frequency, double current = 1.0) {
    tag_face(m, 0, 1, 2, 10);
    tag_face(m, 0, 2, 3, 10);
    tag_face(m, 4, 5, 6, 11);
    tag_face(m, 4, 6, 7, 11);

    Problem p;
    if (frequency == 0.0) {
        p.type = AnalysisType::DC;
    } else {
        p.type = AnalysisType::Frequency;
        p.frequencies = {frequency};
        p.formulation = Formulation::FullWave;
    }
    Body b;
    b.name = "B1";
    b.volume = "metal";
    b.sigma = 5.8e7;
    b.line = 1;
    p.bodies = {b};
    Port p1;
    p1.name = "P1";
    p1.type = PortType::BoundaryCurrent;
    p1.surface = {"bottom"};
    p1.amplitude = current;
    p1.line = 2;
    Port p2;
    p2.name = "P2";
    p2.type = PortType::BoundaryVoltage;
    p2.surface = {"top"};
    // NON-ZERO, and off the real axis. An earlier draft of these tests
    // prescribed 0 V, which made every wrong treatment of a prescribed
    // potential look right: it hid that the value was never imposed at all
    // (leaving constant Phi in the null space) and that ScaledPhi has to
    // divide it by c.
    p2.amplitude = 2.5;
    p2.phase_deg = 30.0;
    p2.line = 3;
    p.ports = {p1, p2};
    return bind_to_mesh(p, m);
}

Complex entry(const SparseMatrixZ& m, int row, int col) {
    const auto& rp = m.row_ptr();
    const auto& ci = m.col_index();
    for (int k = rp[static_cast<std::size_t>(row)]; k < rp[static_cast<std::size_t>(row) + 1]; ++k) {
        if (ci[static_cast<std::size_t>(k)] == col) return m.values()[static_cast<std::size_t>(k)];
    }
    return Complex(0.0, 0.0);
}

double max_abs(const SparseMatrixZ& m) {
    double worst = 0.0;
    for (const Complex& v : m.values()) worst = std::max(worst, std::abs(v));
    return worst;
}

/// Values equal to within rounding. Used where the two paths could legally
/// differ in the last bits.
bool same_values(const AssembledSystem& a, const AssembledSystem& b) {
    if (a.matrix.nnz() != b.matrix.nnz() || a.rhs.size() != b.rhs.size()) return false;
    const double scale = std::max(max_abs(a.matrix), 1e-300);
    for (std::size_t k = 0; k < a.matrix.values().size(); ++k) {
        if (std::abs(a.matrix.values()[k] - b.matrix.values()[k]) > 1e-12 * scale) return false;
    }
    double rhs_scale = 1e-300;
    for (const Complex& v : a.rhs) rhs_scale = std::max(rhs_scale, std::abs(v));
    for (std::size_t k = 0; k < a.rhs.size(); ++k) {
        if (std::abs(a.rhs[k] - b.rhs[k]) > 1e-12 * rhs_scale) return false;
    }
    return true;
}

/// Exactly equal, bit for bit. The ScatterMap changes only WHERE the slots
/// come from, not the values or the order they are accumulated in, so it has
/// no licence to differ at all -- and demanding that catches a subtle wrong
/// slot that a tolerance would hide.
bool bitwise_equal(const AssembledSystem& a, const AssembledSystem& b) {
    return a.matrix.values() == b.matrix.values() && a.rhs == b.rhs;
}

// ---------------------------------------------------------------------------

void test_from_pattern() {
    // A 3x3 with a deliberate hole at (1,0).
    const std::vector<int> row_ptr = {0, 2, 3, 5};
    const std::vector<int> col_index = {0, 2, 1, 0, 2};
    const SparseMatrixZ m = SparseMatrixZ::from_pattern(3, 3, row_ptr, col_index);

    check(m.compressed(), "from_pattern returns an already-compressed matrix");
    check(m.rows() == 3 && m.cols() == 3, "shape is what was asked for");
    check(m.nnz() == 5, "nnz matches the pattern");

    bool all_zero = true;
    for (const Complex& v : m.values()) {
        if (v != Complex(0.0, 0.0)) all_zero = false;
    }
    check(all_zero, "every value starts at zero");
    check(m.row_ptr() == row_ptr && m.col_index() == col_index,
          "the structure is adopted unchanged");

    // Malformed structures are refused rather than producing a matrix whose
    // operations would quietly misbehave.
    const auto rejects = [](const std::vector<int>& rp, const std::vector<int>& ci, int rows,
                            int cols) {
        try {
            SparseMatrixZ::from_pattern(rows, cols, rp, ci);
        } catch (const std::invalid_argument&) {
            return true;
        } catch (...) {
            return false;
        }
        return false;
    };
    check(rejects({0, 2}, {0, 1}, 3, 3), "row_ptr of the wrong length is refused");
    check(rejects({0, 2, 3, 5}, {0, 1, 1, 2}, 3, 3), "row_ptr not ending at nnz is refused");
    // The other side of that check: {0, 2, 3, 4} over four columns is a
    // perfectly good matrix (rows of 2, 1 and 1). An earlier draft of this
    // test expected it to be refused, so the boundary is worth pinning.
    check(!rejects({0, 2, 3, 4}, {0, 1, 1, 2}, 3, 3),
          "a valid row_ptr ending exactly at nnz is accepted");
    check(!rejects({0, 0, 0, 0}, {}, 3, 3), "a structurally empty matrix is accepted");
    check(rejects({0, 2, 2, 4}, {0, 5, 1, 2}, 3, 3), "an out-of-range column is refused");
    check(rejects({0, 2, 3, 5}, {2, 0, 1, 0, 2}, 3, 3),
          "a row whose columns descend is refused -- find_slot binary-searches them");
    check(rejects({0, 2, 3, 5}, {0, 0, 1, 0, 2}, 3, 3), "a duplicate column is refused");
}

// The assembled matrix must agree with the triplet path the project already
// had. They share no code: one writes into slots a symbolic pass provided,
// the other sorts a pile of triplets.
void test_agrees_with_the_triplet_path() {
    Mesh m = make_cube();
    const BoundProblem b = bind_cube(m, 0.0);
    const DofMap d = build_dof_map(b, m);
    const SparsityPattern sp = build_sparsity(d, b, m);
    const AssembledSystem sys = assemble(b, m, d, sp, 0.0, Conditioning::Natural);

    // The same assembly, by hand, into a triplet matrix.
    SparseMatrixZ ref(d.num_total, d.num_total);
    std::array<Complex, 36> aa{};
    std::array<Complex, 60> ap{};
    std::array<Complex, 100> pp{};
    for (const BoundBody& body : b.bodies) {
        const ElementCoefficients c = element_coefficients(body, 0.0);
        for (int t : body.tets) {
            const TetGeometry g = compute_tet_geometry(m, t);
            const TetDofs td = d.local_dofs(t, m, b);
            kernel_AA(g, c, aa.data());
            kernel_APhi(g, c, ap.data());
            kernel_PhiPhi(g, c, pp.data());
            for (int p = 0; p < 6; ++p) {
                if (td.edge[p].index < 0) continue;
                for (int q = 0; q < 6; ++q) {
                    if (td.edge[q].index < 0) continue;
                    ref.add(td.edge[p].index, td.edge[q].index,
                            td.edge[p].coeff * td.edge[q].coeff * aa[p * 6 + q]);
                }
            }
            for (int p = 0; p < 6; ++p) {
                if (td.edge[p].index < 0) continue;
                for (int q = 0; q < 10; ++q) {
                    if (td.phi[q].index < 0) continue;
                    ref.add(td.edge[p].index, td.phi[q].index,
                            td.edge[p].coeff * td.phi[q].coeff * ap[p * 10 + q]);
                    ref.add(td.phi[q].index, td.edge[p].index, Complex(0.0, 0.0));
                }
            }
            for (int p = 0; p < 10; ++p) {
                if (td.phi[p].index < 0) continue;
                for (int q = 0; q < 10; ++q) {
                    if (td.phi[q].index < 0) continue;
                    ref.add(td.phi[p].index, td.phi[q].index,
                            td.phi[p].coeff * td.phi[q].coeff * pp[p * 10 + q]);
                }
            }
        }
    }
    // The same constraint row assembly writes for a voltage port.
    for (std::size_t k = 0; k < d.port_is_fixed.size(); ++k) {
        if (d.port_is_fixed[k]) ref.add(d.port_index[k], d.port_index[k], Complex(1.0, 0.0));
    }
    ref.compress();

    // At DC omega is zero, so the (Phi,A) block is zero and the reference
    // above (which adds explicit zeros there) matches entry for entry.
    double worst = 0.0;
    for (int r = 0; r < d.num_total; ++r) {
        for (int k = sys.matrix.row_ptr()[static_cast<std::size_t>(r)];
             k < sys.matrix.row_ptr()[static_cast<std::size_t>(r) + 1]; ++k) {
            const int col = sys.matrix.col_index()[static_cast<std::size_t>(k)];
            worst = std::max(worst, std::abs(sys.matrix.values()[static_cast<std::size_t>(k)] -
                                             entry(ref, r, col)));
        }
    }
    check(worst < 1e-9 * max_abs(sys.matrix),
          "the pattern-based assembly matches the triplet path entry for entry (worst " +
              std::to_string(worst) + ")");
}

void test_dc_structure() {
    Mesh m = make_cube();
    const BoundProblem b = bind_cube(m, 0.0);
    const DofMap d = build_dof_map(b, m);
    const SparsityPattern sp = build_sparsity(d, b, m);
    const AssembledSystem sys = assemble(b, m, d, sp, 0.0, Conditioning::Natural);

    check(sys.matrix.rows() == d.num_total, "the matrix is num_total square");
    check(sys.rhs.size() == static_cast<std::size_t>(d.num_total), "the RHS has one entry per DOF");

    // At DC the system decouples: alpha = 0 kills both the mass term and
    // the whole (Phi,A) block, so the matrix is block diagonal and real.
    bool real_only = true;
    for (const Complex& v : sys.matrix.values()) {
        if (v.imag() != 0.0) real_only = false;
    }
    check(real_only, "at DC every entry is real -- alpha is exactly zero");

    double worst_coupling = 0.0;
    for (int r = d.num_a; r < d.num_a + d.num_phi; ++r) {
        for (int k = sys.matrix.row_ptr()[static_cast<std::size_t>(r)];
             k < sys.matrix.row_ptr()[static_cast<std::size_t>(r) + 1]; ++k) {
            if (sys.matrix.col_index()[static_cast<std::size_t>(k)] < d.num_a) {
                worst_coupling = std::max(
                    worst_coupling, std::abs(sys.matrix.values()[static_cast<std::size_t>(k)]));
            }
        }
    }
    check(worst_coupling == 0.0,
          "the (Phi,A) block is exactly zero at DC -- the decoupling of FORMULATION.md Sec. 2");

    // The current port's amplitude reaches its own row and nothing else.
    // Asserting rhs == 1 there would be wrong: that row is a Phi equation,
    // so it also receives the eliminated column of the prescribed terminal
    // next door. Differencing against the same problem at 0 A isolates what
    // the current itself contributed.
    Mesh m0 = make_cube();
    const BoundProblem b0 = bind_cube(m0, 0.0, 0.0);
    const DofMap d0 = build_dof_map(b0, m0);
    const SparsityPattern sp0 = build_sparsity(d0, b0, m0);
    const AssembledSystem quiet = assemble(b0, m0, d0, sp0, 0.0, Conditioning::Natural);

    double worst_elsewhere = 0.0;
    for (int i = 0; i < d.num_total; ++i) {
        const Complex delta =
            sys.rhs[static_cast<std::size_t>(i)] - quiet.rhs[static_cast<std::size_t>(i)];
        if (i == d.port_index[0]) {
            check(near(delta, Complex(1.0, 0.0), 1e-12),
                  "the 1 A current port contributes exactly 1 to its own row of the RHS");
        } else {
            worst_elsewhere = std::max(worst_elsewhere, std::abs(delta));
        }
    }
    check(worst_elsewhere == 0.0, "and touches no other row");

    // The voltage port's row is the constraint `V = V_given` and nothing
    // else: its terminal is a prescribed Phi, so the tets scatter nothing
    // into that row and the port has no column anywhere.
    const int vrow = d.port_index[1];
    // 2.5 * exp(j * 30 deg), written out rather than read back from the
    // DofMap, so the check does not agree with the code by construction.
    const Complex v_given(2.5 * 0.86602540378443864676, 1.25);
    check(near(sys.rhs[static_cast<std::size_t>(vrow)], v_given, 1e-14),
          "the prescribed 2.5 V at 30 deg reaches its own row of the RHS");
    check(sys.matrix.row_ptr()[static_cast<std::size_t>(vrow) + 1] -
                  sys.matrix.row_ptr()[static_cast<std::size_t>(vrow)] ==
              1,
          "that row holds exactly one entry");
    check(near(entry(sys.matrix, vrow, vrow), Complex(1.0, 0.0), 1e-15),
          "and it is the unit diagonal of the constraint");

    // The point of all that: a prescribed potential is what makes the DC
    // system solvable. Without it a constant Phi costs nothing, and the
    // matrix is singular -- which is what the code did before the voltage
    // port's value was imposed at all.
    std::vector<Complex> constant_phi(static_cast<std::size_t>(d.num_total), Complex(0.0, 0.0));
    for (int i = d.num_a; i < d.num_total; ++i) {
        constant_phi[static_cast<std::size_t>(i)] = Complex(1.0, 0.0);
    }
    const std::vector<Complex> image = sys.matrix.matvec(constant_phi);
    double largest = 0.0;
    for (int i = d.num_a; i < d.num_total; ++i) {
        largest = std::max(largest, std::abs(image[static_cast<std::size_t>(i)]));
    }
    check(largest > 1e-6 * max_abs(sys.matrix),
          "a constant Phi is NOT in the null space -- the prescribed terminal is the reference");
}

// The three formulations describe the same problem with the Phi rows scaled
// by r and the Phi columns by c. Undoing both must reproduce the natural
// system exactly -- which is checkable without a solver.
void test_formulations_are_row_and_column_scalings() {
    Mesh m = make_cube();
    const double f = 1e6;
    const double omega = 2.0 * M_PI * f;
    const BoundProblem b = bind_cube(m, f);
    const DofMap d = build_dof_map(b, m);
    const SparsityPattern sp = build_sparsity(d, b, m);

    const AssembledSystem natural = assemble(b, m, d, sp, omega, Conditioning::Natural);
    const auto is_phi_row = [&](int i) { return i >= d.num_a; };

    for (const Conditioning which : {Conditioning::RowScaled, Conditioning::ScaledPhi}) {
        const FormulationScales s = formulation_scales(which, omega);
        const AssembledSystem scaled = assemble(b, m, d, sp, omega, which);

        double worst = 0.0;
        for (int r = 0; r < d.num_total; ++r) {
            for (int k = scaled.matrix.row_ptr()[static_cast<std::size_t>(r)];
                 k < scaled.matrix.row_ptr()[static_cast<std::size_t>(r) + 1]; ++k) {
                const int col = scaled.matrix.col_index()[static_cast<std::size_t>(k)];
                Complex undone = scaled.matrix.values()[static_cast<std::size_t>(k)];
                if (is_phi_row(r)) undone /= s.row;
                if (is_phi_row(col)) undone /= s.column;
                worst = std::max(worst,
                                 std::abs(undone - natural.matrix.values()[static_cast<std::size_t>(k)]));
            }
        }
        check(worst < 1e-6 * max_abs(natural.matrix),
              "undoing the row and column scales reproduces the natural system (worst " +
                  std::to_string(worst) + ")");

        // The right-hand side has to follow the same rule, and it is where
        // the column scale is easiest to get wrong. Every row's RHS is its
        // own row scale times the natural one:
        //
        //   an A row              1        (a prescribed Phi's column
        //                                  contributes c * coeff * Phi/c,
        //                                  so the c's cancel)
        //   a Phi row             r
        //   a current port row    r        (its row IS a Phi equation)
        //   a voltage port row    1/c      (not a physical equation at all
        //                                  -- it is V' = V_given / c)
        double rhs_worst = 0.0;
        double rhs_scale = 0.0;
        for (const Complex& v : natural.rhs) rhs_scale = std::max(rhs_scale, std::abs(v));
        for (int i = 0; i < d.num_total; ++i) {
            Complex row_scale = Complex(1.0, 0.0);
            if (i >= d.num_a) row_scale = s.row;
            for (std::size_t k = 0; k < d.port_is_fixed.size(); ++k) {
                if (d.port_is_fixed[k] && d.port_index[k] == i) {
                    row_scale = Complex(1.0, 0.0) / s.column;
                }
            }
            const Complex undone = scaled.rhs[static_cast<std::size_t>(i)] / row_scale;
            rhs_worst = std::max(rhs_worst,
                                 std::abs(undone - natural.rhs[static_cast<std::size_t>(i)]));
        }
        check(rhs_worst < 1e-6 * rhs_scale,
              "and the same for the right-hand side (worst " + std::to_string(rhs_worst) + " of " +
                  std::to_string(rhs_scale) + ")");
    }
}

void test_symmetry() {
    Mesh m = make_cube();
    const double omega = 2.0 * M_PI * 1e6;
    const BoundProblem b = bind_cube(m, 1e6);
    const DofMap d = build_dof_map(b, m);
    const SparsityPattern sp = build_sparsity(d, b, m);

    const auto worst_asymmetry = [&](const AssembledSystem& s) {
        double worst = 0.0;
        for (int r = 0; r < d.num_total; ++r) {
            for (int k = s.matrix.row_ptr()[static_cast<std::size_t>(r)];
                 k < s.matrix.row_ptr()[static_cast<std::size_t>(r) + 1]; ++k) {
                const int col = s.matrix.col_index()[static_cast<std::size_t>(k)];
                worst = std::max(worst, std::abs(s.matrix.values()[static_cast<std::size_t>(k)] -
                                                 entry(s.matrix, col, r)));
            }
        }
        return worst;
    };

    const AssembledSystem natural = assemble(b, m, d, sp, omega, Conditioning::Natural);
    const AssembledSystem rows = assemble(b, m, d, sp, omega, Conditioning::RowScaled);
    const AssembledSystem cols = assemble(b, m, d, sp, omega, Conditioning::ScaledPhi);

    const double scale = max_abs(natural.matrix);
    check(worst_asymmetry(rows) < 1e-9 * max_abs(rows.matrix),
          "the row-scaled formulation is symmetric");
    check(worst_asymmetry(cols) < 1e-9 * max_abs(cols.matrix),
          "the scaled-potential formulation is symmetric");
    // And the natural one is NOT -- without this the symmetry checks could
    // pass on an assembly that symmetrised everything by accident.
    check(worst_asymmetry(natural) > 1e-6 * scale,
          "the natural formulation is visibly NOT symmetric, so the two above are real");
}

void test_dc_refuses_the_scaled_formulations() {
    const auto refuses = [](Conditioning f) {
        try {
            formulation_scales(f, 0.0);
        } catch (const std::invalid_argument&) {
            return true;
        } catch (...) {
            return false;
        }
        return false;
    };
    check(refuses(Conditioning::RowScaled), "row scaling refuses DC rather than dividing by zero");
    check(refuses(Conditioning::ScaledPhi), "the scaled potential refuses DC too");

    bool natural_ok = true;
    try {
        formulation_scales(Conditioning::Natural, 0.0);
    } catch (...) {
        natural_ok = false;
    }
    check(natural_ok, "the natural formulation works at any frequency");

    // The symmetry condition, asserted rather than claimed.
    const double omega = 12345.0;
    const Complex jw(0.0, omega);
    for (const Conditioning f : {Conditioning::RowScaled, Conditioning::ScaledPhi}) {
        const FormulationScales s = formulation_scales(f, omega);
        check(std::abs(s.column - s.row * jw) < 1e-9 * std::abs(s.column),
              "c == r * j*omega for a symmetric formulation");
        check(s.symmetric, "and it says so");
    }
    const FormulationScales n = formulation_scales(Conditioning::Natural, omega);
    check(std::abs(n.column - n.row * jw) > 1.0, "and NOT for the natural one");
    check(!n.symmetric, "which it also says");
}

// The real mesh, and the first time a full-size system exists.
void test_cylinder() {
#ifdef APHI_MESH_DIR
    const std::string dir = APHI_MESH_DIR;
    Mesh m;
    try {
        m = read_gmsh_msh(dir + "/cylinder_box.msh");
    } catch (const std::exception& e) {
        check(false, std::string("cylinder_box.msh unreadable: ") + e.what());
        return;
    }
    scale_mesh_to_metres(m, LengthUnit::Millimetre);

    Problem p;
    p.type = AnalysisType::DC;
    p.length_unit = LengthUnit::Millimetre;
    Body wire;
    wire.name = "B1";
    wire.volume = "wire";
    wire.sigma = 5.8e7;
    wire.line = 1;
    Body air;
    air.name = "B2";
    air.volume = "air";
    air.sigma = 0.0;
    air.line = 2;
    p.bodies = {wire, air};
    Port p1;
    p1.name = "P1";
    p1.type = PortType::BoundaryCurrent;
    p1.surface = {"wire_bottom"};
    p1.amplitude = 1.0;
    p1.line = 3;
    Port p2;
    p2.name = "P2";
    p2.type = PortType::BoundaryVoltage;
    p2.surface = {"wire_top"};
    p2.amplitude = 0.0;
    p2.line = 4;
    p.ports = {p1, p2};

    const BoundProblem b = bind_to_mesh(p, m);
    const DofMap d = build_dof_map(b, m);
    const SparsityPattern sp = build_sparsity(d, b, m);
    const AssembledSystem sys = assemble(b, m, d, sp, 0.0, Conditioning::Natural);

    check(sys.matrix.nnz() == sp.nnz(), "the assembled matrix uses every slot the pattern gave");
    check(max_abs(sys.matrix) > 0.0, "and it is not all zeros");

    // No row may be empty of actual values: a row of zeros is a singular
    // matrix with no other symptom.
    int empty_rows = 0;
    for (int r = 0; r < d.num_total; ++r) {
        double row_max = 0.0;
        for (int k = sys.matrix.row_ptr()[static_cast<std::size_t>(r)];
             k < sys.matrix.row_ptr()[static_cast<std::size_t>(r) + 1]; ++k) {
            row_max = std::max(row_max, std::abs(sys.matrix.values()[static_cast<std::size_t>(k)]));
        }
        if (row_max == 0.0) ++empty_rows;
    }
    check(empty_rows == 0, "no row of the assembled matrix is entirely zero");

    // The A block is symmetric at DC (curl-curl alone), which is a check on
    // the edge orientation signs across the whole mesh -- get one wrong and
    // the two triangles disagree.
    double worst = 0.0;
    for (int r = 0; r < d.num_a; ++r) {
        for (int k = sys.matrix.row_ptr()[static_cast<std::size_t>(r)];
             k < sys.matrix.row_ptr()[static_cast<std::size_t>(r) + 1]; ++k) {
            const int col = sys.matrix.col_index()[static_cast<std::size_t>(k)];
            if (col >= d.num_a) continue;
            worst = std::max(worst, std::abs(sys.matrix.values()[static_cast<std::size_t>(k)] -
                                             entry(sys.matrix, col, r)));
        }
    }
    check(worst < 1e-9 * max_abs(sys.matrix),
          "the A block is symmetric across the whole mesh (worst " + std::to_string(worst) + ")");

    std::cout << "  cylinder: " << sys.matrix.rows() << " x " << sys.matrix.rows() << ", "
              << sys.matrix.nnz() << " nonzeros, largest |entry| " << max_abs(sys.matrix) << "\n";
#endif
}

// The ScatterMap, validated against `find_slot` -- the independent lookup it
// exists to replace.
//
// Without this the map's failure modes are caught only because a wrong index
// happens to leave the value array and crash. A wrong index that landed IN
// bounds would give a silently wrong matrix, so every slot is checked
// against what the pattern says, through the map's public accessors.
void test_scatter_map_is_consistent() {
    Mesh m = make_cube();
    const BoundProblem b = bind_cube(m, 1e6);
    const DofMap d = build_dof_map(b, m);
    const SparsityPattern sp = build_sparsity(d, b, m);
    const ScatterMap map = ScatterMap::build(sp, d, m, b);

    check(map.num_tets() == m.num_tets(), "the map covers every tet");

    int checked = 0;
    bool counts_sane = true, positions_sane = true, slots_right = true, live_ascending = true;
    for (int t = 0; t < m.num_tets(); ++t) {
        const TetDofs td = d.local_dofs(t, m, b);
        const int n = map.live_count(t);
        if (n < 0 || n > 16) counts_sane = false;

        // Rebuild the tet's live list independently, the way build_sparsity
        // does, and require the map to agree about its size.
        std::vector<int> live;
        for (const DofEntry& e : td.edge) {
            if (e.index >= 0) live.push_back(e.index);
        }
        for (const DofEntry& e : td.phi) {
            if (e.index >= 0) live.push_back(e.index);
        }
        std::sort(live.begin(), live.end());
        live.erase(std::unique(live.begin(), live.end()), live.end());
        if (static_cast<int>(live.size()) != n) counts_sane = false;

        // Each local DOF's position must point at its own global index.
        for (int i = 0; i < 6; ++i) {
            const int at = map.edge_at(t)[i];
            if (td.edge[static_cast<std::size_t>(i)].index < 0) {
                if (at != -1) positions_sane = false;
            } else if (at < 0 || at >= n ||
                       live[static_cast<std::size_t>(at)] !=
                           td.edge[static_cast<std::size_t>(i)].index) {
                positions_sane = false;
            }
        }
        for (int i = 0; i < 10; ++i) {
            const int at = map.phi_at(t)[i];
            if (td.phi[static_cast<std::size_t>(i)].index < 0) {
                if (at != -1) positions_sane = false;
            } else if (at < 0 || at >= n ||
                       live[static_cast<std::size_t>(at)] !=
                           td.phi[static_cast<std::size_t>(i)].index) {
                positions_sane = false;
            }
        }

        // And every slot must be the one find_slot would return. This is the
        // check that makes a wrong index a failure rather than a crash.
        for (int r = 0; r < n; ++r) {
            for (int c = 0; c < n; ++c) {
                const int got = map.slots(t)[r * n + c];
                const int want = sp.find_slot(live[static_cast<std::size_t>(r)],
                                              live[static_cast<std::size_t>(c)]);
                if (got != want || got < 0 || got >= static_cast<int>(sp.nnz())) {
                    slots_right = false;
                }
                ++checked;
            }
        }
        for (int i = 1; i < n; ++i) {
            if (live[static_cast<std::size_t>(i)] <= live[static_cast<std::size_t>(i - 1)]) {
                live_ascending = false;
            }
        }
    }
    check(counts_sane, "live_count matches an independently built live list for every tet");
    check(positions_sane, "edge_at and phi_at point at each local DOF's own global index");
    check(live_ascending, "the live list is strictly ascending -- what lets the sweep resume");
    check(slots_right, "every one of the " + std::to_string(checked) +
                           " slots is the one find_slot returns");
}

// `make_system` + `refill` split, and the reusable `ScatterMap`. Together
// these are what makes a frequency sweep cheap, and both have a failure mode
// that produces a plausible-looking wrong matrix rather than a crash:
// `refill` forgetting to clear, and the two scatter paths disagreeing.
void test_refill_and_scatter_map() {
    Mesh m = make_cube();
    const double omega = 2.0 * M_PI * 1e6;
    const BoundProblem b = bind_cube(m, 1e6);
    const DofMap d = build_dof_map(b, m);
    const SparsityPattern sp = build_sparsity(d, b, m);

    const AssembledSystem once = assemble(b, m, d, sp, omega, Conditioning::Natural);

    // --- make_system on its own -----------------------------------------
    AssembledSystem sys = make_system(sp, d.num_total);
    check(sys.matrix.rows() == sp.rows && sys.matrix.nnz() == sp.nnz(),
          "make_system adopts the pattern's shape");
    check(max_abs(sys.matrix) == 0.0, "and starts at zero");
    bool rhs_zero = true;
    for (const Complex& v : sys.rhs) {
        if (v != Complex(0.0, 0.0)) rhs_zero = false;
    }
    check(rhs_zero && sys.rhs.size() == static_cast<std::size_t>(d.num_total),
          "with one zero RHS entry per unknown");

    // --- refill reproduces assemble --------------------------------------
    refill(sys, b, m, d, sp, omega, Conditioning::Natural);
    check(same_values(sys, once), "make_system + refill equals assemble exactly");

    // --- refill CLEARS ----------------------------------------------------
    // The scatter accumulates, so a refill that forgot to zero would double
    // every value. The result would still be symmetric, still have the right
    // pattern, and still scale correctly between formulations -- it would
    // look right to every other test here.
    refill(sys, b, m, d, sp, omega, Conditioning::Natural);
    check(same_values(sys, once), "a second refill at the same omega gives the same system, "
                                  "not twice it -- it clears first");

    // --- a sweep: no residue from the previous frequency -------------------
    const double other = 2.0 * M_PI * 1e9;
    refill(sys, b, m, d, sp, other, Conditioning::Natural);
    const AssembledSystem fresh = assemble(b, m, d, sp, other, Conditioning::Natural);
    check(same_values(sys, fresh),
          "refilling at 1 GHz after 1 MHz equals a fresh assembly at 1 GHz");
    check(!same_values(sys, once),
          "and differs from the 1 MHz system, so the frequency really did change it");

    // --- the ScatterMap path is the same computation ----------------------
    const ScatterMap map = ScatterMap::build(sp, d, m, b);
    check(map.num_tets() == m.num_tets(), "the map covers every tet");
    check(map.bytes() > 0, "and holds something");

    AssembledSystem mapped = make_system(sp, d.num_total);
    refill(mapped, b, m, d, sp, omega, Conditioning::Natural, &map);
    check(bitwise_equal(mapped, once),
          "the ScatterMap path is BITWISE identical to working the slots out per tet -- "
          "same slots, same order, same arithmetic");

    // Across all three formulations, since the map is shared between them.
    for (const Conditioning f : {Conditioning::RowScaled, Conditioning::ScaledPhi}) {
        AssembledSystem a = assemble(b, m, d, sp, omega, f);
        AssembledSystem c = make_system(sp, d.num_total);
        refill(c, b, m, d, sp, omega, f, &map);
        check(bitwise_equal(a, c), "one ScatterMap serves every formulation");
    }

    // --- the map really is frequency-independent ---------------------------
    // Rebuilding it is not a function of omega at all; assert that using one
    // built before any assembly still works at a different frequency.
    AssembledSystem swept = make_system(sp, d.num_total);
    refill(swept, b, m, d, sp, other, Conditioning::Natural, &map);
    check(bitwise_equal(swept, fresh), "and one map serves every frequency");

    // --- the guards -------------------------------------------------------
    // Passing a system or a map that belongs to a different problem would
    // write into slots that mean something else, so both are refused.
    //
    // Worth recording why this does not use the DC cube as the "other"
    // problem, which was the first attempt: the cube is a single CONDUCTOR,
    // so Phi lives on all of it at DC and at AC alike and the two patterns
    // are identical. `phi_on_conductors_only` only changes the structure
    // when there is a non-conducting body to exclude.
    SparsityPattern tiny;
    tiny.rows = 3;
    tiny.cols = 3;
    tiny.row_ptr = {0, 1, 2, 3};
    tiny.col_index = {0, 1, 2};

    bool refused_system = false;
    try {
        AssembledSystem wrong = make_system(tiny, 3);
        refill(wrong, b, m, d, sp, omega, Conditioning::Natural);
    } catch (const std::invalid_argument&) {
        refused_system = true;
    }
    check(refused_system, "a system made from another pattern is refused");

    bool refused_rhs = false;
    try {
        AssembledSystem short_rhs = make_system(sp, d.num_total - 1);
        refill(short_rhs, b, m, d, sp, omega, Conditioning::Natural);
    } catch (const std::invalid_argument&) {
        refused_rhs = true;
    }
    check(refused_rhs, "so is a RHS of the wrong length");

    bool refused_map = false;
    try {
        ScatterMap wrong_map = ScatterMap::build(sp, d, m, b);
        Mesh bigger = make_cube();
        bigger.tets.push_back(bigger.tets.front());  // one more tet, nothing else
        bigger.build_topology();
        refill(sys, b, bigger, d, sp, omega, Conditioning::Natural, &wrong_map);
    } catch (const std::invalid_argument&) {
        refused_map = true;
    }
    check(refused_map, "and a ScatterMap built for a different mesh");

    bool accepted_matching = true;
    try {
        AssembledSystem right = make_system(sp, d.num_total);
        refill(right, b, m, d, sp, omega, Conditioning::Natural);
    } catch (...) {
        accepted_matching = false;
    }
    check(accepted_matching, "and a matching one is not -- the guard is not simply always firing");
}

}  // namespace

int main() {
    test_from_pattern();
    test_agrees_with_the_triplet_path();
    test_dc_structure();
    test_formulations_are_row_and_column_scalings();
    test_symmetry();
    test_dc_refuses_the_scaled_formulations();
    test_scatter_map_is_consistent();
    test_refill_and_scatter_map();
    test_cylinder();

    std::cout << (g_checks - g_failures) << "/" << g_checks << " checks passed\n";
    return g_failures == 0 ? 0 : 1;
}
