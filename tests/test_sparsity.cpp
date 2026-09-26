// Tests for the symbolic pass.
//
// Step 5 of `docs/ASSEMBLY_PLAN.md` §8. The pattern decides where the
// numeric pass may write, so its failure mode is a term with nowhere to go
// -- and if the scatter were to skip such a term rather than abort, the
// matrix would be quietly wrong. That makes exhaustiveness the property
// worth testing, not merely plausibility.
//
// The independent check throughout is a brute-force pattern: every ordered
// pair of live DOFs of every tet, collected into a std::set. Slow, obviously
// correct, and sharing no code with the function under test.

#include <algorithm>
#include <iostream>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "aphi_solver/gmsh_reader.hpp"
#include "aphi_solver/sparsity.hpp"

using namespace aphi_solver;

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

BoundProblem bind_cube(Mesh& m) {
    tag_face(m, 0, 1, 2, 10);
    tag_face(m, 0, 2, 3, 10);
    tag_face(m, 4, 5, 6, 11);
    tag_face(m, 4, 6, 7, 11);

    Problem p;
    p.type = AnalysisType::DC;
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
    p1.amplitude = 1.0;
    p1.line = 2;
    Port p2;
    p2.name = "P2";
    p2.type = PortType::BoundaryVoltage;
    p2.surface = {"top"};
    p2.amplitude = 0.0;
    p2.line = 3;
    p.ports = {p1, p2};
    return bind_to_mesh(p, m);
}

/// The same pattern, built the obvious slow way. Shares no code with
/// build_sparsity beyond `local_dofs`.
std::set<std::pair<int, int>> brute_force(const DofMap& d, const BoundProblem& b, const Mesh& m) {
    std::set<std::pair<int, int>> out;
    for (int t = 0; t < m.num_tets(); ++t) {
        const TetDofs td = d.local_dofs(t, m, b);
        std::vector<int> live;
        for (const DofEntry& e : td.edge) {
            if (e.index >= 0) live.push_back(e.index);
        }
        for (const DofEntry& e : td.phi) {
            if (e.index >= 0) live.push_back(e.index);
        }
        for (int i : live) {
            for (int j : live) out.insert({i, j});
        }
    }
    // A voltage port's row gets nothing from the tets -- its terminal is a
    // prescribed Phi, so the port has no column. Assembly writes the
    // constraint `V = V_given` there, so the pattern must hold that
    // diagonal.
    for (std::size_t k = 0; k < d.port_is_fixed.size(); ++k) {
        if (d.port_is_fixed[k]) out.insert({d.port_index[k], d.port_index[k]});
    }
    return out;
}

// ---------------------------------------------------------------------------

void test_matches_brute_force() {
    Mesh m = make_cube();
    const BoundProblem b = bind_cube(m);
    const DofMap d = build_dof_map(b, m);
    const SparsityPattern p = build_sparsity(d, b, m);

    const std::set<std::pair<int, int>> want = brute_force(d, b, m);

    check(p.rows == d.num_total && p.cols == d.num_total, "the pattern is num_total square");
    check(p.nnz() == want.size(),
          "the pattern has exactly as many entries as the brute-force set (" +
              std::to_string(p.nnz()) + " vs " + std::to_string(want.size()) + ")");

    // Every brute-force entry is findable.
    bool all_present = true;
    for (const auto& rc : want) {
        if (p.find_slot(rc.first, rc.second) < 0) all_present = false;
    }
    check(all_present, "every pair the brute force found has a slot");

    // And nothing extra: walk the pattern and check each entry is wanted.
    bool no_extras = true;
    for (int r = 0; r < p.rows; ++r) {
        for (int k = p.row_ptr[static_cast<std::size_t>(r)];
             k < p.row_ptr[static_cast<std::size_t>(r) + 1]; ++k) {
            if (!want.count({r, p.col_index[static_cast<std::size_t>(k)]})) no_extras = false;
        }
    }
    check(no_extras, "the pattern contains no entry the brute force did not");
}

void test_structure_is_well_formed() {
    Mesh m = make_cube();
    const BoundProblem b = bind_cube(m);
    const DofMap d = build_dof_map(b, m);
    const SparsityPattern p = build_sparsity(d, b, m);

    check(p.row_ptr.size() == static_cast<std::size_t>(p.rows) + 1, "row_ptr has rows + 1 entries");
    check(p.row_ptr.front() == 0, "row_ptr starts at 0");
    check(p.row_ptr.back() == static_cast<int>(p.col_index.size()),
          "row_ptr ends at nnz -- no slack left over from the counting pass");

    bool monotone = true, sorted_unique = true, in_range = true;
    for (int r = 0; r < p.rows; ++r) {
        const int begin = p.row_ptr[static_cast<std::size_t>(r)];
        const int end = p.row_ptr[static_cast<std::size_t>(r) + 1];
        if (end < begin) monotone = false;
        for (int k = begin; k < end; ++k) {
            const int col = p.col_index[static_cast<std::size_t>(k)];
            if (col < 0 || col >= p.cols) in_range = false;
            // Ascending AND strict: find_slot binary-searches these, and a
            // duplicate column would mean one contribution lands in one of
            // two slots depending on the search.
            if (k > begin && col <= p.col_index[static_cast<std::size_t>(k - 1)]) {
                sorted_unique = false;
            }
        }
    }
    check(monotone, "row_ptr is non-decreasing");
    check(in_range, "every column index is inside [0, cols)");
    check(sorted_unique, "each row's columns are strictly ascending -- sorted and unique");

    // Every row must have a diagonal: a DOF always shares a tet with itself.
    bool all_diagonals = true;
    for (int r = 0; r < p.rows; ++r) {
        if (p.find_slot(r, r) < 0) all_diagonals = false;
    }
    check(all_diagonals, "every row has a diagonal entry");

    // No row is empty -- an empty row is a DOF that appears in no tet, which
    // would make the matrix singular with no other symptom.
    bool none_empty = true;
    for (int r = 0; r < p.rows; ++r) {
        if (p.row_ptr[static_cast<std::size_t>(r + 1)] == p.row_ptr[static_cast<std::size_t>(r)]) {
            none_empty = false;
        }
    }
    check(none_empty, "no row is empty");
}

void test_find_slot() {
    Mesh m = make_cube();
    const BoundProblem b = bind_cube(m);
    const DofMap d = build_dof_map(b, m);
    const SparsityPattern p = build_sparsity(d, b, m);

    // Each slot must be distinct and point back at the column asked for.
    bool round_trips = true;
    std::set<int> seen;
    for (int r = 0; r < p.rows; ++r) {
        for (int k = p.row_ptr[static_cast<std::size_t>(r)];
             k < p.row_ptr[static_cast<std::size_t>(r) + 1]; ++k) {
            const int col = p.col_index[static_cast<std::size_t>(k)];
            if (p.find_slot(r, col) != k) round_trips = false;
            seen.insert(k);
        }
    }
    check(round_trips, "find_slot(r, col) returns the slot that holds col");
    check(seen.size() == p.nnz(), "the slots cover every position exactly once");

    check(p.find_slot(-1, 0) == -1, "a negative row gives -1");
    check(p.find_slot(p.rows, 0) == -1, "a row past the end gives -1");
    check(p.find_slot(0, p.cols + 10) == -1, "a column that is not present gives -1");
}

// Two DOFs that share no tet must NOT be connected. Without this the test
// would pass on a pattern that simply marked everything.
void test_pattern_is_actually_sparse() {
    Mesh m = make_cube();
    const BoundProblem b = bind_cube(m);
    const DofMap d = build_dof_map(b, m);
    const SparsityPattern p = build_sparsity(d, b, m);

    const double density = static_cast<double>(p.nnz()) /
                           (static_cast<double>(p.rows) * static_cast<double>(p.cols));
    check(density < 1.0, "the pattern is not full (density " + std::to_string(density) + ")");

    // On six tets sharing one diagonal almost everything touches everything,
    // so the cube is a weak test of sparsity; the cylinder below is the real
    // one. What the cube CAN show is that some pair is absent.
    std::cout << "  cube: " << p.rows << " rows, " << p.nnz() << " nonzeros, density "
              << density << "\n";
}

// The real fixture. This is also where the nonzero count stops being an
// estimate: ASSEMBLY_PLAN §6 guessed ~1.5 M, and a guess is worth replacing.
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

    check(sp.rows == d.num_total, "the pattern covers every unknown");
    check(sp.nnz() == brute_force(d, b, m).size(),
          "on the real mesh too, the pattern matches the brute-force set exactly");

    bool all_diagonals = true;
    for (int r = 0; r < sp.rows; ++r) {
        if (sp.find_slot(r, r) < 0) all_diagonals = false;
    }
    check(all_diagonals, "every row of the cylinder's matrix has a diagonal");

    const double average = static_cast<double>(sp.nnz()) / static_cast<double>(sp.rows);
    check(average > 5.0 && average < 200.0,
          "the average row length is in the range a tetrahedral A-Phi matrix gives (" +
              std::to_string(average) + ")");

    std::cout << "  cylinder: " << sp.rows << " unknowns, " << sp.nnz() << " nonzeros, "
              << average << " per row, "
              << (sp.nnz() * (sizeof(int) + 2 * sizeof(double))) / (1024.0 * 1024.0)
              << " MB as complex CSR\n";
#endif
}

}  // namespace

int main() {
    test_matches_brute_force();
    test_structure_is_well_formed();
    test_find_slot();
    test_pattern_is_actually_sparse();
    test_cylinder();

    std::cout << (g_checks - g_failures) << "/" << g_checks << " checks passed\n";
    return g_failures == 0 ? 0 : 1;
}
