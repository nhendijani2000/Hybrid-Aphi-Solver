// Tests for nodal potential output -- step 1 of `docs/POSTPROCESSING_PLAN.md` §6.
//
// Two things are being checked, and they fail in different ways.
//
// The POTENTIAL has to come back out of the solution vector correctly, which
// means honouring five different reasons a node might have the value it has,
// and applying `Phi = j*omega*Phi'` exactly once under ScaledPhi. Getting that
// wrong gives a smooth, plausible, wrong field.
//
// The WRITER has to be fast and lossless. `std::to_chars` is used for both
// reasons: it is much faster than a stream, and its shortest representation
// round-trips exactly, so the file is smaller AND loses nothing. A test reads
// every number back and requires bitwise equality -- a fixed precision would
// fail that.

#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdio>
#include <fstream>
#include <iterator>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "aphi_solver/basis_functions.hpp"
#include "aphi_solver/constants.hpp"
#include "aphi_solver/postprocess.hpp"

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

std::string temp_path(const char* name) {
    const char* dir = std::getenv("TEMP");
    return (dir != nullptr ? std::string(dir) : std::string(".")) + "/aphi_test_" + name;
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

BoundProblem bind_cube(Mesh& m, double frequency) {
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
    p1.amplitude = 1.0;
    p1.line = 2;
    Port p2;
    p2.name = "P2";
    p2.type = PortType::BoundaryVoltage;
    p2.surface = {"top"};
    p2.amplitude = 2.5;
    p2.line = 3;
    p.ports = {p1, p2};
    return bind_to_mesh(p, m);
}

// ---------------------------------------------------------------------------

void test_node_layout_and_positions() {
    Mesh m = make_cube();
    const BoundProblem b = bind_cube(m, 1e6);
    const DofMap d = build_dof_map(b, m);

    Solution s;
    s.x.assign(static_cast<std::size_t>(d.num_total), Complex(0.0, 0.0));
    s.conditioning = Conditioning::Natural;
    s.omega = 2.0 * M_PI * 1e6;

    const NodalPotential p = potential_at_nodes(m, b, d, s);

    check(p.size() == d.num_p2_nodes, "one entry per P2 node");
    check(p.num_vertices == m.num_nodes() && p.num_edges == m.num_edges(),
          "split into vertices then edge midpoints");
    check(p.size() == m.num_nodes() + m.num_edges(), "and they account for all of it");

    // Vertices keep their own coordinates; midpoints are the mean of the two
    // endpoints. An off-by-one in the numbering would scramble these.
    bool vertices_ok = true;
    for (int v = 0; v < m.num_nodes(); ++v) {
        const Vec3& a = m.nodes[static_cast<std::size_t>(v)];
        const Vec3& g = p.position[static_cast<std::size_t>(v)];
        if (a.x != g.x || a.y != g.y || a.z != g.z) vertices_ok = false;
    }
    check(vertices_ok, "every vertex keeps its own position");

    bool midpoints_ok = true;
    for (int e = 0; e < m.num_edges(); ++e) {
        const std::pair<int, int>& ends = m.edges[static_cast<std::size_t>(e)];
        const Vec3& x = m.nodes[static_cast<std::size_t>(ends.first)];
        const Vec3& y = m.nodes[static_cast<std::size_t>(ends.second)];
        const Vec3& g = p.position[static_cast<std::size_t>(m.num_nodes() + e)];
        if (std::abs(g.x - 0.5 * (x.x + y.x)) > 1e-15 ||
            std::abs(g.y - 0.5 * (x.y + y.y)) > 1e-15 ||
            std::abs(g.z - 0.5 * (x.z + y.z)) > 1e-15) {
            midpoints_ok = false;
        }
    }
    check(midpoints_ok, "and every edge node sits at its edge's midpoint");
}

void test_every_status_is_honoured() {
    Mesh m = make_cube();
    const BoundProblem b = bind_cube(m, 1e6);
    const DofMap d = build_dof_map(b, m);

    // A solution with a distinct value per index, so a misrouted read shows up
    // as the wrong number rather than as plausible noise.
    Solution s;
    s.conditioning = Conditioning::Natural;
    s.omega = 2.0 * M_PI * 1e6;
    s.x.resize(static_cast<std::size_t>(d.num_total));
    for (int i = 0; i < d.num_total; ++i) {
        s.x[static_cast<std::size_t>(i)] = Complex(100.0 + i, -0.5 * i);
    }

    const NodalPotential p = potential_at_nodes(m, b, d, s);

    int free_nodes = 0, terminals = 0, prescribed = 0, absent = 0, cuts = 0;
    bool free_ok = true, terminal_ok = true, fixed_ok = true, absent_ok = true;
    for (int i = 0; i < d.num_p2_nodes; ++i) {
        const std::size_t u = static_cast<std::size_t>(i);
        switch (d.phi_state[u]) {
            case PhiDof::Free:
                ++free_nodes;
                if (p.status[u] != PhiStatus::Free) free_ok = false;
                if (p.value[u] != s.x[static_cast<std::size_t>(d.phi_index[u])]) free_ok = false;
                break;
            case PhiDof::Port: {
                ++terminals;
                if (p.status[u] != PhiStatus::Terminal) terminal_ok = false;
                const std::size_t k = static_cast<std::size_t>(d.phi_port[u]);
                const Complex want = d.port_is_fixed[k]
                                         ? d.port_value[k]
                                         : s.x[static_cast<std::size_t>(d.port_index[k])];
                if (p.value[u] != want) terminal_ok = false;
                break;
            }
            case PhiDof::Cut:
                ++cuts;
                break;
            case PhiDof::Fixed:
                ++prescribed;
                if (p.status[u] != PhiStatus::Fixed) fixed_ok = false;
                if (p.value[u] != d.phi_fixed_value[u]) fixed_ok = false;
                break;
            case PhiDof::Absent:
                ++absent;
                if (p.status[u] != PhiStatus::Absent) absent_ok = false;
                if (p.value[u] != Complex(0.0, 0.0)) absent_ok = false;
                break;
        }
    }

    check(free_nodes > 0 && terminals > 0,
          "the fixture exercises free and terminal nodes (" + std::to_string(free_nodes) +
              " free, " + std::to_string(terminals) + " terminal)");
    check(free_ok, "a free node reads its own unknown");
    check(terminal_ok,
          "a terminal node reads its PORT's value -- every node of one terminal shares it");
    check(fixed_ok, "a prescribed node reads its prescribed value");
    check(absent_ok, "an absent node reads zero and says it is absent");

    // A voltage port is prescribed, so its terminal must carry exactly the
    // value from the input file and nothing from the solve.
    bool carries_prescribed = false;
    for (int i = 0; i < d.num_p2_nodes; ++i) {
        const std::size_t u = static_cast<std::size_t>(i);
        if (p.status[u] == PhiStatus::Terminal && p.value[u] == Complex(2.5, 0.0)) {
            carries_prescribed = true;
        }
    }
    check(carries_prescribed,
          "the voltage port's terminal carries exactly the 2.5 V from the input file");
}

// The one that catches the mistake SOLVER_PLAN §8 lists as a control.
void test_scaled_phi_is_multiplied_by_j_omega() {
    Mesh m = make_cube();
    const double f = 1e6;
    const double omega = 2.0 * M_PI * f;
    const BoundProblem b = bind_cube(m, f);
    const DofMap d = build_dof_map(b, m);

    Solution natural;
    natural.conditioning = Conditioning::Natural;
    natural.omega = omega;
    natural.x.resize(static_cast<std::size_t>(d.num_total));
    for (int i = 0; i < d.num_total; ++i) {
        natural.x[static_cast<std::size_t>(i)] = Complex(1.0 + 0.01 * i, 0.5 - 0.002 * i);
    }

    // Under ScaledPhi the stored unknown is Phi' = Phi/(j*omega), so the same
    // physical field is stored as x / (j*omega).
    Solution scaled = natural;
    scaled.conditioning = Conditioning::ScaledPhi;
    for (int i = 0; i < d.num_total; ++i) {
        scaled.x[static_cast<std::size_t>(i)] /= Complex(0.0, omega);
    }

    check(natural.phi_scale() == Complex(1.0, 0.0), "Natural does not rescale Phi");
    check(scaled.phi_scale() == Complex(0.0, omega), "ScaledPhi rescales by j*omega");

    const NodalPotential pn = potential_at_nodes(m, b, d, natural);
    const NodalPotential ps = potential_at_nodes(m, b, d, scaled);

    double worst = 0.0, scale = 0.0;
    for (int i = 0; i < pn.size(); ++i) {
        const std::size_t u = static_cast<std::size_t>(i);
        if (pn.status[u] == PhiStatus::Absent || pn.status[u] == PhiStatus::Fixed) continue;
        scale = std::max(scale, std::abs(pn.value[u]));
        worst = std::max(worst, std::abs(pn.value[u] - ps.value[u]));
    }
    check(scale > 0.0, "the comparison covers some non-trivial nodes");
    check(worst < 1e-9 * scale,
          "the same physical field stored two ways reads back the same (worst " +
              std::to_string(worst) + " against " + std::to_string(scale) +
              ") -- forgetting the j*omega would be a factor of " + std::to_string(omega));

    // A prescribed value was never scaled by the substitution, so it must not be
    // scaled back either. This is the half of the rule that is easy to miss.
    bool prescribed_unscaled = true;
    for (int i = 0; i < ps.size(); ++i) {
        const std::size_t u = static_cast<std::size_t>(i);
        if (ps.status[u] != PhiStatus::Terminal) continue;
        const std::size_t k = static_cast<std::size_t>(d.phi_port[u]);
        if (d.port_is_fixed[k] && ps.value[u] != d.port_value[k]) prescribed_unscaled = false;
    }
    check(prescribed_unscaled,
          "and a prescribed terminal keeps its physical value under ScaledPhi -- it was never "
          "substituted, so it must not be un-substituted");
}

void test_round_trip_is_lossless() {
    Mesh m = make_cube();
    const BoundProblem b = bind_cube(m, 1e6);
    const DofMap d = build_dof_map(b, m);
    Solution s;
    s.conditioning = Conditioning::RowScaled;
    s.omega = 2.0 * M_PI * 1e6;
    s.x.resize(static_cast<std::size_t>(d.num_total));
    for (int i = 0; i < d.num_total; ++i) {
        // Awkward values on purpose: a fixed-precision writer loses these.
        s.x[static_cast<std::size_t>(i)] =
            Complex(1.0 / (3.0 + i) * 1e-17, -std::sqrt(2.0) * (i + 1) * 1e13);
    }
    const NodalPotential p = potential_at_nodes(m, b, d, s);

    const std::string path = temp_path("potential.out");
    const WriteStats stats = write_potential(path, p, s);
    check(stats.nodes == p.size(), "the writer reports every node");
    check(stats.bytes > 0, "and some bytes");

    std::ifstream in(path);
    check(in.good(), "the file can be reopened");
    std::string line;
    int header = 0, read_back = 0;
    bool exact = true, statuses_ok = true;
    while (std::getline(in, line)) {
        if (!line.empty() && line[0] == '#') {
            ++header;
            continue;
        }
        std::istringstream row(line);
        int index = 0;
        double x = 0, y = 0, z = 0, re = 0, im = 0;
        char status = 0;
        row >> index >> x >> y >> z >> re >> im >> status;
        const std::size_t u = static_cast<std::size_t>(index);
        // BITWISE: to_chars's shortest form round-trips exactly, so anything
        // less than equality here means the writer is lossy.
        if (re != p.value[u].real() || im != p.value[u].imag()) exact = false;
        if (x != p.position[u].x || y != p.position[u].y || z != p.position[u].z) exact = false;
        if (std::string("FTCPA").find(status) == std::string::npos) statuses_ok = false;
        ++read_back;
    }
    check(header >= 5, "the file carries a header describing itself");
    check(read_back == p.size(),
          "every node is in the file (" + std::to_string(read_back) + " of " +
              std::to_string(p.size()) + ")");
    check(exact,
          "and every value reads back BITWISE identical -- to_chars's shortest form round-trips, "
          "where a fixed precision would not");
    check(statuses_ok, "with a recognised status letter on every line");
    std::remove(path.c_str());
}

// Buffering, measured rather than asserted.
void test_buffered_writing_is_faster() {
    const int n = 200000;
    std::vector<double> values(static_cast<std::size_t>(n));
    for (int i = 0; i < n; ++i) {
        values[static_cast<std::size_t>(i)] = std::sqrt(2.0) * (i + 1) / 7.0;
    }

    const std::string fast_path = temp_path("fast.txt");
    const std::string slow_path = temp_path("slow.txt");

    auto t0 = std::chrono::steady_clock::now();
    {
        std::FILE* f = std::fopen(fast_path.c_str(), "wb");
        TextBuffer buf(f);
        for (int i = 0; i < n; ++i) {
            buf.put(i);
            buf.put(' ');
            buf.put(values[static_cast<std::size_t>(i)]);
            buf.put('\n');
        }
        buf.flush();
        std::fclose(f);
    }
    auto t1 = std::chrono::steady_clock::now();

    // The obvious way, for comparison: a stream, one insertion per value.
    {
        std::ofstream out(slow_path);
        out << std::setprecision(17);
        for (int i = 0; i < n; ++i) {
            out << i << ' ' << values[static_cast<std::size_t>(i)] << '\n';
        }
    }
    auto t2 = std::chrono::steady_clock::now();

    const double fast = std::chrono::duration<double, std::milli>(t1 - t0).count();
    const double slow = std::chrono::duration<double, std::milli>(t2 - t1).count();
    std::cout << "  " << n << " values: TextBuffer " << fast << " ms, ofstream << " << slow
              << " ms  ->  " << slow / fast << "x\n";

    check(fast < slow,
          "the buffered writer beats a stream (" + std::to_string(fast) + " ms against " +
              std::to_string(slow) + " ms)");

    // And it is smaller, because the shortest round-tripping form beats
    // setprecision(17) on length as well as on accuracy.
    std::ifstream a(fast_path, std::ios::binary | std::ios::ate);
    std::ifstream b(slow_path, std::ios::binary | std::ios::ate);
    const long long fast_bytes = a.tellg();
    const long long slow_bytes = b.tellg();
    std::cout << "  bytes: " << fast_bytes << " against " << slow_bytes << "\n";
    check(fast_bytes < slow_bytes,
          "and produces a smaller file while losing nothing (" + std::to_string(fast_bytes) +
              " against " + std::to_string(slow_bytes) + " bytes)");

    std::remove(fast_path.c_str());
    std::remove(slow_path.c_str());
}

void test_buffer_handles_awkward_sizes() {
    // A capacity smaller than one line, so make_room has to flush mid-line and
    // then grow. An off-by-one in the buffer bookkeeping shows up here and
    // nowhere else.
    const std::string path = temp_path("tiny.txt");
    {
        std::FILE* f = std::fopen(path.c_str(), "wb");
        // A genuinely tiny buffer, so `make_room` hits its boundary on almost
        // every write. With a comfortable buffer an off-by-one there is reached
        // once in thousands of calls and a test can pass through luck.
        TextBuffer buf(f, 17);
        for (int i = 0; i < 1000; ++i) {
            buf.put(i);
            buf.put(' ');
            buf.put(1.0 / (i + 1));
            buf.put('\n');
        }
        buf.flush();
        std::fclose(f);
    }
    std::ifstream in(path);
    int lines = 0;
    std::string line;
    bool values_ok = true;
    while (std::getline(in, line)) {
        std::istringstream row(line);
        int i = 0;
        double v = 0;
        row >> i >> v;
        if (v != 1.0 / (i + 1)) values_ok = false;
        ++lines;
    }
    check(lines == 1000, "a small buffer still writes every line (" + std::to_string(lines) + ")");
    check(values_ok, "and every value survives the flushes intact");

    // Single characters against an awkward capacity: this is the pattern that
    // drives `used_` exactly to the end of the buffer, over and over.
    const std::string edge = temp_path("edge.txt");
    {
        std::FILE* f = std::fopen(edge.c_str(), "wb");
        TextBuffer buf(f, 16);
        for (int i = 0; i < 5000; ++i) buf.put(static_cast<char>('a' + (i % 26)));
        buf.flush();
        std::fclose(f);
    }
    std::ifstream ein(edge, std::ios::binary);
    std::string content((std::istreambuf_iterator<char>(ein)),
                        std::istreambuf_iterator<char>());
    bool chars_ok = content.size() == 5000;
    for (std::size_t i = 0; i < content.size() && chars_ok; ++i) {
        if (content[i] != static_cast<char>('a' + (i % 26))) chars_ok = false;
    }
    check(chars_ok,
          "5000 single characters through a 16-byte buffer come out exactly right (" +
              std::to_string(content.size()) + " bytes) -- this is what exercises the"
              " flush boundary");
    std::remove(edge.c_str());
    std::remove(path.c_str());
}


// The VTK node order, checked geometrically -- the only way to catch it.
//
// VTK's quadratic tet wants its mid-edge nodes in the order (0,1), (1,2), (0,2),
// (0,3), (1,3), (2,3); ours are numbered (0,1), (0,2), (0,3), (1,2), (2,3),
// (1,3). Four of six disagree. A wrong permutation still produces a file
// ParaView opens and renders smoothly, with the quadratic nodes silently
// attached to the wrong edges -- so this checks each slot's POSITION against the
// midpoint VTK expects to find there.
void test_vtk_quadratic_node_order() {
    Mesh m = make_cube();
    const BoundProblem b = bind_cube(m, 1e6);
    const DofMap d = build_dof_map(b, m);
    Solution s;
    s.conditioning = Conditioning::Natural;
    s.omega = 2.0 * M_PI * 1e6;
    s.x.assign(static_cast<std::size_t>(d.num_total), Complex(1.0, 0.0));
    const NodalPotential p = potential_at_nodes(m, b, d, s);

    // VTK's own definition, written out here rather than derived from ours, so
    // the two cannot agree by construction.
    const int vtk_pairs[6][2] = {{0, 1}, {1, 2}, {0, 2}, {0, 3}, {1, 3}, {2, 3}};

    bool positions_match = true;
    int checked = 0;
    for (int t = 0; t < m.num_tets(); ++t) {
        const std::size_t ut = static_cast<std::size_t>(t);
        for (int slot = 0; slot < 6; ++slot) {
            const int local_edge = kVtkQuadraticTetEdgeOrder[static_cast<std::size_t>(slot)];
            const int e = m.tet_edges[ut][static_cast<std::size_t>(local_edge)];
            const Vec3& got = p.position[static_cast<std::size_t>(m.num_nodes() + e)];

            const Vec3& ca = m.nodes[static_cast<std::size_t>(
                m.tets[ut][static_cast<std::size_t>(vtk_pairs[slot][0])])];
            const Vec3& cb = m.nodes[static_cast<std::size_t>(
                m.tets[ut][static_cast<std::size_t>(vtk_pairs[slot][1])])];
            const Vec3 want{0.5 * (ca.x + cb.x), 0.5 * (ca.y + cb.y), 0.5 * (ca.z + cb.z)};

            if (std::abs(got.x - want.x) > 1e-15 || std::abs(got.y - want.y) > 1e-15 ||
                std::abs(got.z - want.z) > 1e-15) {
                positions_match = false;
            }
            ++checked;
        }
    }
    check(positions_match,
          "every VTK quadratic slot holds the midpoint of the corner pair VTK expects there (" +
              std::to_string(checked) + " slots) -- a wrong permutation renders smoothly and is "
              "wrong");

    // The permutation is a permutation: all six of our edges used, once each.
    std::vector<int> seen(6, 0);
    for (int slot = 0; slot < 6; ++slot) {
        ++seen[static_cast<std::size_t>(kVtkQuadraticTetEdgeOrder[static_cast<std::size_t>(slot)])];
    }
    bool bijective = true;
    for (int c : seen) {
        if (c != 1) bijective = false;
    }
    check(bijective, "and it uses each of our six local edges exactly once");

    // The file parses as the structure it claims.
    const std::string path = temp_path("out.vtk");
    const WriteStats stats = write_vtk(path, m, p, s);
    check(stats.bytes > 0, "the VTK file is written");
    std::ifstream in(path);
    std::string line;
    int points = 0, cells = 0, cell_types = 0, wrong_type = 0, wrong_count = 0;
    while (std::getline(in, line)) {
        std::istringstream row(line);
        std::string word;
        row >> word;
        if (word == "POINTS") {
            row >> points;
        } else if (word == "CELLS") {
            row >> cells;
            for (int c = 0; c < cells; ++c) {
                std::getline(in, line);
                std::istringstream cell(line);
                int n = 0;
                cell >> n;
                if (n != 10) ++wrong_count;
            }
        } else if (word == "CELL_TYPES") {
            row >> cell_types;
            for (int c = 0; c < cell_types; ++c) {
                std::getline(in, line);
                if (line != "24") ++wrong_type;
            }
        }
    }
    check(points == p.size(), "with a point for every P2 node");
    check(cells == m.num_tets(), "a cell for every tet");
    check(wrong_count == 0, "each cell listing ten nodes");
    check(wrong_type == 0, "and each typed 24 (VTK_QUADRATIC_TETRA), so the mid-edge values are "
                           "used rather than discarded");
    std::remove(path.c_str());
}


// A DOF map in which nothing is eliminated: every edge and every P2 node is a
// free unknown, numbered [edges | P2 nodes].
//
// The real map always pins a spanning tree and the Dirichlet edges, so no real
// solution can carry an arbitrary A. That makes the exact-reproduction tests
// below impossible against it -- not because the reconstruction is wrong, but
// because the field being asked for is not in the constrained space. Building
// the map by hand separates the two questions: these tests ask only whether
// compute_fields evaluates the basis correctly, which is what they can answer.
DofMap all_free_dofs(const Mesh& m) {
    DofMap d;
    d.num_nodes = m.num_nodes();
    d.num_p2_nodes = m.num_nodes() + m.num_edges();
    d.num_a = m.num_edges();
    d.num_phi = d.num_p2_nodes;
    d.num_ports = 0;
    d.num_total = d.num_a + d.num_phi;

    d.edge_state.assign(static_cast<std::size_t>(m.num_edges()), EdgeDof::Free);
    d.edge_index.resize(static_cast<std::size_t>(m.num_edges()));
    for (int e = 0; e < m.num_edges(); ++e) d.edge_index[static_cast<std::size_t>(e)] = e;

    d.phi_state.assign(static_cast<std::size_t>(d.num_p2_nodes), PhiDof::Free);
    d.phi_index.resize(static_cast<std::size_t>(d.num_p2_nodes));
    d.phi_port.assign(static_cast<std::size_t>(d.num_p2_nodes), -1);
    d.phi_fixed_value.assign(static_cast<std::size_t>(d.num_p2_nodes), Complex(0.0, 0.0));
    for (int i = 0; i < d.num_p2_nodes; ++i) {
        d.phi_index[static_cast<std::size_t>(i)] = d.num_a + i;
    }
    return d;
}

Vec3 p2_node_position(const Mesh& m, int p2) {
    if (p2 < m.num_nodes()) return m.nodes[static_cast<std::size_t>(p2)];
    const std::pair<int, int>& e = m.edges[static_cast<std::size_t>(p2 - m.num_nodes())];
    const Vec3& a = m.nodes[static_cast<std::size_t>(e.first)];
    const Vec3& b = m.nodes[static_cast<std::size_t>(e.second)];
    return Vec3{0.5 * (a.x + b.x), 0.5 * (a.y + b.y), 0.5 * (a.z + b.z)};
}

// The one structural invariant the edge reconstruction rests on: the global
// basis function is the local one times the tet's edge sign. compute_fields
// uses the global form with the RAW unknown, while assembly uses the local form
// and multiplies by the same sign when it scatters. If those two ever disagree,
// every reconstructed A is negated on roughly 58 % of edges and no residual
// would show it.
void test_edge_sign_convention_matches_assembly() {
    Mesh m = make_cube();
    bool value_ok = true, curl_ok = true;
    int flipped = 0, total = 0;
    const std::array<double, 4> L = {0.17, 0.31, 0.29, 0.23};

    for (int t = 0; t < m.num_tets(); ++t) {
        const TetGeometry g = compute_tet_geometry(m, t);
        for (int le = 0; le < 6; ++le) {
            const double s = static_cast<double>(
                m.tet_edge_signs[static_cast<std::size_t>(t)][static_cast<std::size_t>(le)]);
            const Vec3 vg = whitney_edge_value_global(m, t, g, le, L);
            const Vec3 vl = whitney_edge_value(g, le, L);
            const Vec3 cg = whitney_edge_curl_global(m, t, g, le);
            const Vec3 cl = whitney_edge_curl(g, le);
            if (std::abs(vg.x - s * vl.x) > 1e-12 || std::abs(vg.y - s * vl.y) > 1e-12 ||
                std::abs(vg.z - s * vl.z) > 1e-12) {
                value_ok = false;
            }
            if (std::abs(cg.x - s * cl.x) > 1e-12 || std::abs(cg.y - s * cl.y) > 1e-12 ||
                std::abs(cg.z - s * cl.z) > 1e-12) {
                curl_ok = false;
            }
            if (s < 0.0) ++flipped;
            ++total;
        }
    }
    check(value_ok, "the global Whitney value is the local one times the tet edge sign");
    check(curl_ok, "and so is the curl");
    check(flipped > 0, "and the sign is actually negative somewhere (" + std::to_string(flipped) +
                           " of " + std::to_string(total) +
                           " pairs) -- otherwise this test proves nothing");
}

// B = curl A, reproduced exactly.
//
// A(r) = (1/2) B x r is linear, so its edge integrals are represented exactly by
// the Whitney space, and curl A = B exactly. Setting each unknown to the edge
// line integral must therefore give back B in every tet, to round-off. This is
// the end-to-end check of the edge reconstruction including the orientation.
void test_uniform_b_is_reproduced_exactly() {
    Mesh m = make_cube();
    const BoundProblem bnd = bind_cube(m, 1e6);
    const DofMap d = all_free_dofs(m);

    const Vec3 B{0.3, -0.7, 1.1};
    Solution s;
    s.conditioning = Conditioning::Natural;
    s.omega = 2.0 * M_PI * 1e6;
    s.x.assign(static_cast<std::size_t>(d.num_total), Complex(0.0, 0.0));

    for (int e = 0; e < m.num_edges(); ++e) {
        const std::pair<int, int>& ends = m.edges[static_cast<std::size_t>(e)];
        const Vec3& p = m.nodes[static_cast<std::size_t>(ends.first)];
        const Vec3& q = m.nodes[static_cast<std::size_t>(ends.second)];
        const Vec3 mid{0.5 * (p.x + q.x), 0.5 * (p.y + q.y), 0.5 * (p.z + q.z)};
        // A = (1/2) B x r, evaluated at the midpoint: exact for a linear A.
        const Vec3 A{0.5 * (B.y * mid.z - B.z * mid.y), 0.5 * (B.z * mid.x - B.x * mid.z),
                     0.5 * (B.x * mid.y - B.y * mid.x)};
        const Vec3 dl{q.x - p.x, q.y - p.y, q.z - p.z};
        s.x[static_cast<std::size_t>(e)] = Complex(A.x * dl.x + A.y * dl.y + A.z * dl.z, 0.0);
    }

    const FieldOutput f = compute_fields(m, bnd, d, s);
    double worst = 0.0;
    for (int t = 0; t < m.num_tets(); ++t) {
        const Vec3C& bt = f.b_tet[static_cast<std::size_t>(t)];
        worst = std::max(worst, std::abs(bt[0].real() - B.x));
        worst = std::max(worst, std::abs(bt[1].real() - B.y));
        worst = std::max(worst, std::abs(bt[2].real() - B.z));
        worst = std::max(worst, std::abs(bt[0].imag()));
        worst = std::max(worst, std::abs(bt[1].imag()));
        worst = std::max(worst, std::abs(bt[2].imag()));
    }
    check(worst < 1e-12, "a uniform B is reproduced exactly in every tet from A = (1/2) B x r");

    double worst_v = 0.0;
    for (int v = 0; v < m.num_nodes(); ++v) {
        const Vec3C& bv = f.b_node[static_cast<std::size_t>(v)];
        worst_v = std::max(worst_v, std::abs(bv[0].real() - B.x));
        worst_v = std::max(worst_v, std::abs(bv[2].real() - B.z));
    }
    check(worst_v < 1e-12, "and the volume-weighted per-vertex B agrees, since averaging a "
                           "constant cannot change it");
}

// E = -grad(Phi), reproduced exactly for a linear Phi with A = 0.
//
// A linear function lies in the P2 space, so its gradient is exact. With every
// A unknown zero, E must be exactly -grad(Phi) -- a constant. This is the check
// that the P2 gradient reconstruction and the sign of E are right.
void test_linear_phi_gives_exact_uniform_e() {
    Mesh m = make_cube();
    const BoundProblem bnd = bind_cube(m, 1e6);
    const DofMap d = all_free_dofs(m);

    const Vec3 c{2.5, -1.25, 0.75};
    Solution s;
    s.conditioning = Conditioning::Natural;
    s.omega = 2.0 * M_PI * 1e6;
    s.x.assign(static_cast<std::size_t>(d.num_total), Complex(0.0, 0.0));
    for (int i = 0; i < d.num_p2_nodes; ++i) {
        const Vec3 p = p2_node_position(m, i);
        s.x[static_cast<std::size_t>(d.num_a + i)] =
            Complex(c.x * p.x + c.y * p.y + c.z * p.z + 4.0, 0.0);
    }

    const FieldOutput f = compute_fields(m, bnd, d, s);
    double worst = 0.0;
    for (int t = 0; t < m.num_tets(); ++t) {
        const Vec3C& et = f.e_tet[static_cast<std::size_t>(t)];
        worst = std::max(worst, std::abs(et[0].real() + c.x));
        worst = std::max(worst, std::abs(et[1].real() + c.y));
        worst = std::max(worst, std::abs(et[2].real() + c.z));
    }
    check(worst < 1e-12,
          "with A = 0 a linear Phi gives E = -grad(Phi) exactly, constant over every tet");

    double worst_v = 0.0;
    for (int v = 0; v < m.num_nodes(); ++v) {
        worst_v =
            std::max(worst_v, std::abs(f.e_node[static_cast<std::size_t>(v)][1].real() + c.y));
    }
    check(worst_v < 1e-12, "and at the vertices too");

    // The constant offset must not appear anywhere: E depends on Phi only
    // through its gradient.
    double worst_b = 0.0;
    for (int t = 0; t < m.num_tets(); ++t) {
        for (int i = 0; i < 3; ++i) {
            worst_b = std::max(worst_b, std::abs(f.b_tet[static_cast<std::size_t>(t)]
                                                        [static_cast<std::size_t>(i)]));
        }
    }
    check(worst_b < 1e-14, "and B stays exactly zero, since no A unknown was set");
}

// The per-vertex averaging, checked against an independent recomputation.
void test_vertex_averaging_is_volume_weighted() {
    Mesh m = make_cube();
    const BoundProblem bnd = bind_cube(m, 1e6);
    const DofMap d = all_free_dofs(m);

    Solution s;
    s.conditioning = Conditioning::Natural;
    s.omega = 2.0 * M_PI * 1e6;
    s.x.resize(static_cast<std::size_t>(d.num_total));
    for (int i = 0; i < d.num_total; ++i) {
        s.x[static_cast<std::size_t>(i)] = Complex(std::sin(0.7 * i + 0.3), std::cos(0.4 * i));
    }
    const FieldOutput f = compute_fields(m, bnd, d, s);

    // Total accumulated weight must be exactly four times the mesh volume:
    // every tet gives its volume to each of its four vertices. This catches a
    // signed volume leaking through, which would silently cancel.
    double total_w = 0.0, total_vol = 0.0;
    for (int v = 0; v < m.num_nodes(); ++v) total_w += f.vertex_weight[static_cast<std::size_t>(v)];
    for (int t = 0; t < m.num_tets(); ++t) {
        total_vol += std::abs(compute_tet_geometry(m, t).volume);
    }
    check(std::abs(total_w - 4.0 * total_vol) < 1e-12 * total_vol,
          "the accumulated vertex weight is exactly 4x the mesh volume");

    // B is constant per tet, so the per-vertex value is a plain volume-weighted
    // mean of the incident tets -- recomputed here from b_tet alone.
    std::vector<Complex> num(static_cast<std::size_t>(m.num_nodes()), Complex(0.0, 0.0));
    std::vector<double> den(static_cast<std::size_t>(m.num_nodes()), 0.0);
    for (int t = 0; t < m.num_tets(); ++t) {
        const double vol = std::abs(compute_tet_geometry(m, t).volume);
        for (int k = 0; k < 4; ++k) {
            const std::size_t v = static_cast<std::size_t>(
                m.tets[static_cast<std::size_t>(t)][static_cast<std::size_t>(k)]);
            num[v] += vol * f.b_tet[static_cast<std::size_t>(t)][2];
            den[v] += vol;
        }
    }
    double worst = 0.0;
    for (int v = 0; v < m.num_nodes(); ++v) {
        const std::size_t uv = static_cast<std::size_t>(v);
        if (den[uv] <= 0.0) continue;
        worst = std::max(worst, std::abs(num[uv] / den[uv] - f.b_node[uv][2]));
    }
    check(worst < 1e-14, "and b_node at a vertex is that mean, recomputed independently from b_tet");
    check(f.num_orphan_vertices == 0, "the cube has no orphan vertices");
}

// Phi is copied, not rebuilt: compute_fields must agree with potential_at_nodes
// exactly, or the VTK file and the field output would disagree about the same
// quantity.
void test_phi_matches_potential_at_nodes() {
    Mesh m = make_cube();
    const BoundProblem bnd = bind_cube(m, 1e6);
    const DofMap d = build_dof_map(bnd, m);

    Solution s;
    s.conditioning = Conditioning::ScaledPhi;  // the one where the jw matters
    s.omega = 2.0 * M_PI * 1e6;
    s.x.resize(static_cast<std::size_t>(d.num_total));
    for (int i = 0; i < d.num_total; ++i) {
        s.x[static_cast<std::size_t>(i)] = Complex(std::cos(0.9 * i), std::sin(0.6 * i + 1.1));
    }

    const NodalPotential p = potential_at_nodes(m, bnd, d, s);
    const FieldOutput f = compute_fields(m, bnd, d, s);
    check(f.phi_node == p.value,
          "compute_fields reports exactly the Phi that potential_at_nodes does, bit for bit, "
          "including the ScaledPhi jw");
    check(f.num_nodes() == d.num_p2_nodes && f.num_vertices == m.num_nodes() &&
              f.num_edges == m.num_edges(),
          "and the field arrays cover every P2 node, indexed as Phi is");
    check(f.a_node.size() == f.phi_node.size() && f.e_node.size() == f.phi_node.size() &&
              f.b_node.size() == f.phi_node.size() && f.h_node.size() == f.phi_node.size(),
          "so a field and the potential can be read at the same node index");
}

// The mid-edge rule, stated as its own test because it is the one place the
// treatment of Phi and of the fields deliberately differs.
//
// A, B, H, E have no mid-edge degree of freedom, so a mid-edge value is the
// mean of the two endpoints. Phi does have one, and averaging it would replace
// the quadratic with its linear interpolant.
void test_mid_edge_fields_are_endpoint_means_but_phi_is_not() {
    Mesh m = make_cube();
    const BoundProblem bnd = bind_cube(m, 1e6);
    const DofMap d = all_free_dofs(m);

    // A Phi that is quadratic, so its mid-edge values are genuinely NOT the
    // endpoint means -- otherwise this test could not tell the two rules apart.
    Solution s;
    s.conditioning = Conditioning::Natural;
    s.omega = 2.0 * M_PI * 1e6;
    s.x.assign(static_cast<std::size_t>(d.num_total), Complex(0.0, 0.0));
    for (int i = 0; i < d.num_p2_nodes; ++i) {
        const Vec3 p = p2_node_position(m, i);
        s.x[static_cast<std::size_t>(d.num_a + i)] = Complex(p.x * p.x + 2.0 * p.y * p.z, 0.0);
    }
    for (int e = 0; e < m.num_edges(); ++e) {
        s.x[static_cast<std::size_t>(e)] = Complex(0.3 * std::sin(1.7 * e), 0.1 * std::cos(e));
    }

    const FieldOutput f = compute_fields(m, bnd, d, s);

    double worst_field = 0.0, worst_phi = 0.0;
    for (int e = 0; e < m.num_edges(); ++e) {
        const std::pair<int, int>& ends = m.edges[static_cast<std::size_t>(e)];
        const std::size_t a0 = static_cast<std::size_t>(ends.first);
        const std::size_t a1 = static_cast<std::size_t>(ends.second);
        const std::size_t mid = static_cast<std::size_t>(d.edge_p2(e));
        for (int i = 0; i < 3; ++i) {
            const std::size_t ui = static_cast<std::size_t>(i);
            worst_field = std::max(
                worst_field, std::abs(f.e_node[mid][ui] -
                                      0.5 * (f.e_node[a0][ui] + f.e_node[a1][ui])));
            worst_field = std::max(
                worst_field, std::abs(f.b_node[mid][ui] -
                                      0.5 * (f.b_node[a0][ui] + f.b_node[a1][ui])));
        }
        worst_phi = std::max(worst_phi, std::abs(f.phi_node[mid] -
                                                 0.5 * (f.phi_node[a0] + f.phi_node[a1])));
    }
    // Mid-edge FIELDS are evaluated there from the incident tets, not averaged
    // from the endpoints -- the endpoint mean cannot honour the material-side
    // selection at an interface. On a single-material mesh the two agree
    // closely, which is what this checks: near, but not identical.
    check(worst_field > 0.0,
          "a mid-edge field is evaluated there, not averaged from its endpoints, so it does "
          "not match the endpoint mean exactly on a rough field");

    // Against a field whose answer is known, the direct evaluation must be
    // exact. A linear Phi with A = 0 gives E = -grad(Phi), a CONSTANT, and a
    // constant is reproduced by either reconstruction -- so this pins the
    // direct evaluation rather than merely bounding it.
    Solution uniform;
    uniform.conditioning = Conditioning::Natural;
    uniform.omega = s.omega;
    uniform.x.assign(static_cast<std::size_t>(d.num_total), Complex(0.0, 0.0));
    const Vec3 grad{1.5, -0.5, 2.25};
    for (int i = 0; i < d.num_p2_nodes; ++i) {
        const Vec3 p = p2_node_position(m, i);
        uniform.x[static_cast<std::size_t>(d.num_a + i)] =
            Complex(grad.x * p.x + grad.y * p.y + grad.z * p.z, 0.0);
    }
    const FieldOutput uf = compute_fields(m, bnd, d, uniform);
    double worst_uniform = 0.0;
    for (int e = 0; e < m.num_edges(); ++e) {
        const std::size_t mid = static_cast<std::size_t>(d.edge_p2(e));
        worst_uniform = std::max(worst_uniform, std::abs(uf.e_node[mid][0] + grad.x));
        worst_uniform = std::max(worst_uniform, std::abs(uf.e_node[mid][2] + grad.z));
    }
    check(worst_uniform < 1e-12,
          "and on a field whose answer is known -- linear Phi, A = 0, so E is constant -- the "
          "mid-edge value is exact");

    // Phi, by contrast, is never averaged at all: its mid-edge value is a
    // solved P2 unknown, and the endpoint mean would discard the term that
    // makes the space second order.
    check(worst_phi > 1e-6,
          "while the mid-edge Phi is NOT that mean -- it is the solved quadratic value "
          "(largest difference " + std::to_string(worst_phi) + ")");
}


// The material-interface flag, which is what keeps a meaningless nodal average
// from being read as a field.
//
// E's normal component genuinely jumps across a conductor/insulator interface,
// so averaging over tets on both sides gives a value that is neither. Measured
// on the 50 Hz cylinder: thresholding this flag to 0 takes the worst transverse
// E in the wire from 2705 V/m down to 9.3e-07.
void test_material_interface_flag() {
    Mesh m = make_cube();
    BoundProblem bnd = bind_cube(m, 1e6);
    const DofMap d = all_free_dofs(m);

    Solution s;
    s.conditioning = Conditioning::Natural;
    s.omega = 2.0 * M_PI * 1e6;
    s.x.assign(static_cast<std::size_t>(d.num_total), Complex(0.3, -0.2));

    // One body: nothing can straddle anything.
    const FieldOutput one = compute_fields(m, bnd, d, s);
    bool all_clear = true;
    for (unsigned char f : one.on_material_interface) {
        if (f != 0) all_clear = false;
    }
    check(all_clear && one.num_interface_nodes == 0,
          "with a single body no node is on a material interface");

    // Split the mesh in two by fiat. The geometry is untouched, so any node the
    // flag picks up is one whose incident tets really do span both.
    BoundBody second = bnd.bodies[0];
    second.name = "B2";
    second.sigma = 0.0;
    bnd.bodies.push_back(second);
    const int half = m.num_tets() / 2;
    for (int t = half; t < m.num_tets(); ++t) {
        bnd.body_of_tet[static_cast<std::size_t>(t)] = 1;
    }

    const FieldOutput two = compute_fields(m, bnd, d, s);
    check(two.num_interface_nodes > 0,
          "splitting the mesh into two bodies marks " +
              std::to_string(two.num_interface_nodes) + " nodes");

    // Recomputed independently. A node -- vertex OR mid-edge -- is on an
    // interface exactly when sigma or mu is not the same across all the tets
    // that touch IT. Note this is not the same as propagating a vertex's flag
    // along its edges: an edge running from a conductor-surface vertex out into
    // the air touches only air tets, so it is uniform and NOT flagged, even
    // though one of its endpoints is.
    std::vector<int> seen(static_cast<std::size_t>(d.num_p2_nodes), -1);
    std::vector<unsigned char> want(static_cast<std::size_t>(d.num_p2_nodes), 0);
    for (int t = 0; t < m.num_tets(); ++t) {
        const std::size_t ut = static_cast<std::size_t>(t);
        const int b = bnd.body_of_tet[ut];
        for (int k = 0; k < 10; ++k) {
            const std::size_t u = static_cast<std::size_t>(
                k < 4 ? m.tets[ut][static_cast<std::size_t>(k)]
                      : d.edge_p2(m.tet_edges[ut][static_cast<std::size_t>(k - 4)]));
            if (seen[u] < 0) {
                seen[u] = b;
            } else if (seen[u] != b) {
                want[u] = 1;
            }
        }
    }
    check(two.on_material_interface == want,
          "and the flag is exactly that set, vertices and mid-edge nodes alike");

    // sigma is carried per tet, which is the only place it is single-valued.
    bool sigma_ok = true;
    for (int t = 0; t < m.num_tets(); ++t) {
        const double want_s =
            bnd.bodies[static_cast<std::size_t>(bnd.body_of_tet[static_cast<std::size_t>(t)])].sigma;
        if (two.sigma_tet[static_cast<std::size_t>(t)] != want_s) sigma_ok = false;
    }
    check(sigma_ok, "and sigma_tet is each tet's own body's sigma, so J = sigma E is formed "
                    "where sigma is single-valued");
}


// fields.out: every node, every field, as text.
//
// The risk in a 31-column file is that the header and the rows drift apart --
// a column added to one and not the other produces a file that parses and is
// silently mislabelled from that column on. So this counts them against each
// other rather than trusting either.
void test_write_solution() {
    Mesh m = make_cube();
    const BoundProblem bnd = bind_cube(m, 1e6);
    const DofMap d = build_dof_map(bnd, m);

    Solution s;
    s.conditioning = Conditioning::RowScaled;
    s.omega = 2.0 * M_PI * 1e6;
    s.x.resize(static_cast<std::size_t>(d.num_total));
    for (int i = 0; i < d.num_total; ++i) {
        s.x[static_cast<std::size_t>(i)] = Complex(std::sin(0.5 * i), std::cos(0.3 * i + 0.2));
    }

    const NodalPotential p = potential_at_nodes(m, bnd, d, s);
    const FieldOutput f = compute_fields(m, bnd, d, s);

    const std::string path = temp_path("fields.out");
    const WriteStats stats = write_solution(path, p, f, s);
    check(stats.nodes == p.size(), "one entry per P2 node");
    check(stats.bytes > 0, "the file is written");

    std::ifstream in(path);
    std::string line, header;
    int rows = 0, wrong_width = 0;
    std::vector<std::string> first_row;
    while (std::getline(in, line)) {
        if (!line.empty() && line[0] == '#') {
            if (line.rfind("# index ", 0) == 0) header = line;
            continue;
        }
        if (line.empty()) continue;
        std::istringstream row(line);
        std::vector<std::string> cols;
        std::string tok;
        while (row >> tok) cols.push_back(tok);
        if (rows == 0) first_row = cols;
        if (cols.size() != first_row.size()) ++wrong_width;
        ++rows;
    }
    check(rows == p.size(), "a row per node, and no more");
    check(wrong_width == 0, "every row has the same number of columns");

    // The header names its columns. Count them and require the same number.
    std::istringstream hs(header.substr(1));
    int header_cols = 0;
    std::string h;
    while (hs >> h) ++header_cols;
    check(header_cols == static_cast<int>(first_row.size()),
          "the header names exactly as many columns as the rows carry (" +
              std::to_string(header_cols) + " vs " + std::to_string(first_row.size()) +
              ") -- otherwise every column past the mismatch is silently mislabelled");
    check(first_row.size() == 4u + 2u + 30u + 1u,
          "index, position, Phi, five complex vectors (A B H E J), and the flag");

    // The values are the ones compute_fields produced, not a re-derivation.
    // Column 4 is Re(Phi); the last is the interface flag.
    std::ifstream again(path);
    int checked = 0;
    bool values_match = true;
    while (std::getline(again, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream row(line);
        int idx = 0;
        double x = 0, y = 0, z = 0, pr = 0, pi = 0;
        row >> idx >> x >> y >> z >> pr >> pi;
        const std::size_t u = static_cast<std::size_t>(idx);
        if (pr != f.phi_node[u].real() || pi != f.phi_node[u].imag()) values_match = false;
        if (x != p.position[u].x || z != p.position[u].z) values_match = false;
        ++checked;
    }
    check(values_match && checked == p.size(),
          "positions and Phi read back bit-for-bit, so to_chars round-trips them");

    std::remove(path.c_str());
}


// One file per field: each must carry its own and NOT the others.
//
// The failure mode a selector invites is a block guarded by the wrong flag, so
// B_field.vtk quietly ships E as well, or ships nothing. Checking only that the
// wanted array is present would miss both halves of that, so this checks
// presence AND absence, for every field against every other.
void test_field_set_selection() {
    Mesh m = make_cube();
    const BoundProblem bnd = bind_cube(m, 1e6);
    const DofMap d = build_dof_map(bnd, m);

    Solution s;
    s.conditioning = Conditioning::RowScaled;
    s.omega = 2.0 * M_PI * 1e6;
    s.x.resize(static_cast<std::size_t>(d.num_total));
    for (int i = 0; i < d.num_total; ++i) {
        s.x[static_cast<std::size_t>(i)] = Complex(std::sin(0.4 * i), std::cos(0.7 * i));
    }
    const NodalPotential p = potential_at_nodes(m, bnd, d, s);
    const FieldOutput f = compute_fields(m, bnd, d, s);

    struct Case {
        FieldSet set;
        const char* name;
        const char* array;   // the nodal VECTORS this file must have
        int text_columns;    // index + x,y,z + 6 + iface
    };
    const Case cases[] = {
        {FieldSet::A, "A", "A_real", 11},
        {FieldSet::B, "B", "B_real", 11},
        {FieldSet::H, "H", "H_real", 11},
        {FieldSet::E, "E", "E_real", 11},
        {FieldSet::J, "J", "J_real", 11},
    };
    const char* every[] = {"A_real", "B_real", "H_real", "E_real", "J_real"};

    for (const Case& c : cases) {
        check(std::string(field_set_name(c.set)) == c.name,
              std::string("field_set_name gives \"") + c.name + "\"");

        // --- the VTK ---
        const std::string vpath = temp_path((std::string(c.name) + "_field.vtk").c_str());
        const WriteStats vs = write_vtk(vpath, m, p, s, &f, c.set);
        check(vs.bytes > 0, std::string(c.name) + "_field.vtk is written");

        std::ifstream vin(vpath);
        std::string line;
        std::vector<std::string> arrays;
        while (std::getline(vin, line)) {
            std::istringstream row(line);
            std::string word, name;
            row >> word >> name;
            if (word == "VECTORS" || word == "SCALARS") arrays.push_back(name);
        }
        vin.close();

        for (const char* a : every) {
            const bool present =
                std::find(arrays.begin(), arrays.end(), std::string(a)) != arrays.end();
            const bool wanted = std::string(a) == c.array;
            check(present == wanted,
                  std::string(c.name) + "_field.vtk " + (wanted ? "has " : "does NOT have ") + a);
        }
        // Phi belongs to the potential file, not to a field file.
        check(std::find(arrays.begin(), arrays.end(), std::string("phi_real")) == arrays.end(),
              std::string(c.name) + "_field.vtk leaves Phi out");
        // But the interface flag travels with every field, because without it
        // the nodal values cannot be read safely.
        check(std::find(arrays.begin(), arrays.end(), std::string("material_interface")) !=
                  arrays.end(),
              std::string(c.name) + "_field.vtk still carries material_interface");
        std::remove(vpath.c_str());

        // --- the text file ---
        const std::string tpath = temp_path((std::string(c.name) + "_field.out").c_str());
        const WriteStats ts = write_solution(tpath, p, f, s, c.set);
        check(ts.nodes == p.size(), std::string(c.name) + "_field.out has a row per node");

        std::ifstream tin(tpath);
        std::string header;
        int rows = 0, header_cols = 0, row_cols = 0;
        while (std::getline(tin, line)) {
            if (!line.empty() && line[0] == '#') {
                if (line.rfind("# index ", 0) == 0) header = line;
                continue;
            }
            if (line.empty()) continue;
            if (rows == 0) {
                std::istringstream row(line);
                std::string tok;
                while (row >> tok) ++row_cols;
            }
            ++rows;
        }
        tin.close();
        std::istringstream hs(header.substr(1));
        std::string h;
        while (hs >> h) ++header_cols;
        check(row_cols == c.text_columns,
              std::string(c.name) + "_field.out carries " + std::to_string(c.text_columns) +
                  " columns, not the combined 31");
        check(header_cols == row_cols,
              std::string(c.name) + "_field.out header names exactly its columns");
        std::remove(tpath.c_str());
    }

    // And FieldSet::All still writes everything, so the combined file did not
    // become a casualty of the split.
    const std::string apath = temp_path("all.vtk");
    write_vtk(apath, m, p, s, &f, FieldSet::All);
    std::ifstream ain(apath);
    std::string line;
    int found = 0;
    while (std::getline(ain, line)) {
        std::istringstream row(line);
        std::string word, name;
        row >> word >> name;
        for (const char* a : every) {
            if (name == a) ++found;
        }
        if (name == "phi_real") ++found;
    }
    ain.close();
    check(found == 6, "FieldSet::All still writes Phi and all five fields");
    std::remove(apath.c_str());
}


// J = sigma E, and what happens at a conductor surface.
//
// J is the one field where averaging across a material interface has a RIGHT
// answer rather than no answer. At a conductor surface J is finite inside and
// zero outside, and the meaningful value -- the one a skin effect makes
// largest -- is the inside limit. So J averages over conducting tets only.
//
// An earlier version zeroed J at every interface node, on the same reasoning
// used for E. On the cylinder that put zero current density on the entire
// conductor surface, where the skin effect makes it maximal: a nodal J plot
// showed a hollow shell, exactly backwards.
void test_current_density() {
    Mesh m = make_cube();
    BoundProblem bnd = bind_cube(m, 1e6);
    const DofMap d = all_free_dofs(m);

    Solution s;
    s.conditioning = Conditioning::Natural;
    s.omega = 2.0 * M_PI * 1e6;
    s.x.resize(static_cast<std::size_t>(d.num_total));
    for (int i = 0; i < d.num_total; ++i) {
        s.x[static_cast<std::size_t>(i)] = Complex(std::cos(0.6 * i), std::sin(0.8 * i));
    }

    // --- one body: every tet conducts, so J's average runs over the same tets
    //     E's does, and at a VERTEX the two must agree to round-off.
    const FieldOutput one = compute_fields(m, bnd, d, s);
    const double sigma = bnd.bodies[0].sigma;
    check(sigma > 0.0, "the cube's body conducts, so this test can see anything at all");

    double worst_vertex = 0.0, scale = 0.0, worst_tet = 0.0;
    for (int v = 0; v < m.num_nodes(); ++v) {
        const std::size_t u = static_cast<std::size_t>(v);
        for (int k = 0; k < 3; ++k) {
            const std::size_t uk = static_cast<std::size_t>(k);
            worst_vertex =
                std::max(worst_vertex, std::abs(one.j_node[u][uk] - sigma * one.e_node[u][uk]));
            scale = std::max(scale, std::abs(sigma * one.e_node[u][uk]));
        }
    }
    for (int t = 0; t < m.num_tets(); ++t) {
        const std::size_t ut = static_cast<std::size_t>(t);
        for (int k = 0; k < 3; ++k) {
            const std::size_t uk = static_cast<std::size_t>(k);
            worst_tet = std::max(worst_tet, std::abs(one.j_tet[ut][uk] - sigma * one.e_tet[ut][uk]));
        }
    }
    check(worst_vertex <= 1e-12 * scale,
          "with one body, J = sigma E at every vertex to round-off");
    check(worst_tet == 0.0, "and exactly in every tet, where no averaging happens at all");

    // Mid-edge nodes too: the relation is applied at the node, so it holds
    // wherever a node is, without exception.
    double worst_mid = 0.0;
    for (int e = 0; e < m.num_edges(); ++e) {
        const std::size_t u = static_cast<std::size_t>(d.edge_p2(e));
        for (int k = 0; k < 3; ++k) {
            const std::size_t uk = static_cast<std::size_t>(k);
            worst_mid = std::max(worst_mid, std::abs(one.j_node[u][uk] - sigma * one.e_node[u][uk]));
        }
    }
    check(worst_mid == 0.0, "and at mid-edge nodes exactly, since sigma is applied at the node");

    // B = mu H at every node, the same way.
    const double mu = kMu0 * bnd.bodies[0].mu_r;
    double worst_h = 0.0, h_scale = 0.0;
    for (int i = 0; i < d.num_p2_nodes; ++i) {
        const std::size_t u = static_cast<std::size_t>(i);
        for (int k = 0; k < 3; ++k) {
            const std::size_t uk = static_cast<std::size_t>(k);
            worst_h = std::max(worst_h, std::abs(one.b_node[u][uk] - mu * one.h_node[u][uk]));
            h_scale = std::max(h_scale, std::abs(one.b_node[u][uk]));
        }
    }
    check(worst_h <= 1e-12 * h_scale, "and B = mu H at every node, to round-off");

    // --- two bodies, one an insulator. This is the case that matters.
    BoundBody insulator = bnd.bodies[0];
    insulator.name = "B2";
    insulator.sigma = 0.0;
    bnd.bodies.push_back(insulator);
    const int half = m.num_tets() / 2;
    for (int t = half; t < m.num_tets(); ++t) bnd.body_of_tet[static_cast<std::size_t>(t)] = 1;

    const FieldOutput two = compute_fields(m, bnd, d, s);
    check(two.num_interface_nodes > 0, "the split creates interface nodes");

    // The point of the change: an interface node in contact with the conductor
    // carries the CONDUCTOR-SIDE current density, not zero.
    int iface_live = 0, iface_dead = 0;
    for (int i = 0; i < d.num_p2_nodes; ++i) {
        const std::size_t u = static_cast<std::size_t>(i);
        if (two.on_material_interface[u] == 0) continue;
        double mag = 0.0;
        for (int k = 0; k < 3; ++k) mag += std::abs(two.j_node[u][static_cast<std::size_t>(k)]);
        if (mag > 0.0) {
            ++iface_live;
        } else {
            ++iface_dead;
        }
    }
    check(iface_live > 0,
          "an interface node touching the conductor carries the conductor-side J (" +
              std::to_string(iface_live) + " of " + std::to_string(iface_live + iface_dead) +
              ") -- zeroing these put a hollow shell where the skin effect is strongest");

    // In the insulator's interior, J is zero because sigma is -- a different
    // reason from the interface one, and both have to hold.
    int insulating_tets = 0;
    bool insulator_is_zero = true;
    for (int t = 0; t < m.num_tets(); ++t) {
        const std::size_t ut = static_cast<std::size_t>(t);
        if (two.sigma_tet[ut] != 0.0) continue;
        ++insulating_tets;
        for (int k = 0; k < 3; ++k) {
            if (two.j_tet[ut][static_cast<std::size_t>(k)] != Complex(0.0, 0.0)) {
                insulator_is_zero = false;
            }
        }
    }
    check(insulating_tets > 0 && insulator_is_zero,
          "J is zero in all " + std::to_string(insulating_tets) + " insulating tets");

    // And E there is NOT zero, so that zero is sigma's doing rather than a
    // collapsed solution.
    double e_in_insulator = 0.0;
    for (int t = 0; t < m.num_tets(); ++t) {
        const std::size_t ut = static_cast<std::size_t>(t);
        if (two.sigma_tet[ut] != 0.0) continue;
        for (int k = 0; k < 3; ++k) {
            e_in_insulator =
                std::max(e_in_insulator, std::abs(two.e_tet[ut][static_cast<std::size_t>(k)]));
        }
    }
    check(e_in_insulator > 0.0,
          "while E in the insulator is not zero, so J vanishing there is sigma's doing");

    // A node reached by no conducting tet at all must be zero: that is the
    // insulator's interior, and it is an answer rather than a gap.
    int deep = 0;
    bool deep_is_zero = true;
    for (int i = 0; i < d.num_p2_nodes; ++i) {
        const std::size_t u = static_cast<std::size_t>(i);
        if (two.on_material_interface[u] != 0) continue;
        bool touches_conductor = false;
        for (int t = 0; t < m.num_tets(); ++t) {
            if (two.sigma_tet[static_cast<std::size_t>(t)] <= 0.0) continue;
            for (int c = 0; c < 4; ++c) {
                if (static_cast<std::size_t>(
                        m.tets[static_cast<std::size_t>(t)][static_cast<std::size_t>(c)]) == u) {
                    touches_conductor = true;
                }
            }
        }
        if (touches_conductor || i >= m.num_nodes()) continue;
        ++deep;
        for (int k = 0; k < 3; ++k) {
            if (two.j_node[u][static_cast<std::size_t>(k)] != Complex(0.0, 0.0)) {
                deep_is_zero = false;
            }
        }
    }
    check(deep_is_zero, "a vertex no conducting tet reaches carries exactly zero J (" +
                            std::to_string(deep) + " of them)");
}

}  // namespace

// ---------------------------------------------------------------------------
// write_postprocess_manifest
//
// The manifest is the entire interface between the solver and the renderer, so
// the two things it must get right are the ones no caller can check for itself:
// the resolved body tag, and METRES. A millimetre offset left unconverted would
// slice a metre mesh 1000x too far out and produce an empty picture that looks
// like a modelling mistake rather than a unit bug.

std::string read_all(const std::string& path) {
    std::FILE* f = std::fopen(path.c_str(), "rb");
    if (f == nullptr) return {};
    std::string out;
    char buf[4096];
    std::size_t n = 0;
    while ((n = std::fread(buf, 1, sizeof buf, f)) > 0) out.append(buf, n);
    std::fclose(f);
    return out;
}

void test_postprocess_manifest() {
    using aphi_solver::PostComponent;
    using aphi_solver::PostDisplay;
    using aphi_solver::PostField;
    using aphi_solver::PostGeometry;
    using aphi_solver::PostPlane;
    using aphi_solver::PostprocessRequest;

    aphi_solver::Problem p;
    p.length_unit = aphi_solver::LengthUnit::Millimetre;

    PostprocessRequest a;
    a.name = "PP1";
    a.geometry = PostGeometry::Body;
    a.body = "ring";
    a.field = PostField::J;
    a.display = PostDisplay::Vector;
    a.phase_deg = 90.0;
    p.postprocess.push_back(a);

    PostprocessRequest b;
    b.name = "PP2";
    b.geometry = PostGeometry::Plane;
    b.plane = PostPlane::YZ;
    b.offset = 2.5;  // millimetres
    b.field = PostField::B;
    p.postprocess.push_back(b);

    PostprocessRequest c;
    c.name = "PP3";
    c.geometry = PostGeometry::Points;
    c.points.push_back(aphi_solver::Vec3(6.5, 0.0, 0.0));  // millimetres
    c.field = PostField::Phi;
    c.display = PostDisplay::Phase;
    p.postprocess.push_back(c);

    aphi_solver::BoundProblem bound;
    aphi_solver::BoundBody bb;
    bb.name = "B1";
    bb.volume = "ring";
    bb.tag = 7;
    bound.bodies.push_back(bb);

    const std::string path = "test_manifest.json";
    const aphi_solver::WriteStats st =
        aphi_solver::write_postprocess_manifest(path, p, bound, 50.0, "");
    check(st.nodes == 3, "manifest reports three requests");
    check(st.bytes > 0, "manifest wrote something");

    const std::string text = read_all(path);

    // The resolved tag, which only binding knows.
    check(text.find("\"body_tag\": 7") != std::string::npos,
          "the Physical Volume name resolved to its tag");
    check(text.find("\"body\": \"ring\"") != std::string::npos, "and the name is carried too");

    // METRES. 2.5 mm is 0.0025 m, and 6.5 mm is 0.0065 m.
    check(text.find("0.0025") != std::string::npos,
          "a plane offset is converted from the file's length_unit to metres");
    check(text.find("2.5,") == std::string::npos && text.find(" 2.5\n") == std::string::npos,
          "and the unconverted millimetre value is NOT present");
    check(text.find("0.0065") != std::string::npos, "probe points converted to metres as well");

    // Which file, and what the arrays there are called.
    check(text.find("\"vtk\": \"J_field.vtk\"") != std::string::npos, "J lives in J_field.vtk");
    check(text.find("\"cell_prefix\": \"J_cell\"") != std::string::npos, "with per-cell arrays");
    check(text.find("\"vtk\": \"potential.vtk\"") != std::string::npos,
          "Phi lives in potential.vtk");
    check(text.find("\"point_prefix\": \"phi\"") != std::string::npos, "under phi_*");
    check(text.find("\"cell_prefix\": \"\"") != std::string::npos,
          "and has no per-cell companion, being a scalar");
    check(text.find("\"is_vector\": false") != std::string::npos, "phi is flagged as a scalar");

    // phase_deg appears only where the display reads it.
    const std::size_t pp1 = text.find("\"name\": \"PP1\"");
    const std::size_t pp2 = text.find("\"name\": \"PP2\"");
    const std::size_t ph = text.find("\"phase_deg\"");
    check(ph != std::string::npos && ph > pp1 && ph < pp2,
          "phase_deg is written for the vector display and only there");
    check(text.find("\"phase_deg\"", pp2) == std::string::npos,
          "and not for the phase-independent one");

    std::remove(path.c_str());

    // A file that asked for no pictures must not grow an empty manifest.
    aphi_solver::Problem none;
    const aphi_solver::WriteStats empty =
        aphi_solver::write_postprocess_manifest("should_not_exist.json", none, bound, 50.0, "");
    check(empty.bytes == 0, "no requests writes nothing");
    check(read_all("should_not_exist.json").empty(), "and leaves no file behind");

    // A sweep's suffix reaches the filenames, so two solves do not read each
    // other's fields.
    const std::string sp = "test_manifest_2.json";
    aphi_solver::write_postprocess_manifest(sp, p, bound, 1000.0, "_f2");
    const std::string swept = read_all(sp);
    check(swept.find("\"vtk\": \"J_field_f2.vtk\"") != std::string::npos,
          "the suffix reaches the field filenames");
    check(swept.find("\"frequency_hz\": 1000") != std::string::npos, "and the frequency is its own");
    std::remove(sp.c_str());
}

int main() {
    test_node_layout_and_positions();
    test_every_status_is_honoured();
    test_scaled_phi_is_multiplied_by_j_omega();
    test_round_trip_is_lossless();
    test_buffered_writing_is_faster();
    test_buffer_handles_awkward_sizes();
    test_vtk_quadratic_node_order();
    test_edge_sign_convention_matches_assembly();
    test_uniform_b_is_reproduced_exactly();
    test_linear_phi_gives_exact_uniform_e();
    test_vertex_averaging_is_volume_weighted();
    test_phi_matches_potential_at_nodes();
    test_mid_edge_fields_are_endpoint_means_but_phi_is_not();
    test_material_interface_flag();
    test_write_solution();
    test_field_set_selection();
    test_current_density();
    test_postprocess_manifest();

    std::cout << (g_checks - g_failures) << "/" << g_checks << " checks passed\n";
    return g_failures == 0 ? 0 : 1;
}
