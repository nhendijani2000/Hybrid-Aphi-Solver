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

}  // namespace

int main() {
    test_node_layout_and_positions();
    test_every_status_is_honoured();
    test_scaled_phi_is_multiplied_by_j_omega();
    test_round_trip_is_lossless();
    test_buffered_writing_is_faster();
    test_buffer_handles_awkward_sizes();

    std::cout << (g_checks - g_failures) << "/" << g_checks << " checks passed\n";
    return g_failures == 0 ? 0 : 1;
}
