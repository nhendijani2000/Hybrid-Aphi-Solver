// tools/compare_gauges.cpp -- a small CLI for comparing the two tree-cotree
// gauge variants (docs/TREE_COTREE_GAUGE.md) on a real mesh, rather than
// only the tiny hand-built meshes in tests/test_gauge_variants.cpp.
//
// Usage:
//   compare_gauges <mesh.msh> [--gauge=A|D|both] [--pec-nodes=i,j,k,...] [--skip-kappa]
//
// M = C^T * C (vacuum, nu = 1) is used as the test matrix, exactly as
// tests/test_gauge_variants.cpp already does -- there is no real assembly
// pipeline yet (Phase 04), so this is the honest, physically legitimate
// stand-in available today, not a placeholder invented just for this tool.
// See docs/CONDITIONING.md and docs/LINEAR_SOLVER.md for how this tool's
// output feeds into the (gauge, frequency_scaling) input-file design.

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include "aphi_solver/gauge_variants.hpp"
#include "aphi_solver/gmsh_reader.hpp"
#include "aphi_solver/incidence.hpp"
#include "aphi_solver/mesh.hpp"
#include "aphi_solver/tree_cotree.hpp"

using namespace aphi_solver;

namespace {

void print_usage(const char* argv0) {
    std::cerr << "Usage: " << argv0 << " <mesh.msh> [--gauge=A|D|both] [--pec-nodes=i,j,k,...] [--skip-kappa]\n\n"
              << "  --gauge=A|D|both   Which gauge variant(s) to build and report (default: both)\n"
              << "  --pec-nodes=LIST   Comma-separated 0-based node indices to tag as PEC\n"
              << "                     (default: none -- a single reference node is picked\n"
              << "                     automatically, per build_tree_cotree's no-PEC fallback)\n"
              << "  --skip-kappa       Skip the condition-number estimate (CG-based power\n"
              << "                     iteration) and report only structural stats -- reduced\n"
              << "                     size, nonzeros, fill-in. Use this on large meshes where\n"
              << "                     you only need to check fill-in/size quickly and don't\n"
              << "                     need the (slower) conditioning comparison right now.\n"
              << "  --help             Show this message\n\n"
              << "Note: M = C^T*C (vacuum, nu=1) is used as the test matrix -- there is no real\n"
              << "assembly pipeline yet (Phase 04), so this is the same honest stand-in\n"
              << "tests/test_gauge_variants.cpp already uses, not a placeholder invented for\n"
              << "this tool. See tools/generate_cube_mesh.py for a way to make bigger test\n"
              << "meshes than the ones checked into meshes/.\n";
}

struct Args {
    std::string mesh_path;
    bool want_a = true;
    bool want_d = true;
    bool skip_kappa = false;
    std::vector<int> pec_nodes;
};

std::vector<int> parse_int_list(const std::string& s) {
    std::vector<int> result;
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, ',')) {
        if (!item.empty()) result.push_back(std::stoi(item));
    }
    return result;
}

Args parse_args(int argc, char** argv) {
    Args args;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--help" || a == "-h") {
            print_usage(argv[0]);
            std::exit(0);
        } else if (a.rfind("--gauge=", 0) == 0) {
            const std::string v = a.substr(8);
            if (v == "A" || v == "a") {
                args.want_a = true;
                args.want_d = false;
            } else if (v == "D" || v == "d") {
                args.want_a = false;
                args.want_d = true;
            } else if (v == "both") {
                args.want_a = true;
                args.want_d = true;
            } else {
                throw std::runtime_error("unrecognized --gauge value: " + v + " (expected A, D, or both)");
            }
        } else if (a.rfind("--pec-nodes=", 0) == 0) {
            args.pec_nodes = parse_int_list(a.substr(12));
        } else if (a == "--skip-kappa") {
            args.skip_kappa = true;
        } else if (a.rfind("--", 0) == 0) {
            throw std::runtime_error("unrecognized option: " + a);
        } else if (args.mesh_path.empty()) {
            args.mesh_path = a;
        } else {
            throw std::runtime_error("unexpected extra argument: " + a);
        }
    }
    if (args.mesh_path.empty()) {
        throw std::runtime_error("missing required <mesh.msh> argument");
    }
    return args;
}

double elapsed_ms(std::chrono::steady_clock::time_point t0, std::chrono::steady_clock::time_point t1) {
    return std::chrono::duration<double, std::milli>(t1 - t0).count();
}

// Returns the estimated condition number, or NaN when `skip_kappa` is set --
// in which case the (expensive, CG-based) estimate is never computed at all,
// not just hidden from the printout. See --skip-kappa in print_usage() for
// why: on large meshes the conditioning estimate is the slow part, and this
// lets you check structural stats (size/nnz/fill-in) without paying for it.
double report_variant(const GaugeVariant& variant, bool skip_kappa) {
    const int n = variant.reduced_matrix.rows;
    const long long nnz = static_cast<long long>(variant.reduced_matrix.entries.size());
    const double density = (n > 0) ? static_cast<double>(nnz) / (static_cast<double>(n) * static_cast<double>(n)) : 0.0;

    std::cout << "  " << variant.name << ":\n"
              << "    reduced size: " << n << " x " << n << "\n"
              << "    nonzeros: " << nnz << " (density " << (density * 100.0) << "%)\n";

    if (skip_kappa) {
        std::cout << "    condition number estimate: skipped (--skip-kappa)\n";
        return std::numeric_limits<double>::quiet_NaN();
    }

    const auto t0 = std::chrono::steady_clock::now();
    const double kappa = estimate_condition_number(variant.reduced_matrix);
    const auto t1 = std::chrono::steady_clock::now();
    std::cout << "    condition number estimate: " << kappa << " (" << elapsed_ms(t0, t1) << " ms)\n";
    return kappa;
}

}  // namespace

int main(int argc, char** argv) {
    Args args;
    try {
        args = parse_args(argc, argv);
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n\n";
        print_usage(argv[0]);
        return 1;
    }

    Mesh mesh;
    try {
        mesh = read_gmsh_msh(args.mesh_path);
    } catch (const GmshReadError& e) {
        std::cerr << "error reading mesh: " << e.what() << "\n";
        return 1;
    }

    std::vector<bool> is_pec(static_cast<std::size_t>(mesh.num_nodes()), false);
    for (int idx : args.pec_nodes) {
        if (idx < 0 || idx >= mesh.num_nodes()) {
            std::cerr << "error: --pec-nodes index " << idx << " is out of range (mesh has "
                      << mesh.num_nodes() << " nodes, 0-based)\n";
            return 1;
        }
        is_pec[static_cast<std::size_t>(idx)] = true;
    }

    std::cout << "Mesh: " << args.mesh_path << "\n"
              << "  nodes=" << mesh.num_nodes() << " edges=" << mesh.num_edges()
              << " tets=" << mesh.num_tets() << " pec_nodes=" << args.pec_nodes.size() << "\n";

    const TreeCotreeResult tc = build_tree_cotree(mesh, is_pec);
    std::cout << "  tree-cotree: num_groups=" << tc.num_groups
              << " num_reference_groups=" << tc.num_reference_groups
              << " tree_edge_count=" << tc.tree_edge_count << "\n\n";

    // M = CT * C (vacuum stand-in) -- see the file header comment for why
    // this, not a real assembled A-Phi matrix, is used.
    const SparseMatrix C = build_curl_matrix(mesh);
    const SparseMatrix M = multiply(transpose(C), C);

    GaugeVariant gauge_a, gauge_d;
    double kappa_a = 0.0, kappa_d = 0.0;
    if (args.want_a) {
        gauge_a = build_albanese_rubinacci_gauge(M, tc);
        kappa_a = report_variant(gauge_a, args.skip_kappa);
    }
    if (args.want_d) {
        const EssentialIncidenceMatrix F = compute_essential_incidence_matrix(mesh, tc);
        gauge_d = build_munteanu_unsymmetric_gauge(M, tc, F);
        kappa_d = report_variant(gauge_d, args.skip_kappa);
    }

    if (args.want_a && args.want_d) {
        std::cout << "\n";
        if (args.skip_kappa) {
            std::cout << "  kappa_D / kappa_A = skipped (--skip-kappa)\n";
        } else {
            std::cout << "  kappa_D / kappa_A = " << (kappa_d / kappa_a)
                      << (kappa_d < kappa_a ? "  (D better-conditioned)" : "  (A better-conditioned or equal)") << "\n";
        }
        const long long nnz_a = static_cast<long long>(gauge_a.reduced_matrix.entries.size());
        const long long nnz_d = static_cast<long long>(gauge_d.reduced_matrix.entries.size());
        std::cout << "  nnz_D / nnz_A = " << (static_cast<double>(nnz_d) / static_cast<double>(nnz_a))
                  << "  (Method D's fill-in relative to Method A's exact submatrix)\n";
    }

    return 0;
}
