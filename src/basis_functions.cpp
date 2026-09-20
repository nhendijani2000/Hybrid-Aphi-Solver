#include "aphi_solver/basis_functions.hpp"

#include <cmath>
#include <stdexcept>
#include <utility>

namespace aphi_solver {

namespace {

/// Inverts a 4x4 matrix via Gauss-Jordan elimination with partial pivoting.
/// Throws std::runtime_error if the matrix is (numerically) singular -- for
/// the M built from a tet's vertex coordinates in compute_tet_geometry below,
/// that means a degenerate (zero-volume) tet.
std::array<std::array<double, 4>, 4> invert4x4(std::array<std::array<double, 4>, 4> m) {
    std::array<std::array<double, 4>, 4> inv{};
    for (int i = 0; i < 4; ++i) {
        inv[i][i] = 1.0;
    }

    for (int col = 0; col < 4; ++col) {
        int pivot_row = col;
        double pivot_val = std::abs(m[col][col]);
        for (int row = col + 1; row < 4; ++row) {
            if (std::abs(m[row][col]) > pivot_val) {
                pivot_row = row;
                pivot_val = std::abs(m[row][col]);
            }
        }
        if (pivot_val < 1e-14) {
            throw std::runtime_error("compute_tet_geometry: degenerate (zero-volume) tet");
        }
        if (pivot_row != col) {
            std::swap(m[pivot_row], m[col]);
            std::swap(inv[pivot_row], inv[col]);
        }

        double pivot = m[col][col];
        for (int k = 0; k < 4; ++k) {
            m[col][k] /= pivot;
            inv[col][k] /= pivot;
        }

        for (int row = 0; row < 4; ++row) {
            if (row == col) continue;
            double factor = m[row][col];
            if (factor == 0.0) continue;
            for (int k = 0; k < 4; ++k) {
                m[row][k] -= factor * m[col][k];
                inv[row][k] -= factor * inv[col][k];
            }
        }
    }
    return inv;
}

}  // namespace

TetGeometry compute_tet_geometry(const Mesh& mesh, int t) {
    const TetVerts& verts = mesh.tets[t];

    // L_i(x,y,z) = a[i] + grad_L[i].dot((x,y,z)) must satisfy L_i(p_j) =
    // delta_ij for the tet's 4 vertices p_0..p_3. Writing M's row j as
    // [1, x_j, y_j, z_j], that is M * v_i = e_i for each i, i.e. v_i (the
    // stacked [a_i, grad_L_i] coefficients) is column i of M^{-1}. This is
    // generic and correct regardless of vertex-ordering orientation.
    std::array<std::array<double, 4>, 4> M{};
    for (int j = 0; j < 4; ++j) {
        const Vec3& p = mesh.nodes[verts[j]];
        M[j][0] = 1.0;
        M[j][1] = p.x;
        M[j][2] = p.y;
        M[j][3] = p.z;
    }

    std::array<std::array<double, 4>, 4> Minv = invert4x4(M);

    TetGeometry g;
    for (int i = 0; i < 4; ++i) {
        g.a[i] = Minv[0][i];
        g.grad_L[i] = Vec3(Minv[1][i], Minv[2][i], Minv[3][i]);
    }
    g.volume = std::abs(mesh.signed_tet_volume(t));
    return g;
}

std::array<double, 4> evaluate_barycentric(const TetGeometry& g, const Vec3& p) {
    std::array<double, 4> L{};
    for (int i = 0; i < 4; ++i) {
        L[i] = g.a[i] + g.grad_L[i].dot(p);
    }
    return L;
}

Vec3 whitney_edge_value(const TetGeometry& g, int local_edge, const std::array<double, 4>& L) {
    const auto& [vi, vj] = kTetLocalEdgeVerts[local_edge];
    // N_ij(p) = L_i(p) * grad(L_j) - L_j(p) * grad(L_i)
    return g.grad_L[vj] * L[vi] - g.grad_L[vi] * L[vj];
}

Vec3 whitney_edge_curl(const TetGeometry& g, int local_edge) {
    const auto& [vi, vj] = kTetLocalEdgeVerts[local_edge];
    // curl(N_ij) = 2 * grad(L_i) x grad(L_j) -- constant over the tet.
    return g.grad_L[vi].cross(g.grad_L[vj]) * 2.0;
}

double p2_nodal_value(int local_node, const std::array<double, 4>& L) {
    if (local_node < 4) {
        double li = L[local_node];
        return (2.0 * li - 1.0) * li;
    }
    const auto& [v0, v1] = kTetLocalEdgeVerts[local_node - 4];
    return 4.0 * L[v0] * L[v1];
}

Vec3 p2_nodal_gradient(const TetGeometry& g, int local_node, const std::array<double, 4>& L) {
    if (local_node < 4) {
        double li = L[local_node];
        // N_i = (2*Li - 1) * Li = 2*Li^2 - Li, so d/dLi = 4*Li - 1, and
        // grad(N_i) = (4*Li - 1) * grad(Li) by the chain rule.
        return g.grad_L[local_node] * (4.0 * li - 1.0);
    }
    const auto& [v0, v1] = kTetLocalEdgeVerts[local_node - 4];
    // N_{4+e} = 4 * L_v0 * L_v1, so grad = 4*(L_v1*grad(L_v0) + L_v0*grad(L_v1)).
    return g.grad_L[v0] * (4.0 * L[v1]) + g.grad_L[v1] * (4.0 * L[v0]);
}

}  // namespace aphi_solver
