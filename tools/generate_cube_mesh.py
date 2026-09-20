#!/usr/bin/env python3
"""Generates a structured NxNxN-cube tetrahedral mesh as a Gmsh 2.2 ASCII
file (gmsh_reader.cpp's supported format), for exercising
tools/compare_gauges on something bigger than this project's 1-2-tet unit
test meshes.

Each unit cube is split into 6 tets via the standard Kuhn triangulation --
all six share the cube's main diagonal, and the triangulation is consistent
across a shared face between adjacent cubes (both cubes agree on how that
face's diagonal is drawn), so the result is a proper conforming tet mesh, not
a set of overlapping/mismatched tets.

Usage: python3 generate_cube_mesh.py [n] [out.msh]
  n         divisions per axis (default 4) -- (n+1)^3 nodes, 6*n^3 tets.
  out.msh   output path (default cube_<n>.msh)
"""
import sys


def generate(n, out_path):
    node_id = {}
    nodes = []
    for k in range(n + 1):
        for j in range(n + 1):
            for i in range(n + 1):
                idx = len(nodes) + 1  # Gmsh node ids are 1-based
                node_id[(i, j, k)] = idx
                nodes.append((idx, i, j, k))

    # Kuhn triangulation of the unit cube with corners (0,0,0)-(1,1,1):
    # 6 tets, each containing the main diagonal (0,0,0)-(1,1,1).
    kuhn_tets = [
        (0, 0, 0), (1, 0, 0), (1, 1, 0), (1, 1, 1),
        (0, 0, 0), (1, 1, 0), (0, 1, 0), (1, 1, 1),
        (0, 0, 0), (0, 1, 0), (0, 1, 1), (1, 1, 1),
        (0, 0, 0), (0, 1, 1), (0, 0, 1), (1, 1, 1),
        (0, 0, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1),
        (0, 0, 0), (1, 0, 1), (1, 0, 0), (1, 1, 1),
    ]
    tet_offsets = [kuhn_tets[t * 4:(t + 1) * 4] for t in range(6)]

    tets = []
    for k in range(n):
        for j in range(n):
            for i in range(n):
                for offsets in tet_offsets:
                    verts = [node_id[(i + dx, j + dy, k + dz)] for (dx, dy, dz) in offsets]
                    tets.append(verts)

    with open(out_path, "w") as f:
        f.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")
        f.write("$Nodes\n")
        f.write(f"{len(nodes)}\n")
        for idx, i, j, k in nodes:
            f.write(f"{idx} {i} {j} {k}\n")
        f.write("$EndNodes\n")
        f.write("$Elements\n")
        f.write(f"{len(tets)}\n")
        for eid, verts in enumerate(tets, start=1):
            f.write(f"{eid} 4 2 0 1 {verts[0]} {verts[1]} {verts[2]} {verts[3]}\n")
        f.write("$EndElements\n")

    print(f"wrote {out_path}: {len(nodes)} nodes, {len(tets)} tets (n={n})")


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    out_path = sys.argv[2] if len(sys.argv) > 2 else f"cube_{n}.msh"
    generate(n, out_path)
