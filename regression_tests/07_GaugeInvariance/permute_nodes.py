"""Rewrite a gmsh 4.1 .msh with the nodes listed in a different ORDER.

    python permute_nodes.py in.msh out.msh [seed]

The mesh that comes out is geometrically IDENTICAL: same nodes, same
coordinates, same tags, same elements. Only the order in which the nodes appear
in the file changes.

THAT IS THE WHOLE TRICK, and it is worth saying why it produces a different
spanning tree, because nothing about it is obvious from the file:

  * `Accumulator::add_node` in src/gmsh_reader.cpp assigns each node its
    internal index from `mesh.nodes.size()` -- that is, IN FILE ORDER. The gmsh
    tag is only a lookup key; it does not set the index.
  * `Mesh::build` sorts edges by their canonical (sorted) vertex pair, so edge
    numbering is a function of the node numbering.
  * `build_tree_cotree` roots at the lowest-numbered node of Dirichlet
    component 0 and walks neighbours in `NodeAdjacency` order, which is edge
    order.

So permuting the file permutes the node indices, which permutes the edge
indices, which moves the root and the traversal order, which yields a different
spanning tree -- over a problem that is unchanged in every physical respect.

The permutation is applied to the (tag, coordinate) PAIRS within each entity
block, so every node keeps its own tag and its own coordinates. `$Elements`
refers to nodes by tag, so it needs no edit at all.

WHY NOT JUST RENUMBER THE TAGS: because the reader ignores tag values for
ordering. Renumbering tags in place would produce a file that looks shuffled and
yields the identical tree -- a test that silently proves nothing.
"""
import random
import sys


def permute(src, dst, seed=12345):
    with open(src, "r") as f:
        lines = f.read().split("\n")

    try:
        i0 = lines.index("$Nodes")
        i1 = lines.index("$EndNodes")
    except ValueError:
        raise SystemExit("%s: no $Nodes section -- is this a 4.1 .msh?" % src)

    rng = random.Random(seed)
    head = lines[i0 + 1].split()
    num_blocks, num_nodes = int(head[0]), int(head[1])

    out = [lines[i0 + 1]]
    p = i0 + 2
    moved = 0
    for _ in range(num_blocks):
        bh = lines[p].split()
        n_in_block = int(bh[3])
        out.append(lines[p])
        p += 1
        tags = lines[p:p + n_in_block]
        p += n_in_block
        coords = lines[p:p + n_in_block]
        p += n_in_block

        order = list(range(n_in_block))
        rng.shuffle(order)
        moved += sum(1 for k, o in enumerate(order) if k != o)
        out.extend(tags[o] for o in order)
        out.extend(coords[o] for o in order)

    if p != i1:
        raise SystemExit("%s: node section did not parse cleanly (%d vs %d)"
                         % (src, p, i1))
    with open(dst, "w", newline="\n") as f:
        f.write("\n".join(lines[:i0 + 1] + out + lines[i1:]))
    return num_nodes, moved


if __name__ == "__main__":
    n, m = permute(sys.argv[1], sys.argv[2],
                   int(sys.argv[3]) if len(sys.argv) > 3 else 12345)
    print("  %s -> %s: %d nodes, %d moved (%.1f %%)"
          % (sys.argv[1], sys.argv[2], n, m, 100.0 * m / n))
