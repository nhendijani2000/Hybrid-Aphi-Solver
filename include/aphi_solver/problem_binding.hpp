#pragma once

#include <complex>
#include <string>
#include <vector>

#include "aphi_solver/input_file.hpp"  // InputError
#include "aphi_solver/mesh.hpp"
#include "aphi_solver/problem.hpp"
#include "aphi_solver/tree_cotree.hpp"

namespace aphi_solver {

/// Binding: resolving a parsed `Problem` against a real mesh.
///
/// `Mesh` gives a complete edge/face/tet topology but **no bodies and no
/// surfaces** -- only `tet_tags` (one integer per tet) and
/// `tagged_boundary_faces` (a flat list of node triples in file order, not
/// grouped and not resolved to face indices). This stage turns those into
/// the containers the rest of the solver wants, and checks everything that
/// needs the mesh to be checkable at all.
///
/// The division of labour with `parse_input_file` is deliberate: a typo is a
/// different class of mistake from a name that does not exist in the mesh,
/// and reporting them at separate stages is what lets each message name the
/// right thing.

/// One body, resolved: its material plus the tets it owns.
struct BoundBody {
    std::string name;    ///< the `[Body ...]` handle
    std::string volume;  ///< as written in the input file
    int tag = -1;        ///< the Physical Volume tag it resolved to
    double sigma = 0.0;
    double eps_r = 1.0;
    double mu_r = 1.0;
    std::vector<int> tets;  ///< ascending

    bool is_conductor() const { return sigma > 0.0; }
};

/// One tagged surface, resolved to topology. `faces` holds global face
/// indices; `nodes` and `edges` are the vertices and edges those faces use,
/// each ascending and unique.
///
/// `nodes` is the P1 vertex set. Phi's P2 element also has an unknown at
/// every edge midpoint, and `edges` is exactly that set -- so a terminal's
/// full P2 node set is `nodes` plus one per entry of `edges`. Keeping them
/// separate rather than pre-merged lets the DOF map number them in whatever
/// order it chooses.
struct BoundSurface {
    std::string name;
    int tag = -1;
    std::vector<int> faces;
    std::vector<int> nodes;
    std::vector<int> edges;
};

/// One port, resolved: its faces, its orientation, and which side of a cut
/// is which.
struct BoundPort {
    std::string name;
    PortType type = PortType::BoundaryCurrent;

    /// The excitation, phase already applied. Feeds the port row's
    /// right-hand side (current port) or the prescribed potential on the
    /// live side (voltage port) -- with no sign factor, per problem.hpp's
    /// orientation note.
    std::complex<double> amplitude;

    BoundSurface surface;  ///< merged over every name the port listed

    /// **d**, resolved. For an internal port this runs from the grounded
    /// minus side to the live plus side; for a boundary port it is the
    /// inward normal of the terminal (positive current enters the domain).
    Vec3 direction;

    /// True if `d` came from the input file's `current_direction`. False
    /// means it was derived, and the port's reported I and V then carry an
    /// arbitrary -- though mesh-stable -- sign.
    bool direction_from_hint = false;

    /// Internal ports only, parallel to `surface.faces`: which adjacent tet
    /// of each cut face lies on the plus side (the side `direction` points
    /// into). -1 throughout for a boundary port, which has only one side.
    std::vector<int> plus_side_tet;

    /// Internal ports only: edges belonging to exactly one cut face -- the
    /// rim. Kept because the rim is where an incomplete cut shows up, and
    /// where the full-wave delta gap lives.
    std::vector<int> rim_edges;

    /// Internal ports only: every tet sharing at least one node with the
    /// cut, ascending, with `tet_side` giving +1 (plus side) or -1 (minus
    /// side) for each.
    ///
    /// `plus_side_tet` above covers only tets with a *face* on the cut.
    /// This covers the rest -- a tet touching at one node or one edge still
    /// has to know which side it sees that node from, because Phi is 0 on
    /// the minus side and the port's unknown on the plus side, and a tet
    /// reads one or the other.
    ///
    /// The classification is a plane test, not a flood fill: the cut is
    /// planar (binding rejects it otherwise) and is built from mesh faces,
    /// so no tet straddles it and every tet's centroid is strictly on one
    /// side. Exact, and no tolerance to tune.
    std::vector<int> touching_tets;
    std::vector<signed char> tet_side;

    bool is_internal() const { return aphi_solver::is_internal(type); }

    /// +1 (plus side), -1 (minus side), or 0 if tet `t` does not touch this
    /// cut at all. Binary search over `touching_tets`.
    int side_of_tet(int t) const;
};

/// One electrically connected conducting region: the tets a current could
/// flow through without leaving the conductors.
///
/// A path is **not** the same as a body. Two bodies that touch -- a via
/// meeting a trace -- are one path made of two materials, and it is the
/// path, not the body, that needs a potential reference. Conversely one
/// body meshed as two separate islands is two paths.
///
/// Built from `sigma > 0` alone, never from `phi_tet`: a conduction path is
/// a property of the material, not of which formulation is being solved.
/// (An earlier draft built it from `phi_tet`, which is the conductors only
/// at DC -- at full wave it would have merged the entire domain into one
/// "path".)
struct ConductionPath {
    std::vector<int> tets;            ///< ascending
    std::vector<int> bodies;          ///< indices into BoundProblem::bodies
    std::vector<std::string> ports;   ///< ports touching this path
    int reference_count = 0;          ///< ports on it that fix a potential

    /// True when nothing fixes this path's potential and no port touches
    /// it. At DC its Phi would be determined only up to a constant, so one
    /// node is pinned; the value is physically meaningless, since a
    /// floating conductor carries no DC current.
    bool is_floating = false;
    int pin_node = -1;  ///< the node to pin when floating, else -1
};

/// A problem bound to a mesh: everything the DOF map needs, and nothing it
/// would have to re-derive.
struct BoundProblem {
    Problem problem;  ///< the source description, kept for its settings

    std::vector<BoundBody> bodies;
    std::vector<int> body_of_tet;  ///< size num_tets; index into `bodies`
    std::vector<BoundPort> ports;

    /// Which tets Phi lives on, per `phi_on_conductors_only`. All true in
    /// the full-wave regime; the conductors only at DC and in the reduced
    /// variant, where Phi's equation is identically empty elsewhere.
    std::vector<bool> phi_tet;

    /// The conducting regions, and which one each tet belongs to (-1 where
    /// sigma == 0). Always built, in every regime.
    std::vector<ConductionPath> conduction_paths;
    std::vector<int> path_of_tet;

    /// The tree-cotree gauge for this mesh and boundary condition, run here
    /// rather than left to the caller because the DOF map cannot classify a
    /// single edge without it: every edge is FREE, TREE (a = 0, the gauge)
    /// or DIRICHLET (a = 0, from n x A = 0), and two of the three come from
    /// this result.
    ///
    /// One tree serves both regimes. The mask below comes from the boundary
    /// condition, which does not depend on frequency, and at DC the
    /// magnetostatic stage still needs a gauge because curl-curl is
    /// singular without one.
    TreeCotreeResult gauge;

    /// The n x A = 0 edge mask the gauge was built from. Kept because the
    /// DOF map needs it too, and re-deriving it is how the two could drift
    /// apart.
    std::vector<bool> dirichlet_edge;

    std::vector<std::string> warnings;
};

/// Scales every node coordinate to metres. Call **once**, after reading the
/// mesh and before binding: a Gmsh file carries bare numbers, and
/// `[mesh] length_unit` is what says which unit they are in.
///
/// Kept separate from `bind_to_mesh` rather than folded into it, because a
/// function that silently rescales its input is surprising, and calling it
/// twice would be a quiet factor-of-1000 error.
void scale_mesh_to_metres(Mesh& mesh, LengthUnit unit);

/// Resolves `problem` against `mesh` and checks everything that needs the
/// mesh. See `docs/INPUT_FILE_PLAN.md` Sec. 3.2 for the complete
/// list; in outline, it rejects a name that is not in the mesh, a volume no
/// body claims (or two do), a boundary port whose faces are interior (or an
/// internal port whose faces are not), a non-planar cut, two terminals
/// sharing a node, an incomplete cut at DC, and a conductor with no
/// potential reference or with two.
///
/// Throws InputError, carrying the line of the `[Body ...]` or `[port ...]`
/// section at fault so the message points into the input file rather than
/// merely describing the problem.
BoundProblem bind_to_mesh(const Problem& problem, const Mesh& mesh);

}  // namespace aphi_solver
