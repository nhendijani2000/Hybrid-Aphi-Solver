#pragma once

#include <complex>
#include <string>
#include <vector>

#include "aphi_solver/input_file.hpp"  // InputError
#include "aphi_solver/mesh.hpp"
#include "aphi_solver/problem.hpp"

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

    bool is_internal() const { return aphi_solver::is_internal(type); }
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
/// mesh. See `Claude outputs/input_file_plan.md` Sec. 3.2 for the complete
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
