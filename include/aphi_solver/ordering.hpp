#pragma once

#include <string>
#include <vector>

#include "aphi_solver/sparsity.hpp"

namespace aphi_solver {

/// Fill-reducing reordering of the unknowns, before factorization.
/// `docs/SOLVER_PLAN.md` §3 and step 1 of §9.
///
/// Factorizing in the DOF map's natural order -- `[a | phi | V]`, which is edge
/// and node numbering inherited from the mesh file -- produces far more
/// nonzeros in the factor than the matrix has, and both memory and flops scale
/// with that. A permutation is the largest single lever on a direct solve.
///
/// Three of them, because a measurement needs a baseline.
enum class Ordering {
    /// No permutation. **The control**: if a real ordering does not beat this
    /// substantially, it is not working.
    Natural,

    /// Reverse Cuthill-McKee (Cuthill & McKee 1969; the reversal is George's).
    /// Reduces bandwidth and profile rather than fill directly, which is
    /// usually worse than minimum degree for a 3D factorization -- it is here
    /// as a cheap sanity check that ordering helps at all, and as a fallback.
    ReverseCuthillMcKee,

    /// Approximate minimum degree (Amestoy, Davis & Duff 1996, 2004). The
    /// standard fill-reducing ordering and the one to use.
    /// **Not implemented yet** -- step 3 of `SOLVER_PLAN.md` §9.
    ApproximateMinimumDegree
};

/// The input-file keyword for an ordering, and back. `from` returns false for
/// anything unrecognised rather than guessing.
const char* ordering_keyword(Ordering ordering);
bool ordering_from_keyword(const std::string& word, Ordering& out);

/// A permutation of the unknowns, kept in both directions because both are
/// needed and confusing them is the classic error in this area: the matrix is
/// permuted with one and the right-hand side with the other.
///
/// `perm[i]` is the OLD index that now sits at NEW position `i`.
/// `iperm[j]` is the NEW position of OLD index `j`.
/// So `iperm[perm[i]] == i` and `perm[iperm[j]] == j`.
struct Permutation {
    std::vector<int> perm;
    std::vector<int> iperm;

    int size() const { return static_cast<int>(perm.size()); }

    /// True when this really is a permutation: both arrays the same length,
    /// every entry in range, and the two mutually inverse. Cheap, and worth
    /// asserting -- a non-bijective "permutation" loses or duplicates
    /// unknowns, which downstream shows up as a plausible wrong answer rather
    /// than a crash.
    bool is_valid() const;

    /// The permutation that undoes this one.
    Permutation inverse() const { return Permutation{iperm, perm}; }

    /// The identity on `n` unknowns.
    static Permutation identity(int n);
};

/// Computes the permutation for `pattern`.
///
/// The adjacency is built by treating every stored entry `(i, j)` as an
/// undirected edge, so this accepts a full pattern, an upper-triangle one
/// (`SparsityStorage::UpperTriangle`), or an unsymmetric one -- in the last
/// case it orders the pattern of `A + Aᵀ`, which is what an unsymmetric
/// factorization wants anyway. One implementation therefore serves both solver
/// paths.
///
/// Isolated unknowns are handled rather than assumed away: a voltage port's
/// constraint row holds nothing but its own diagonal
/// (`docs/ASSEMBLY_PLAN.md` §10), so the adjacency graph really does have
/// isolated vertices, and more than one connected component is normal.
///
/// Throws std::invalid_argument for `ApproximateMinimumDegree` until step 3 of
/// `SOLVER_PLAN.md` §9 implements it, rather than silently falling back to
/// another ordering and reporting a fill figure that belongs to something else.
Permutation compute_ordering(const SparsityPattern& pattern, Ordering ordering);

/// `P A Pᵀ`: the same nonzeros, at permuted positions. Entry `(i, j)` of the
/// result is entry `(perm[i], perm[j])` of the input.
///
/// An upper-triangle pattern stays an upper-triangle pattern: a permutation
/// generally moves an entry across the diagonal, so each pair is normalised
/// back to `col >= row`. The nonzero COUNT is preserved either way, which is
/// one of the checks in `SOLVER_PLAN.md` §8.
SparsityPattern permute_pattern(const SparsityPattern& pattern, const Permutation& p);

/// What Reverse Cuthill-McKee actually reduces, so its effect can be measured
/// rather than asserted.
///
/// `bandwidth` is the largest `|i - j|` over stored entries; `profile` is the
/// sum over rows of the distance from the row's first stored column to the
/// diagonal. Both are computed over the pattern as given -- for an
/// upper-triangle pattern they describe that triangle, which is the right
/// thing since the other one mirrors it.
struct BandwidthStats {
    int bandwidth = 0;
    long long profile = 0;
    double mean_bandwidth = 0.0;  ///< profile / rows, the quantity that tracks fill better
};

BandwidthStats bandwidth_stats(const SparsityPattern& pattern);

}  // namespace aphi_solver
