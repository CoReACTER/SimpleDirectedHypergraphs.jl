"""
    struct SnodeDistanceKKHeuristic(sources::Set{Int}, targets::Set{Int}) <: AbstractDistance

Represent a distance between a set of *source* nodes and a set of *target* nodes in a of directed hypergraph `h`,
where a distance in this context is defined as the length of the shortest hyperpath that reaches all `t ∈ targets`.

This approach uses the *heuristic* shortest hyperpath algorithm of Krieger & Kececioglu 2022
(DOI: 10.1186/s13015-022-00217-9). The Krieger & Kececioglu heuristic algorithm commonly (but not always) returns the
true optimal (shortest) hyperpath; this algorithm can thus be thought of as an upper bound of the true distance between
the sinks and the source. To calculate the exact (minimal) distance, instead use `SnodeDistanceKKExact`.
"""
struct SnodeDistanceKKHeuristic <: AbstractDistance
    sources::Set{Int}
    targets::Set{Int}
end

"""
    struct SnodeDistanceKKILP(sources::Set{Int}, targets::Set{Int}) <: AbstractDistance

Represent a distance between a set of *source* nodes and a set of *target* nodes in a of directed hypergraph `h`,
where a distance in this context is defined as the length of the shortest hyperpath that reaches all `t ∈ targets`.

This approach uses the *exact* shortest hyperpath algorithm of Krieger & Kececioglu 2023 (DOI: 10.1089/cmb.2023.0242).
The Krieger & Kececioglu algorithm solves an integer linear programming formalism related to hypergraph cutting.

For a lower-cost heuristic algorith, see `SnodeDistanceKKHeuristic`.
"""
struct SnodeDistanceKKILP <: AbstractDistance
    sources::Set{Int}
    targets::Set{Int}
end

function SimpleHypergraphs.distance(
    h::H,
    distance_method::SnodeDistanceKKHeuristic,
    cost_function::Function,
    hyperedge_weights::Vector{T}
) where {H <: AbstractDirectedHypergraph}

    return cost_function()

end