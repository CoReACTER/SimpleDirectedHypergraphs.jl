"""
    SnodeDistanceKKHeuristic <: AbstractDistance

Constructor: SnodeDistanceKKHeuristic(sources::Set{Int}, targets::Set{Int})

Represent a distance between a set of *source* nodes (`sources`) and a set of *target* nodes (`targets`) in a directed
hypergraph, where a distance in this context is defined as the length of the shortest hyperpath that reaches all
`t ∈ targets`.

This approach uses the *heuristic* shortest hyperpath algorithm of Krieger & Kececioglu 2022
(DOI: 10.1186/s13015-022-00217-9). The Krieger & Kececioglu heuristic algorithm commonly (but not always) returns the
true optimal (shortest) hyperpath; this algorithm can thus be thought of as an upper bound of the true distance between
the sinks and the source. To calculate the exact (minimal) distance, instead use `SnodeDistanceKKExact`.

**Arguments**

* `sources` : A `Set{Int}` of vertex indices that the pathfinding algorithm starts with
* `targets` : A `Set{Int}` of vertex indices that must be reached

"""
struct SnodeDistanceKKHeuristic <: AbstractDistance
    sources::Set{Int}
    targets::Set{Int}
end

"""
    SnodeDistanceKKILP <: AbstractDistance

Constructor: SnodeDistanceKKILP(sources::Set{Int}, targets::Set{Int})

Represent a distance between a set of *source* nodes (`sources`) and a set of *target* nodes (`targets`) in a directed
hypergraph, where a distance in this context is defined as the length of the shortest hyperpath that reaches all
`t ∈ targets`.

This approach uses the *exact* shortest hyperpath algorithm of Krieger & Kececioglu 2023 (DOI: 10.1089/cmb.2023.0242).
The Krieger & Kececioglu algorithm solves an integer linear programming formalism related to hypergraph cutting.

For a lower-cost heuristic algorithm, see `SnodeDistanceKKHeuristic`.

**Arguments**

* `sources` : A `Set{Int}` of vertex indices that the pathfinding algorithm starts with
* `targets` : A `Set{Int}` of vertex indices that must be reached

"""
struct SnodeDistanceKKILP <: AbstractDistance
    sources::Set{Int}
    targets::Set{Int}
end

"""
    distance(
        h::H,
        distance_method::SnodeDistanceKKHeuristic,
        dihyperedge_weights::AbstractVector{T};
    ) where {H <: AbstractDirectedHypergraph, T <: Real}

Return the shortest distance between the `distance_method.sources` and the
`distance_method.targets`, assuming tha the cost (i.e., distance) between two vertices is the sum
of the (nonnegative) `hyperedge_weights` of the hyperedges in the shortest path between the two sets
of vertices. If there is no path between the `sources` and `targets`, returns `typemax(T)`.

Here, the heuristic directed hypergraph pathfinding algorithm of Krieger & Kececioglu (2022),
DOI: 10.1186/s13015-022-00217-9 is used to find the shortest path. This algorithm is not guaranteed
to find the optimal pathway, but in practice, it produces the optimal pathway approximately 99% of
the time.

**Arguments**

* `h` : Dihypergraph
* `distance_method` : A `SnodeDistanceKKHeuristic` instance
* `dihyperedge_weights` : `AbstractVector` of real-valued weights with length equal to `nhe(h)`

**Returns**

If no path can be found, returns `typemax(T)`. Otherwise, returns the sum of the weights of the
dihyperedges included in the (estimated) shortest path.

"""
function SimpleHypergraphs.distance(
        h::H,
        distance_method::SnodeDistanceKKHeuristic,
        dihyperedge_weights::AbstractVector{T}
    ) where {H <: AbstractDirectedHypergraph, T <: Real}

    path = shortest_hyperpath_kk_heuristic(
        h,
        distance_method.sources,
        distance_method.targets,
        dihyperedge_weights
    )

    if length(path) == 0
        return typemax(T)
    end

    return sum(
        dihyperedge_weights[e]
            for e in path
    )
end

"""
    distance(
        h::H,
        distance_method::SnodeDistanceKKILP,
        dihyperedge_weights::AbstractVector{T}
    ) where {H <: AbstractDirectedHypergraph}

Return the shortest distance between the `distance_method.sources` and the `distance_method.targets`,
assuming tha the cost (i.e., distance) between two vertices is the sum of the (nonnegative) `hyperedge_weights`
of the hyperedges in the shortest path between the two sets of vertices. If there is no path between the `sources`
and `targets`, returns `typemax(T)`.

Here, the *exact* shortest hyperpath algorithm of Krieger & Kececioglu 2023 (DOI: 10.1089/cmb.2023.0242) is used.
The Krieger & Kececioglu algorithm solves an integer linear programming formalism related to hypergraph cutting.

**Arguments**

* `h` : Dihypergraph
* `distance_method` : A `SnodeDistanceKKILP` instance
* `dihyperedge_weights` : `AbstractVector` of real-valued weights with length equal to `nhe(h)`

**Returns**

If no path can be found, returns `typemax(T)`. Otherwise, returns the sum of the weights of the
dihyperedges included in the shortest path.

"""
function SimpleHypergraphs.distance(
        h::H,
        distance_method::SnodeDistanceKKILP,
        dihyperedge_weights::AbstractVector{T}
    ) where {H <: AbstractDirectedHypergraph, T <: Real}

    path = shortest_hyperpath_kk_ilp(
        h,
        distance_method.sources,
        distance_method.targets,
        dihyperedge_weights
    )

    if length(path) == 0
        return typemax(T)
    end

    return sum(
        dihyperedge_weights[e]
            for e in path
    )
end

"""
    diameter(
        h::H,
        distance_method::SnodeDistanceKKHeuristic,
    ) where {H <: AbstractDirectedHypergraph, T <: Real}

Return the diameter of a hypergraph `hg` (maximum number of hyperedges required to go between any
two vertices) based on Krieger & Kececioglu's heuristic algorithm (DOI: 10.1186/s13015-022-00217-9).
If there exist some vertices not reachable from other vertices (i.e., if `h` is not strongly
connected), then this will return `Inf64`.

TODO: should consider starting from more than one vertex; prerequisite problem

**Arguments**

* `h` : Dihypergraph
* `distance_method` : An `SnodeDistanceKKHeuristic` instance (note: initial `sources` and `targets`
    do not matter)

**Returns**

If `h` is not strongly connected, returns `Inf64`; otherwise, returns the maximum distance from any
one source vertex to any one target vertex.

"""
function Graphs.diameter(
        h::H,
        _::SnodeDistanceKKHeuristic,
    ) where {H <: AbstractDirectedHypergraph}
    if length(get_strongly_connected_components(h)) != 1
        return Inf64
    end

    # Since we're interested in the number of steps, assign each hyperedge a uniform weight of `1`
    edge_weights = ones(nhe(h))

    max_dist = 0.0
    for i in 1:nhv(h)
        for j in (i + 1):nhv(h)
            dist = SimpleHypergraphs.distance(
                h,
                SnodeDistanceKKHeuristic(Set{Int}(i), Set{Int}(j)),
                edge_weights
            )
            if dist > max_dist
                max_dist = dist
            end
        end
    end

    return max_dist
end

"""
    diameter(
        hg::H,
        distance_method::SnodeDistanceKKILP,
    ) where {H <: AbstractDirectedHypergraph, T <: Real}

Return the diameter of a hypergraph `h` (maximum number of hyperedges required to go between any
two vertices) based on Krieger & Kececioglu's exact integer linear programming algorithm
(DOI: 10.1089/cmb.2023.0242). If there exist some vertices not reachable from other vertices
(i.e., if `h` is not strongly connected), then this will return `Inf64`.

TODO: should consider starting from more than one vertex; prerequisite problem

**Arguments**

* `h` : Dihypergraph
* `distance_method` : An `SnodeDistanceKKILP` instance (note: initial `sources` and `targets`
    do not matter)

**Returns**

If `h` is not strongly connected, returns `Inf64`; otherwise, returns the maximum distance from any
one source vertex to any one target vertex.

"""
function Graphs.diameter(
        hg::H,
        _::SnodeDistanceKKILP,
    ) where {H <: AbstractDirectedHypergraph}
    if length(get_strongly_connected_components(hg)) != 1
        return Inf64
    end

    # Since we're interested in the number of steps, assign each hyperedge a uniform weight of `1`
    edge_weights = ones(nhe(hg))

    max_dist = 0
    for i in 1:nhv(hg)
        for j in (i + 1):nhv(hg)
            dist = SimpleHypergraphs.distance(
                hg,
                SnodeDistanceKKILP(Set{Int}(i), Set{Int}(j)),
                edge_weights
            )
            if dist > max_dist
                max_dist = dist
            end
        end
    end

    return max_dist
end
