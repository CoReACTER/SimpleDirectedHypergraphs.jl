"""
    forward_reachable(
        h::H,
        source::Int,
    ) where {H <: AbstractDirectedHypergraph}

Traverses a dihypergraph `h` starting from vertex with index `source` to determine all other vertices
and dihyperedges that are reachable, following dihyperedges along their forward direction (i.e., from
tail to head).

**Arguments**

* `h` : Dihypergraph
* `source` : A vertex index from which to begin the search

**Returns**

`(vs, es)`, where `vs` is a `Set{Int}` of reached vertex indices and `es` is a `Set{Int}` of
reached dihyperedge indices

"""
function forward_reachable(
        h::H,
        source::Int,
    ) where {H <: AbstractDirectedHypergraph}
    # Priority queue of reached vertices
    Q = Queue{Int}()
    enqueue!(Q, source)

    reached_vs = BitVector(falses(nhv(h)))
    reached_vs[source] = true

    hes_tail_count = length.(keys.(h.hg_tail.he2v))

    # Which vertices/hyperedges have been reached?
    vs = Set{Int}()
    es = Set{Int}()

    while length(Q) > 0
        v = dequeue!(Q)
        push!(vs, v)

        for out_e in keys(h.hg_tail.v2he[v])
            # Following pseudocode exactly. This feels awkward; how slow would it be to just query reached_vs?
            hes_tail_count[out_e] -= 1

            if hes_tail_count[out_e] == 0
                push!(es, out_e)

                for w in keys(h.hg_head.he2v[out_e])
                    if !reached_vs[w]
                        enqueue!(Q, w)
                        reached_vs[w] = true
                    end
                end
            end
        end
    end

    # Return reached vertices and dihyperedges
    return (vs, es)
end

"""
    backward_traceable(
        h::H,
        target::Int,
    ) where {H <: AbstractDirectedHypergraph}

Traverses a dihypergraph `h` starting from vertex with index `target` to determine all other
vertices and dihyperedges that are reachable, following dihyperedges along their reverse direction
(i.e., from head to tail).

**Arguments**

* `h` : Dihypergraph
* `source` : A vertex index from which to begin the search

**Returns**

`(vs, es)`, where `vs` is a `Set{Int}` of reached vertex indices and `es` is a `Set{Int}` of
reached dihyperedge indices

"""
function backward_traceable(
        h::H,
        target::Int,
    ) where {H <: AbstractDirectedHypergraph}
    # Priority queue of reached vertices
    Q = Queue{Int}()
    enqueue!(Q, target)

    reached_vs = BitVector(falses(nhv(h)))
    reached_vs[target] = true

    marked_hes = BitVector(falses(nhe(h)))

    # Which vertices/hyperedges have been reached?
    vs = Set{Int}()
    es = Set{Int}()

    while length(Q) > 0
        v = dequeue!(Q)
        push!(vs, v)

        for in_e in keys(h.hg_head.v2he[v])
            if !marked_hes[in_e]
                push!(es, in_e)
                marked_hes[in_e] = true

                for w in keys(h.hg_tail.he2v[in_e])
                    if !reached_vs[w]
                        enqueue!(Q, w)
                        reached_vs[w] = true
                    end
                end
            end
        end
    end

    # Return reached vertices and dihyperedges
    return (vs, es)
end


"""
    shortest_hyperpath_kk_heuristic(
        h::H,
        source::Int,
        target::Int,
        dihyperedge_weights::AbstractVector{T}
    ) where {H <: AbstractDirectedHypergraph, T <: Real}

    shortest_hyperpath_kk_heuristic(
        h::DirectedHypergraph{T, V, E, D},
        source::Int,
        targets::Set{Int},
        dihyperedge_weights::AbstractVector{T}
    ) where {T <: Real, V, E, D <: AbstractDict{Int,T}}

    shortest_hyperpath_kk_heuristic(
        h::DirectedHypergraph{T, V, E, D},
        sources::Set{Int},
        target::Int,
        dihyperedge_weights::AbstractVector{T}
    ) where {T <: Real, V, E, D <: AbstractDict{Int,T}}

    shortest_hyperpath_kk_heuristic(
        h::DirectedHypergraph{T, V, E, D},
        sources::Set{Int},
        targets::Set{Int},
        dihyperedge_weights::AbstractVector{T}
    ) where {T <: Real, V, E, D <: AbstractDict{Int,T}}

Implements the heuristic dihypergraph pathfinding algorithm of Krieger & Kececioglu (2022),
DOI: 10.1186/s13015-022-00217-9. This algorithm is not guaranteed to find the optimal pathway from
`source`/`sources` to `target`/`targets`, but in practice, it produces the optimal pathway
approximately 99% of the time.
    
Note that, ostensibly, this algorithm only works for single-source, single-sink pathfinding (i.e.,
with a single `source` and a single `target`). If the user provides multiple `sources` and/or
multiple `targets`, the multi-source/multi-sink problem will be reformulated as a single-source,
single-sink problem by adding a *metasource* vertex (connected to all source vertices by a single,
0-cost dihyperedge) and/or *metatarget* vertex (connected to all target vertices by a single, 0-cost
hyperedge).

**Arguments**

* `h` : Dihypergraph
* `source`/`sources` : Either the index of a single source vertex or a `Set{Int}` of vertex indices
* `target`/`targets` : Either the index of a single target vertex or a `Set{Int}` of vertex indices
* `dihyperedge_weights` : A `Vector` of nonnegative real-valued weights of length equal to `nhe(h)`

**Returns**

`path`, a `Set{Int}` of dihyperedge vertices constituting an estimated shortest path from
`source`/`sources` to `target`/`targets`

"""
function shortest_hyperpath_kk_heuristic(
        h::H,
        source::Int,
        target::Int,
        dihyperedge_weights::AbstractVector{T}
    ) where {H <: AbstractDirectedHypergraph, T <: Real}

    reached_vs = BitVector(falses(nhv(h)))
    reached_vs[source] = true

    marked_hes = BitVector(falses(nhe(h)))
    removed_hes = BitVector(falses(nhe(h)))

    hyperedge_inedges = [Set{Int}() for _ in 1:nhe(h)]
    hes_tail_count = length.(keys.(h.hg_tail.he2v))

    hyperedge_costs = fill(typemax(T), nhe(h))

    hyperedge_heap_points = Vector{Union{Nothing, Int}}(nothing, nhe(h))

    # Verify that the target can be reached
    fr = forward_reachable(h, source)
    @assert target ∈ fr[1]

    # Doubly reachable dihyperedges
    dr_hes = sort!(
        collect(
            intersect(
                fr[2],
                backward_traceable(h, target)[2]
            )
        )
    )

    # Eliminate non-doubly reachable dihyperedges
    h_copy = deepcopy(h)
    h_copy.hg_tail[:, InvertedIndices.Not(dr_hes)] .= nothing
    h_copy.hg_head[:, InvertedIndices.Not(dr_hes)] .= nothing

    # Min-heap for dihyperedges
    Hmin = MutableBinaryMinHeap{Tuple{T, Int}}()
    for out_e in keys(h.hg_tail.v2he[source])
        # If only the source is needed for this dihyperedge
        if length(h.hg_tail.he2v[out_e]) == 1
            hyperedge_heap_points[out_e] = push!(Hmin, (dihyperedge_weights[out_e], out_e))
        end
    end

    while length(Hmin) > 0
        e = pop!(Hmin)[2]
        removed_hes[e] = true

        path = short_hyperpath_vhe(h, source, e, hyperedge_inedges, hyperedge_costs)
        hyperedge_costs[e] = sum(dihyperedge_weights[x] for x in path)

        out_edges = Set{Int}()
        for v in keys(h.hg_head.he2v[e])
            for f in keys(h.hg_tail.v2he[v])
                if !marked_hes[f]
                    if !reached_vs[v]
                        hes_tail_count[f] -= 1
                    end

                    if hes_tail_count[f] == 0
                        push!(out_edges, f)
                        marked_hes[f] = true
                    end
                end
            end
            reached_vs[v] = true
        end

        for f in out_edges
            marked_hes[f] = false
        end

        for f in out_edges
            push!(hyperedge_inedges[f], e)
            if !isnothing(hyperedge_heap_points[f]) && !removed_hes[f]
                update!(
                    Hmin,
                    hyperedge_heap_points[f],
                    (
                        sum(
                            dihyperedge_weights[x]
                                for x in short_hyperpath_vhe(h, source, f, hyperedge_inedges, hyperedge_costs)
                        ),
                        f,
                    )
                )
            elseif isnothing(hyperedge_heap_points[f]) && hes_tail_count[f] == 0
                hyperedge_heap_points[f] = push!(
                    Hmin,
                    (
                        sum(
                            dihyperedge_weights[x]
                                for x in short_hyperpath_vhe(h, source, f, hyperedge_inedges, hyperedge_costs)
                        ),
                        f,
                    )
                )
            end
        end
    end

    path = Set{Int}()
    cost = typemax(T)
    for in_e in keys(h.hg_head.v2he[target])
        if !isnothing(hyperedge_heap_points[in_e])
            p = short_hyperpath_vhe(h, source, in_e, hyperedge_inedges, hyperedge_costs)
            cost_p = sum(dihyperedge_weights[e] for e in p)
            if cost_p < cost
                path = p
                cost = cost_p
            end
        end
    end

    return path
end

function shortest_hyperpath_kk_heuristic(
        h::DirectedHypergraph{T, V, E, D},
        source::Int,
        targets::Set{Int},
        dihyperedge_weights::AbstractVector{S}
    ) where {S <: Real, T <: Real, V, E, D <: AbstractDict{Int, T}}
    h_copy = deepcopy(h)

    # Add a single "metatarget" vertex to reformulate as single-source, single-sink pathfinding problem
    # The dihyperedge from the targets to the metatarget will have a cost of 0 associated with it
    metatarget = SimpleHypergraphs.add_vertex!(h_copy)
    meta_he = SimpleHypergraphs.add_hyperedge!(
        h_copy;
        vertices_tail = D(x => convert(T, 0) for x in targets),
        vertices_head = D(metatarget => convert(T, 0))
    )

    path = shortest_hyperpath_kk_heuristic(
        h_copy,
        source,
        metatarget,
        vcat(dihyperedge_weights, convert(S, 0))
    )

    # Remove the fictitious dihyperedge from the targets to the metatarget
    return setdiff(path, Set{Int}(meta_he))
end

function shortest_hyperpath_kk_heuristic(
        h::DirectedHypergraph{T, V, E, D},
        sources::Set{Int},
        target::Int,
        dihyperedge_weights::AbstractVector{S}
    ) where {S <: Real, T <: Real, V, E, D <: AbstractDict{Int, T}}
    h_copy = deepcopy(h)

    # Add a single "metasource" vertex to reformulate as single-source, single-sink pathfinding problem
    # The dihyperedge from the metasource to the sources will have a cost of 0 associated with it
    metasource = SimpleHypergraphs.add_vertex!(h_copy)
    meta_he = SimpleHypergraphs.add_hyperedge!(
        h_copy;
        vertices_tail = D(metasource => convert(T, 0)),
        vertices_head = D(x => convert(T, 0) for x in sources)
    )

    path = shortest_hyperpath_kk_heuristic(
        h_copy,
        metasource,
        target,
        vcat(dihyperedge_weights, convert(S, 0))
    )

    # Remove the fictitious dihyperedge from the metasource to the sources
    return setdiff(path, Set{Int}(meta_he))
end

function shortest_hyperpath_kk_heuristic(
        h::DirectedHypergraph{T, V, E, D},
        sources::Set{Int},
        targets::Set{Int},
        dihyperedge_weights::AbstractVector{S}
    ) where {S <: Real, T <: Real, V, E, D <: AbstractDict{Int, T}}
    h_copy = deepcopy(h)

    # Add a single "metasource" vertex to reformulate as single-source, single-sink pathfinding problem
    # The dihyperedge from the metasource to the sources will have a cost of 0 associated with it
    metasource = SimpleHypergraphs.add_vertex!(h_copy)
    meta_he_source = SimpleHypergraphs.add_hyperedge!(
        h_copy;
        vertices_tail = D(metasource => convert(T, 0)),
        vertices_head = D(x => convert(T, 0) for x in sources)
    )

    # Add a single "metatarget" vertex to reformulate as single-source, single-sink pathfinding problem
    # The dihyperedge from the targets to the metatarget will have a cost of 0 associated with it
    metatarget = SimpleHypergraphs.add_vertex!(h_copy)
    meta_he_target = SimpleHypergraphs.add_hyperedge!(
        h_copy;
        vertices_tail = D(x => convert(T, 0) for x in targets),
        vertices_head = D(metatarget => convert(T, 0))
    )

    path = shortest_hyperpath_kk_heuristic(
        h_copy,
        metasource,
        metatarget,
        vcat(dihyperedge_weights, [convert(S, 0), convert(S, 0)])
    )

    # Remove fictitious dihyperedges
    return setdiff(path, Set{Int}([meta_he_source, meta_he_target]))
end

"""
    short_hyperpath_vhe(
        h::H,
        v::Int,
        he::Int,
	he_inedges::Vector{Set{Int}},
	he_costs::AbstractVector{T}
    ) where {H <: AbstractDirectedHypergraph, T <: Real}

Obtain a (relatively, but not necessarily optimally) short hyperpath in dihypergraph `h` from a
vertex `v` to a dihyperedge `he`, `short_hyperpath_vhe` uses a greedy algorithm to first select
dihyperedges for a superpath and then prune unnecessary dihyperedges to achieve a (generally
shorter) hyperpath. 

**Arguments**

* `h` : Dihypergraph
* `v` : Vertex index
* `he` : Dihyperedge index
* `he_inedges` : `Vector` containing the incoming dihyperedges to each dihyperedge in `h` (as
    `Set{Int}`)
* `he_costs` : A `Vector` of dihyperedge costs

**Returns**

`path`, a hyperpath from `v` to `he`, represented as a `Set{Int}` of dihyperedge indices

"""
function short_hyperpath_vhe(
        h::H,
        v::Int,
        he::Int,
        he_inedges::Vector{Set{Int}},
        he_costs::AbstractVector{T}
    ) where {H <: AbstractDirectedHypergraph, T <: Real}
    marked_hes = BitVector(falses(nhe(h)))

    Q = Queue{Int}()
    for e in he_inedges[he]
        enqueue!(Q, e)
        marked_hes[e] = true
    end

    superpath = Set{Int}(he)
    path = Set{Int}(he)

    # Construct (likely redundant) superpath by backtracking from target
    while length(Q) > 0
        e = dequeue!(Q)
        push!(superpath, e)

        for f in he_inedges[e]
            if !marked_hes[f]
                enqueue!(Q, f)
                marked_hes[f] = true
            end
        end
    end

    # TODO: try to be more clever about this
    superpath = sort(collect(superpath), by = x -> he_costs[x], rev = true)
    h_copy = deepcopy(h)
    # Eliminate all dihyperedges not on superpath
    h_copy.hg_tail[:, InvertedIndices.Not(superpath)] .= nothing
    h_copy.hg_head[:, InvertedIndices.Not(superpath)] .= nothing

    # Remove target from superpath; does not make sense to remove target in the next stage
    filter!(x -> x != he, superpath)

    # Try to minimize the size of the path by eliminating unnecessary dihyperedges
    for e in superpath
        h_copy.hg_tail[:, e] .= nothing
        h_copy.hg_head[:, e] .= nothing

        # Only if dihyperedge is essential for reaching target,
        if !is_reachable(h_copy, v, he, :dihyperedge)
            # Restore dihyperedge to hypergraph copy
            h_copy[:, e] .= h[:, e]
            push!(path, e)
        end
    end

    return path
end

"""
    is_reachable(
        h::H,
        source::Int,
        target::Int,
        target_type::Symbol
    ) where {H<:AbstractDirectedHypergraph}

Use `forward_reachable` to determine if `target` (either a vertex index, if
`target_type === :vertex` or a dihyperedge index if `target_type === :dihyperedge`) can be reached
from vertex with index `source` in dihypergraph `h`.

**Arguments**

* `h` : Dihypergraph
* `source` : Vertex index
* `target` : Either a vertex index or a dihyperedge index, depending on `target_type`
* `target_type` : A `Symbol` (either `:vertex` or `dihyperedge`) indicating what type of object
    `target` represents

**Returns**

A `Bool` indicating if `target` is reachable from `source`

"""
function is_reachable(
        h::H,
        source::Int,
        target::Int,
        target_type::Symbol
    ) where {H <: AbstractDirectedHypergraph}
    @assert target_type ∈ [:vertex, :dihyperedge] "`target_type` must be :vertex or :dihyperedge"

    fr = forward_reachable(h, source)

    #TODO: make a short-circuiting version of this. Can't for the life of me figure out why this isn't working...
    if target_type === :vertex
        return target ∈ fr[1]
    else
        return target ∈ fr[2]
    end
end

"""
    get_hyperpath(
        h::H,
        source::Int,
        target::Int,
        out::Set{Int}
    ) where {H <: AbstractDirectedHypergraph}

If one exists, obtain a hyperpath in dihypergraph `h` from a source vertex with index `source` to a
target vertex with index `target`. The hyperpath cannot include any dihyperedge with index included
in the set `out`.

**Arguments**

* `h` : Dihypergraph
* `source` : Vertex index that serves as the start of the path
* `target` : Vertex index that serves as the target of the path
* `out` : A `Set{Int}` of dihyperedge indices to exclude in path-making

**Returns**

`path`, a `Set{Int}` of dihyperedge indices constituting a valid hyperpath from `source` to `target`

"""
function get_hyperpath(h::H, source::Int, target::Int, out::Set{Int}) where {H <: AbstractDirectedHypergraph}
    # Remove excluded dihyperedges
    h_copy = deepcopy(h)
    inds = sort(collect(out))
    h_copy.hg_tail[:, inds] .= nothing
    h_copy.hg_head[:, inds] .= nothing

    reached_vs, reached_es = forward_reachable(h_copy, source)

    # Path does not exist
    if target ∉ reached_vs
        return Set{Int}()
    end

    path = Set{Int}()

    # Try to minimize the size of the path by eliminating unnecessary dihyperedges
    for e in reached_es
        h_copy.hg_tail[:, e] .= nothing
        h_copy.hg_head[:, e] .= nothing

        # Only retain if dihyperedge is essential for reaching target
        if !is_reachable(h_copy, source, target, :vertex)
            # Restore dihyperedge to hypergraph copy
            h_copy[:, e] .= h[:, e]
            push!(path, e)
        end
    end

    return path
end

"""
    all_hyperpaths(h::H, source::Int, target::Int) where {H <: AbstractDirectedHypergraph}

    all_hyperpaths(
        h::DirectedHypergraph{T,V,E,D},
        source::Int,
        targets::Set{Int}
    ) where {T<:Real,V,E,D<:AbstractDict{Int,T}}

    all_hyperpaths(
        h::DirectedHypergraph{T,V,E,D},
        sources::Set{Int},
        target::Int
    ) where {T<:Real,V,E,D<:AbstractDict{Int,T}}

    all_hyperpaths(
        h::DirectedHypergraph{T,V,E,D},
        sources::Set{Int},
        targets::Set{Int}
    ) where {T<:Real,V,E,D<:AbstractDict{Int,T}}

Exhaustively (but efficiently) generate all hyperpaths in dihypergraph `h` from some source(s) to
some target(s), using the algorithm described by Krieger & Kececioglu (2022),
DOI: 10.1186/s13015-022-00217-9. 
    
Note that, ostensibly, this algorithm only works for single-source, single-sink pathfinding (i.e.,
with a single `source` index and a single `target` index). If the user provides multiple `sources`
and/or multiple `targets`, the multi-source/multi-sink problem will be reformulated as a
single-source, single-sink problem by adding a *metasource* vertex (connected to all source
vertices by a single dihyperedge) and/or *metatarget* vertex (connected to all target vertices by a
single dihyperedge).

**Arguments**

* `h` : Dihypergraph
* `source`/`sources` : Either the index of a single source vertex or a `Set{Int}` of vertex indices
* `target`/`targets` : Either the index of a single target vertex or a `Set{Int}` of vertex indices

**Returns**

`paths`, a `Set` of all possible hyperpaths from `source`/`sources` to `target`/`targets`, with
each path represented as a `Set{Int}` of dihyperedge indices

"""
function all_hyperpaths(h::H, source::Int, target::Int) where {H <: AbstractDirectedHypergraph}
    # Queue of subproblems
    # Subproblem is defined as a "out" set of dihyperedges (which must not be present in a path)
    # and "keep" dihyperedges which must be present in the path
    Q = Queue{Tuple{Set{Int}, Set{Int}}}()

    paths = Set{Set{Int}}()

    # Start with no restrictions
    enqueue!(Q, (Set{Int}(), Set{Int}()))

    while length(Q) > 0
        out, keep = dequeue!(Q)

        path = get_hyperpath(h, source, target, out)

        if length(path) > 0 && path ∉ paths
            push!(paths, path)

            k = keep
            for e in setdiff(path, keep)
                enqueue!(Q, (union(out, e), k))
                push!(k, e)
            end
        end
    end

    return paths
end

function all_hyperpaths(
        h::DirectedHypergraph{T, V, E, D},
        source::Int,
        targets::Set{Int}
    ) where {T <: Real, V, E, D <: AbstractDict{Int, T}}
    h_copy = deepcopy(h)

    # Add a single "metatarget" vertex to reformulate as single-source, single-sink pathfinding problem
    # The dihyperedge from the targets to the metatarget will have a cost of 0 associated with it
    metatarget = SimpleHypergraphs.add_vertex!(h_copy)
    meta_he = SimpleHypergraphs.add_hyperedge!(
        h_copy;
        vertices_tail = D(x => convert(T, 0) for x in targets),
        vertices_head = D(metatarget => convert(T, 0))
    )

    paths = all_hyperpaths(
        h_copy,
        source,
        metatarget,
    )

    # Remove the fictitious dihyperedge from the targets to the metatarget
    return Set(setdiff(p, Set{Int}(meta_he)) for p in paths)
end

function all_hyperpaths(
        h::DirectedHypergraph{T, V, E, D},
        sources::Set{Int},
        target::Int
    ) where {T <: Real, V, E, D <: AbstractDict{Int, T}}
    h_copy = deepcopy(h)

    # Add a single "metasource" vertex to reformulate as single-source, single-sink pathfinding problem
    # The dihyperedge from the metasource to the sources will have a cost of 0 associated with it
    metasource = SimpleHypergraphs.add_vertex!(h_copy)
    meta_he = SimpleHypergraphs.add_hyperedge!(
        h_copy;
        vertices_tail = D(metasource => convert(T, 0)),
        vertices_head = D(x => convert(T, 0) for x in sources)
    )

    paths = all_hyperpaths(
        h_copy,
        metasource,
        target,
    )

    # Remove the fictitious dihyperedge from the metasource to the sources
    return Set(setdiff(p, Set{Int}(meta_he)) for p in paths)
end

function all_hyperpaths(
        h::DirectedHypergraph{T, V, E, D},
        sources::Set{Int},
        targets::Set{Int}
    ) where {T <: Real, V, E, D <: AbstractDict{Int, T}}
    h_copy = deepcopy(h)

    # Add a single "metasource" vertex to reformulate as single-source, single-sink pathfinding problem
    # The dihyperedge from the metasource to the sources will have a cost of 0 associated with it
    metasource = SimpleHypergraphs.add_vertex!(h_copy)
    meta_he_source = SimpleHypergraphs.add_hyperedge!(
        h_copy;
        vertices_tail = D(metasource => convert(T, 0)),
        vertices_head = D(x => convert(T, 0) for x in sources)
    )

    # Add a single "metatarget" vertex to reformulate as single-source, single-sink pathfinding problem
    # The dihyperedge from the targets to the metatarget will have a cost of 0 associated with it
    metatarget = SimpleHypergraphs.add_vertex!(h_copy)
    meta_he_target = SimpleHypergraphs.add_hyperedge!(
        h_copy;
        vertices_tail = D(x => convert(T, 0) for x in targets),
        vertices_head = D(metatarget => convert(T, 0))
    )

    paths = all_hyperpaths(
        h_copy,
        metasource,
        metatarget,
    )

    # Remove fictitious dihyperedges
    return Set(setdiff(p, Set{Int}([meta_he_source, meta_he_target])) for p in paths)
end

"""
    initialize_ilp_model(
        h::H,
        source::Int,
        target::Int,
        dihyperedge_weights::AbstractVector{T}
    ) where {H<:AbstractDirectedHypergraph, T<:Real}

Define variables, objective, and initial constraints for integer linear programming-based
optimization of hyperpaths in dihypergraph `h` from vertex `source` to vertex `target`, where
dihyperedges have weights `dihyperedge_weights`.

See Krieger & Kececioglu (2022), DOI: 10.1186/s13015-022-00217-9 for details.

**Arguments**

* `h` : Dihypergraph
* `source` : Vertex index
* `target` : Vertex index
* `dihyperedge_weights` : `AbstractVector` of nonnegative real-valued weights with length equal to
    `nhe(h)`

**Returns**

* `model` : A `JuMP.Model`
* `x` : `JuMP` model variables
* `cuts` : `Vector{Set{Int}}` of `source`-`target` cuts containing groups of vertices not including
    `target`
* `crosses` : `Vector{BitVector}` of length `length(cuts)`. Each entry is a `BitVector` of length
    `nhe(h)` where each entry indicates if a dihyperedge crosses the corresponding cut

"""
function initialize_ilp_model(
        h::H,
        source::Int,
        target::Int,
        dihyperedge_weights::AbstractVector{T}
    ) where {H <: AbstractDirectedHypergraph, T <: Real}

    # First, verify that the problem is well-posed
    # That is, can `target` be reached from `source`
    @assert is_reachable(h, source, target, :vertex)

    # Initialize integer linear programming model
    model = JuMP.Model(GLPK.Optimizer)

    # Create one binary variable for each dihyperedge in `h`
    @JuMP.variable(model, x[1:nhe(h)], Bin)
    JuMP.set_start_value.(x, 1)

    # Verify that all variables are bound to be either 0 (not present in hyperpath) or 1 (present in hyperpath)
    @assert all(JuMP.is_binary.(x))

    # Define initial constraints
    for i in 1:nhe(h)
        # Tail-covering inequalities
        for v in keys(h.hg_tail.he2v[i])
            if v == source
                continue
            end

            in_hes = collect(keys(h.hg_head.v2he[v]))
            @JuMP.constraint(model, sum([x[ih] for ih in in_hes]) >= x[i])
        end

        # Head-hitting inequalities
        if target ∉ keys(h.hg_head.he2v[i])
            hits = Int[]
            for j in 1:nhe(h)
                if i != j && length(intersect(Set(keys(h.hg_head.he2v[i])), Set(keys(h.hg_tail.he2v[j])))) >= 1
                    push!(hits, j)
                end
            end

            @JuMP.constraint(model, sum([x[j] for j in hits]) >= x[i])
        end
    end

    # Target-production inequality
    @JuMP.constraint(model, sum([x[e] for e in keys(h.hg_head.v2he[target])]) >= 1)

    # Distance-based inequalities
    dist_ests = fill(typemax(T), nhv(h))
    for v in forward_reachable(h, source)[1]
        # TODO: this is inefficient
        # Currently, will repeat a lot of work
        # TODO: you are here; don't yet even understand what the problem is...
        dist_ests[v] = sum(
            [dihyperedge_weights[e] for e in shortest_hyperpath_kk_heuristic(h, source, v, dihyperedge_weights)]
        )
    end

    cuts = Set{Int}[]
    crosses = BitVector[]

    for d in unique(values(dist_ests))
        if d > dist_ests[target]
            continue
        end

        # Construct an s,t-cut, where the source is within the cut set of vertices and the target is not in the set
        cut_d = Set(v for (v, dist) in enumerate(dist_ests) if dist < d)
        push!(cut_d, source)
        setdiff!(cut_d, Set(target))

        push!(cuts, cut_d)

        cross = [
            issubset(Set(keys(h.hg_tail.he2v[i])), cut_d) && !issubset(Set(keys(h.hg_head.he2v[i])), cut_d)
                for i in 1:nhe(h)
        ]
        push!(crosses, BitVector(cross))

        @JuMP.constraint(model, dot(x, cross) >= 1)
    end

    # Define objective function
    @JuMP.objective(model, Min, dot(x, dihyperedge_weights))

    return model, x, cuts, crosses
end

"""
    shortest_hyperpath_kk_ilp(
        h::H,
        source::Int,
        target::Int,
        dihyperedge_weights::AbstractVector{T}
    ) where {H<:AbstractDirectedHypergraph, T<:Real}

    shortest_hyperpath_kk_ilp(
        h::DirectedHypergraph{T,V,E,D},
        source::Int,
        targets::Set{Int},
        dihyperedge_weights::AbstractVector{S}
    ) where {S<:Real,T<:Real,V,E,D<:AbstractDict{Int,T}}

    shortest_hyperpath_kk_ilp(
        h::DirectedHypergraph{T,V,E,D},
        sources::Set{Int},
        target::Int,
        dihyperedge_weights::AbstractVector{S}
    ) where {S<:Real,T<:Real,V,E,D<:AbstractDict{Int,T}}

    shortest_hyperpath_kk_ilp(
        h::DirectedHypergraph{T,V,E,D},
        sources::Set{Int},
        targets::Set{Int},
        dihyperedge_weights::AbstractVector{S}
    ) where {S<:Real,T<:Real,V,E,D<:AbstractDict{Int,T}}

Implements the exact dihypergraph pathfinding algorithm of Krieger & Kececioglu (2023),
DOI: 10.1089/cmb.2023.0242. This algorithm is guaranteed to find the optimal pathway from `source`
to `target` based on some nonnegative real-valued `dihyperedge_weights`, if and only if such a path
exists.

Note that, ostensibly, this algorithm only works for single-source, single-sink pathfinding (i.e.,
with a single `source` and a single `target`). If the user provides multiple `sources` and/or
multiple `targets`, the multi-source/multi-sink problem will be reformulated as a single-source,
single-sink problem by adding a *metasource* vertex (connected to all source vertices by a single,
0-cost dihyperedge) and/or *metatarget* vertex (connected to all target vertices by a single, 0-cost
dihyperedge).

**Arguments**

* `h` : Dihypergraph
* `source`/`sources` : Either the index of a single source vertex or a `Set{Int}` of vertex indices
* `target`/`targets` : Either the index of a single target vertex or a `Set{Int}` of vertex indices
* `dihyperedge_weights` : A `Vector` of nonnegative real-valued weights of length equal to `nhe(h)`

**Returns**

`path`, a `Set{Int}` of dihyperedge vertices constituting the shortest path from `source`/`sources`
to `target`/`targets`

"""
function shortest_hyperpath_kk_ilp(
        h::H,
        source::Int,
        target::Int,
        dihyperedge_weights::AbstractVector{T}
    ) where {H <: AbstractDirectedHypergraph, T <: Real}

    # TODO: do I need to carry `x` over like this? Not sure about variable scope
    model, x, cuts, crosses = initialize_ilp_model(h, source, target, dihyperedge_weights)

    JuMP.optimize!(model)

    # Convert floating-point solution into BitVector
    # TODO: Is this necessary w/ JuMP? Or will the output really be binary?
    solution = JuMP.value.(x) .> 0.5

    # Check for s,t-cuts that are not crossed by the current solution
    new_cuts, new_crosses = expand_cuts(h, source, target, cuts, crosses, solution)

    while length(new_cuts) > 0
        # Add new constraints to the model
        for cross in new_crosses
            @JuMP.constraint(model, dot(x, cross) >= 1)
        end

        # Re-optimize model with new cut-constraints
        JuMP.optimize!(model)
        solution = JuMP.value.(x) .> 0.5

        new_cuts, new_crosses = expand_cuts(h, source, target, new_cuts, new_crosses, solution)
    end

    return Set(findall(solution))
end

function shortest_hyperpath_kk_ilp(
        h::DirectedHypergraph{T, V, E, D},
        source::Int,
        targets::Set{Int},
        dihyperedge_weights::AbstractVector{S}
    ) where {S <: Real, T <: Real, V, E, D <: AbstractDict{Int, T}}
    h_copy = deepcopy(h)

    # Add a single "metatarget" vertex to reformulate as single-source, single-sink pathfinding problem
    # The dihyperedge from the targets to the metatarget will have a cost of 0 associated with it
    metatarget = SimpleHypergraphs.add_vertex!(h_copy)
    meta_he = SimpleHypergraphs.add_hyperedge!(
        h_copy;
        vertices_tail = D(x => convert(T, 0) for x in targets),
        vertices_head = D(metatarget => convert(T, 0))
    )

    path = shortest_hyperpath_kk_ilp(
        h_copy,
        source,
        metatarget,
        vcat(dihyperedge_weights, convert(S, 0))
    )

    # Remove the fictitious dihyperedge from the targets to the metatarget
    return setdiff(path, Set{Int}(meta_he))
end

function shortest_hyperpath_kk_ilp(
        h::DirectedHypergraph{T, V, E, D},
        sources::Set{Int},
        target::Int,
        dihyperedge_weights::AbstractVector{S}
    ) where {S <: Real, T <: Real, V, E, D <: AbstractDict{Int, T}}
    h_copy = deepcopy(h)

    # Add a single "metasource" vertex to reformulate as single-source, single-sink pathfinding problem
    # The dihyperedge from the metasource to the sources will have a cost of 0 associated with it
    metasource = SimpleHypergraphs.add_vertex!(h_copy)
    meta_he = SimpleHypergraphs.add_hyperedge!(
        h_copy;
        vertices_tail = D(metasource => convert(T, 0)),
        vertices_head = D(x => convert(T, 0) for x in sources)
    )

    path = shortest_hyperpath_kk_ilp(
        h_copy,
        metasource,
        target,
        vcat(dihyperedge_weights, convert(S, 0))
    )

    # Remove the fictitious dihyperedge from the metasource to the sources
    return setdiff(path, Set{Int}(meta_he))
end

function shortest_hyperpath_kk_ilp(
        h::DirectedHypergraph{T, V, E, D},
        sources::Set{Int},
        targets::Set{Int},
        dihyperedge_weights::AbstractVector{S}
    ) where {S <: Real, T <: Real, V, E, D <: AbstractDict{Int, T}}
    h_copy = deepcopy(h)

    # Add a single "metasource" vertex to reformulate as single-source, single-sink pathfinding problem
    # The dihyperedge from the metasource to the sources will have a cost of 0 associated with it
    metasource = SimpleHypergraphs.add_vertex!(h_copy)
    meta_he_source = SimpleHypergraphs.add_hyperedge!(
        h_copy;
        vertices_tail = D(metasource => convert(T, 0)),
        vertices_head = D(x => convert(T, 0) for x in sources)
    )

    # Add a single "metatarget" vertex to reformulate as single-source, single-sink pathfinding problem
    # The dihyperedge from the targets to the metatarget will have a cost of 0 associated with it
    metatarget = SimpleHypergraphs.add_vertex!(h_copy)
    meta_he_target = SimpleHypergraphs.add_hyperedge!(
        h_copy;
        vertices_tail = D(x => convert(T, 0) for x in targets),
        vertices_head = D(metatarget => convert(T, 0))
    )

    path = shortest_hyperpath_kk_ilp(
        h_copy,
        metasource,
        metatarget,
        vcat(dihyperedge_weights, [convert(T, 0), convert(T, 0)])
    )

    # Remove fictitious dihyperedges
    return setdiff(path, Set{Int}([meta_he_source, meta_he_target]))
end

"""
    expand_cuts(
        h::H,
        source::Int,
        target::Int,
        cuts::Vector{Set{Int}},
        crosses::Vector{BitVector},
        curr_sol::BitVector
    ) where {H<:AbstractDirectedHypergraph}

A helper function for `shortest_hyperpath_kk_ilp`.

Efficiently identify new `source`-`target` cuts of a dihypergraph `h` that are violated by current
solution of the integer linear programming problem `curr_sol`.

**Arguments**

* `h` : Dihypergraph
* `source` : Vertex index
* `target` : Vertex index
* `cuts` : `Vector` of `Set{Int}` of vertex indices constituting `source`-`target` cuts
* `crosses` : `Vector` of `BitVectors`, one for each cut in `cuts`. Each `BitVector` has length
    `nhe(h)` and indicates which dihyperedges cross a particular cut
* `curr_sol` : The current solution to the ILP problem, a `BitVector` of length `nhe(h)` indicating
    which dihyperedges are included in the solution

**Returns**

* `new_cuts` : A collection of new cuts, represented as in `cuts` (see above)
* `new_crosses` : `BitVectors` corresponding to each new cut in `new_cuts`, indicating which
    dihyperedges cross each cut

"""
function expand_cuts(
        h::H,
        source::Int,
        target::Int,
        cuts::Vector{Set{Int}},
        crosses::Vector{BitVector},
        curr_sol::BitVector
    ) where {H <: AbstractDirectedHypergraph}
    new_cuts = Set{Int}[]
    new_crosses = BitVector[]

    for (old_cut, old_cross) in zip(cuts, crosses)
        # Source augmentation
        new_cut = old_cut
        new_cross = old_cross
        for e in 1:nhe(h)
            # If `e` is an active dihyperedge that crosses the current cut, add the head of `e` to the cut so `e` no
            # longer crosses it
            if curr_sol[e] && new_cross[e]
                new_cut = union(new_cut, Set(keys(h.hg_head.he2v[e])))
                new_cross = [
                    issubset(Set(keys(h.hg_tail.he2v[i])), new_cut) && !issubset(Set(keys(h.hg_head.he2v[i])), new_cut)
                        for i in 1:nhe(h)
                ]
            end
        end

        # If this is a valid s,t-cut that no active dihyperedges cross, add it to the new cut list
        if !(any(new_cross) || target ∈ new_cut)
            push!(new_cuts, new_cut)
            push!(new_crosses, new_cross)
        end

        # Sink augmentation
        # TODO: you could (should) combine these two into one loop
        new_cut = old_cut
        new_cross = old_cross
        while any(curr_sol .&& new_cross)
            for e in 1:nhe(h)
                # If `e` is an active dihyperedge that crosses the current cut, remove vertices from the tail of `e` to the
                # cut so `e` no longer crosses it
                if curr_sol[e] && new_cross[e]
                    # Greedily pick the vertex in the tail of `e` that causes the fewest dihyperedges to newly cross this cut
                    new_cut_ev = Dict{Int, Set{Int}}()
                    new_cross_ev = Dict{Int, Vector{Bool}}()
                    greedy_v = 0
                    min_length = typemax(Int)
                    for v in keys(h.hg_tail.he2v[e])
                        new_cut_ev[v] = setdiff(new_cut, Set(v))
                        new_cross_ev[v] = [
                            issubset(Set(keys(h.hg_tail.he2v[i])), new_cut_ev[v]) &&
                                !issubset(Set(keys(h.hg_head.he2v[i])), new_cut_ev[v])
                                for i in 1:nhe(h)
                        ]
                        v_length = length(findall(map(!, old_cross) .&& new_cross_ev[v]))
                        if v_length < min_length
                            min_length = v_length
                            greedy_v = v
                        end
                    end

                    new_cut = new_cut_ev[greedy_v]
                    new_cross = new_cross_ev[greedy_v]
                end
            end
        end

        # If this is still a valid s,t-cut
        if !any(new_cross) && source ∈ new_cut
            push!(new_cuts, new_cut)
            push!(new_crosses, new_cross)
        end

    end

    return new_cuts, new_crosses
end
