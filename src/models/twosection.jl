"""
    Graphs.ne(t::TwoSectionView{H}) where {H <: AbstractDirectedHypergraph}

Return the number of edges in 2-section view `t` of a directed hypergraph.

**Arguments**

* `t` : `TwoSectionView` of a dihypergraph

**Returns**

`s`, the number of edges in `t`

"""
function Graphs.ne(t::TwoSectionView{H}) where {H <: AbstractDirectedHypergraph}
    s = 0
    for x in 1:nhe(t.h)
        s += length(t.h.hg_tail.he2v[x]) * length(t.h.hg_head.he2v[x])
    end
    return s
end

"""
    Graphs.all_neighbors(
        t::TwoSectionView,
        v::Integer;
        incoming::Bool = true,
	outgoing::Bool = true
    ) where {H<:AbstractDirectedHypergraph}

Returns the indices of all vertices forward- and/or backward-adjacent to the vertex with index `v`
in two-section view `t` of a dihypergraph.

**Arguments**

* `t` : `TwoSectionView` of a dihypergraph
* `v` : Vertex index
* `incoming` : If `true` (default `true`), include nodes connected to v by dihyperedges where `v` is in the head.
* `outgoing` : If `true` (default `true`), include nodes connected to v by dihyperedges where `v` is in the tail.

**Returns**

A `Vector` of vertex indices. If both `incoming` and `outgoing` are `false`, returns an empty
`Vector` 

"""
function Graphs.all_neighbors(
        t::TwoSectionView{H},
        v::Integer;
        incoming::Bool = true,
        outgoing::Bool = true
    ) where {H <: AbstractDirectedHypergraph}
    neighbors = Set{Int}()

    if !(incoming || outgoing)
        return collect(neighbors)
    end

    if incoming
        for he in keys(t.h.hg_head.v2he[v])
            union!(neighbors, keys(t.h.hg_tail.he2v[he]))
        end
    end

    if outgoing
        for he in keys(t.h.hg_tail.v2he[v])
            union!(neighbors, keys(t.h.hg_head.he2v[he]))
        end
    end

    delete!(neighbors, v) #remove v from its neighborhood
    return collect(neighbors) #returns the corresponding array
end

"""
    Graphs.has_edge(t::TwoSectionView{H}, s::Int, d::Int) where {H <: AbstractDirectedHypergraph}

**Arguments**

* `t` : `TwoSectionView` of a dihypergraph
* `s` : Source vertex index
* `d` : Destination vertex index

**Returns**

`true` if there is an edge in the two-section view from `s` to `d`, and `false` otherwise

"""
function Graphs.has_edge(t::TwoSectionView{H}, s::Int, d::Int) where {H <: AbstractDirectedHypergraph}
    s == d && return false
    return !isempty(intersect(keys(t.h.hg_tail.v2he[s]), keys(t.h.hg_head.v2he[d])))
end

"""
    Graphs.outneighbors(t::TwoSectionView{H}, v::Int)

**Arguments**

* `t` : `TwoSectionView` of a dihypergraph
* `v` : Vertex index

**Returns**

Returns all outgoing neighbors of vertex `v` in the two-section representation of a dihypergraph;
that is, all vertices `w` in `t` where there exists an edge with `v` as the source and `w` as the
destination

"""
Graphs.outneighbors(t::TwoSectionView{H}, v::Int) where {H <: AbstractDirectedHypergraph} =
    Graphs.all_neighbors(t, v, incoming = false)

"""
    Graphs.inneighbors(t::TwoSectionView{H}, v::Int)

**Arguments**

* `t` : `TwoSectionView` of a dihypergraph
* `v` : Vertex index

**Returns**

Returns all ingoing neighbors of vertex `v` in the two-section representation of a dihypergraph;
that is, all vertices `w` in `b` where there exists an edge with `w` as the source and `v` as the
destination

"""
Graphs.inneighbors(t::TwoSectionView{H}, v::Int) where {H <: AbstractDirectedHypergraph} =
    Graphs.all_neighbors(t, v, outgoing = false)


"""
    Graphs.SimpleDiGraph(t::TwoSectionView{H}) where {H<:AbstractDirectedHypergraph}

Creates a `Graphs.SimpleDiGraph` representation of a `TwoSectionView` t.

This creates a copy of the date. Note that the weights information is not stored
in the created `SimpleDiGraph`.

**Arguments**

* `t` : `TwoSectionView` of a dihypergraph

**Returns**

A `SimpleDiGraph`

"""
function Graphs.SimpleDiGraph(t::TwoSectionView{H}) where {H <: AbstractDirectedHypergraph}
    g = SimpleDiGraph(nv(t))
    for v in Graphs.vertices(t)
        outneighbors_v = Graphs.outneighbors(t, v)

        for neighbor in outneighbors_v
            add_edge!(g, v, neighbor)
        end
    end
    return g
end

"""
    Graphs.is_directed(t::TwoSectionView{H})
    Graphs.is_directed(::Type{TwoSectionView{H}})

**Arguments**

* `t` : `TwoSectionView` of a dihypergraph 

**Returns**

`true`

"""
Graphs.is_directed(t::TwoSectionView{H}) where {H <: AbstractDirectedHypergraph} = true
Graphs.is_directed(::Type{TwoSectionView{H}}) where {H <: AbstractDirectedHypergraph} = true


"""
    SimpleHypergraphs.shortest_path(
	t::TwoSectionView{H},
	source::Int,
	target::Int
    ) where {H<:AbstractDirectedHypergraph}

Finds a single shortest path in a two-section representation of a dihypergraph `t` between vertices
`source` and `target`. Note that if several paths of the same length exist, only one will be
returned.

**Arguments**

* `t` : `TwoSectionView` of a dihypergraph
* `source` : Vertex index
* `target` : Vertex index

**Returns**

A path, represented by a sequence of edge indices

"""
function SimpleHypergraphs.shortest_path(t::TwoSectionView{H}, source::Int, target::Int) where {H <: AbstractDirectedHypergraph}
    checkbounds(t.h.hg_tail.v2he, source)
    checkbounds(t.h.hg_head.v2he, target)
    dj = dijkstra_shortest_paths(t, source)
    return enumerate_paths(dj)[target]
end


"""
    Graphs.SimpleGraphs.fadj(t::TwoSectionView{H}) where {H<:AbstractDirectedHypergraph}

Generates a forward adjacency list for the two-section representation of a dihypergraph.

**Arguments**

* `t` : `TwoSectionView` of a dihypergraph

**Returns**

For each vertex `v` in `t`, a `Vector{Int}` with vertex indices reachable via one edge with `v` as
the source

"""
function Graphs.SimpleGraphs.fadj(t::TwoSectionView{H}) where {H <: AbstractDirectedHypergraph}
    res = [Vector{Int}() for _ in 1:Graphs.nv(t)]

    for he in 1:nhe(t.h)
        vs_tail, vs_head = getvertices(t.h, he)
        for v_tail in keys(vs_tail)
            for v_head in keys(vs_head)
                if v_head != v_tail
                    push!(res[v_tail], v_head)
                end
            end
        end
    end

    return sort!.(res)
end


"""
    Graphs.SimpleGraphs.badj(t::TwoSectionView{H}) where {H<:AbstractDirectedHypergraph}

Generates an adjency list for this view of a hypergraph.
"""
function Graphs.SimpleGraphs.badj(t::TwoSectionView{H}) where {H <: AbstractDirectedHypergraph}
    res = [Vector{Int}() for _ in 1:Graphs.nv(t)]
    for he in 1:nhe(t.h)
        vs_tail, vs_head = getvertices(t.h, he)
        for v_tail in keys(vs_tail)
            for v_head in keys(vs_head)
                if v_head != v_tail
                    push!(res[v_head], v_tail)
                end
            end
        end
    end
    return sort!.(res)
end

"""
    Graphs.SimpleGraphs.fadj(t::TwoSectionView{H}, v::Integer) where {H <: AbstractDirectedHypergraph}

Generates a forward adjacency list for one vertex in the two-section representation of a dihypergraph.

**Arguments**

* `t` : `TwoSectionView` of a dihypergraph
* `v` : Vertex index

**Returns**

A `Vector{Int}` with vertex indices reachable via one edge with `v` as the source

"""
Graphs.SimpleGraphs.fadj(t::TwoSectionView{H}, v::Integer) where {H <: AbstractDirectedHypergraph} = Graphs.outneighbors(t, v)

"""
    Graphs.SimpleGraphs.badj(t::TwoSectionView{H}, v::Integer) where {H <: AbstractDirectedHypergraph}

Generates a backward adjacency list for one vertex in the two-section representation of a dihypergraph.

**Arguments**

* `t` : `TwoSectionView` of a dihypergraph
* `v` : Vertex index

**Returns**

A `Vector{Int}` with vertex indices reachable via one edge with `v` as the destination

"""
Graphs.SimpleGraphs.badj(t::TwoSectionView{H}, v::Integer) where {H <: AbstractDirectedHypergraph} = Graphs.inneighbors(t, v)


"""
    get_twosection_adjacency_mx(
	h::H;
	count_self_loops::Bool,
        replace_weights::Union{Nothing,Real}
    ) where {H<:AbstractDirectedHypergraph}

Returns an adjacency matrix for a two section view of a hypergraph `h`.

Note: If two vertices are connected by more than one dihyperedge in `h`, the associated weights of
the vertices in each relevant dihyperedge will be summed.

**Arguments**

* `h` : Dihypergraph
* `count_self_loops` : If `true` (default `false`), include edges from a vertex `v` to `v`
* `replace_weights` : If not `nothing` (default `nothing`), replace all weights with the provided
    value.

**Returns**

A `m`-by-`m` adjacency matrix, where `m = nhv(h)`

"""
function SimpleHypergraphs.get_twosection_adjacency_mx(
        h::H;
        count_self_loops::Bool = false,
        replace_weights::Union{Nothing, Real} = nothing
    ) where {T <: Real, H <: AbstractDirectedHypergraph{Tuple{Union{T, Nothing}, Union{T, Nothing}}}}
    mx = zeros(replace_weights === nothing ? Vector{T} : typeof(replace_weights), nhv(h), nhv(h))
    for he in 1:nhe(h)
        for vt in keys(h.hg_tail.he2v[he])
            for vh in keys(h.hg_head.he2v[he])
                vt == vh && !count_self_loops && continue
                mx[vt, vh] += replace_weights === nothing ? [h.hg_tail.he2v[he][vt], h.hg_head.he2v[he][vh]] : replace_weights
            end
        end
    end
    return mx
end
