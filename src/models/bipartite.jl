"""
    Graphs.ne(b::BipartiteView{H}) where {H <: AbstractDirectedHypergraph}

**Arguments**

* `b` : `BipartiteView` of a dihypergraph

**Returns**

Returns the number of edges of the bipartite representation of a dihypergraph

"""
Graphs.ne(b::BipartiteView{H}) where {H <: AbstractDirectedHypergraph} = sum(length.(b.h.hg_tail.v2he)) + sum(length.(b.h.hg_head.v2he))

"""
    Graphs.all_neighbors(b::BipartiteView{H}, v::Int) where {H <: AbstractDirectedHypergraph}

**Arguments**

* `b` : `BipartiteView` of a dihypergraph
* `v` : Vertex index

**Returns**

Returns all unique neighbors of vertex `v` in the bipartite representation of a dihypergraph

"""
function Graphs.all_neighbors(b::BipartiteView{H}, v::Int) where {H <: AbstractDirectedHypergraph}
    n1 = nhv(b.h)

    return if v <= n1
        t, h = gethyperedges(b.h, v)
        n1 .+ unique([collect(keys(t)); collect(keys(h))])
    else
        t, h = getvertices(b.h, v - n1)
        unique([collect(keys(t)); collect(keys(h))])
    end
end

"""
    Graphs.has_edge(b::BipartiteView{H}, s::Int, d::Int) where {H <: AbstractDirectedHypergraph}

**Arguments**

* `b` : `BipartiteView` of a dihypergraph
* `s` : Source vertex index
* `d` : Destination vertex index

**Returns**

`Bool`; `true` if, in the original dihypergraph (`h`), `s` is the index of a vertex in the tail of
the dihyperedge with index `d - nhv(h)` or if `d` is the index of a vertex in the head of the
dihyperedge with index `s - nhv(h)`, and `false` otherwise

"""

function Graphs.has_edge(b::BipartiteView{H}, s, d) where {H <: AbstractDirectedHypergraph}
    n1 = nhv(b.h)

    return if s <= n1
        d > n1 && haskey(b.h.hg_tail.v2he[s], d - n1)
    else
        d <= n1 && haskey(b.h.hg_head.he2v[s - n1], d)
    end
end

"""
    Graphs.outneighbors(b::BipartiteView{H}, v::Int) where {H <: AbstractDirectedHypergraph}

**Arguments**

* `b` : `BipartiteView` of a dihypergraph
* `v` : Vertex index

**Returns**

Returns all outgoing neighbors of vertex `v` in the bipartite representation of a dihypergraph;
that is, all vertices `w` in `b` where there exists an edge with `v` as the source and `w` as the
destination

"""
function Graphs.outneighbors(
        b::BipartiteView{H},
        v::Integer
    ) where {H <: AbstractDirectedHypergraph}

    n1 = nhv(b.h)

    return if v <= n1
        t, _ = gethyperedges(b.h, v)
        n1 .+ collect(keys(t))
    else
        _, h = getvertices(b.h, v - n1)
        collect(keys(h))
    end

end

"""
    Graphs.inneighbors(b::BipartiteView{H}, v::Int) where {H <: AbstractDirectedHypergraph}

**Arguments**

* `b` : `BipartiteView` of a dihypergraph
* `v` : Vertex index

**Returns**

Returns all ingoing neighbors of vertex `v` in the bipartite representation of a dihypergraph;
that is, all vertices `w` in `b` where there exists an edge with `w` as the source and `v` as the
destination

"""
function Graphs.inneighbors(
        b::BipartiteView{H},
        v::Integer
    ) where {H <: AbstractDirectedHypergraph}

    n1 = nhv(b.h)

    return if v <= n1
        _, h = gethyperedges(b.h, v)
        n1 .+ collect(keys(h))
    else
        t, _ = getvertices(b.h, v - n1)
        collect(keys(t))
    end

end


"""
    Graphs.SimpleGraph(b::BipartiteView{H}) where {H<:AbstractDirectedHypergraph}

Creates a `Graphs.SimpleGraph` representation of a `BipartiteView` b.

This creates a copy of the data. Note that the weights information is not stored
in the created `SimpleGraph`.

**Arguments**

* `b` : `BipartiteView` of a dihypergraph

**Returns**

`g`, a `SimpleGraph`

"""
function Graphs.SimpleGraph(b::BipartiteView{H}) where {H <: AbstractDirectedHypergraph}
    g = SimpleGraph(nv(b))

    n1 = nhv(b.h)
    for v in 1:n1
        t, h = gethyperedges(b.h, v)

        for he in unique([collect(keys(t)); collect(keys(h))])
            add_edge!(g, v, n1 + he)
        end
    end
    return g
end


"""
    Graphs.SimpleDiGraph(b::BipartiteView{H}) where {H<:AbstractDirectedHypergraph}

Creates a `Graphs.SimpleDiGraph` representation of a `BipartiteView` b.

This creates a copy of the data. Note that the weights information is not stored
in the created `SimpleDiGraph`.

**Arguments**

* `b` : `BipartiteView` of a dihypergraph

**Returns**

`g`, a `SimpleDiGraph`

"""
function Graphs.SimpleDiGraph(b::BipartiteView{H}) where {H <: AbstractDirectedHypergraph}
    g = SimpleDiGraph(nv(b))

    n1 = nhv(b.h)

    for v in 1:n1
        t, h = gethyperedges(b.h, v)

        for he in keys(t)
            add_edge!(g, v, n1 + he)
        end

        for he in keys(h)
            add_edge!(g, n1 + he, v)
        end

    end

    return g
end

"""
    Graphs.isdirected(b::BipartiteView{H}) where {H <: AbstractDirectedHypergraph}

    Graphs.isdirected(::Type{BipartiteView{H}}) where {H <: AbstractDirectedHypergraph}

**Arguments**

* `b` : `BipartiteView` of a dihypergraph

**Returns**

`true`

"""
Graphs.is_directed(b::BipartiteView{H}) where {H <: AbstractDirectedHypergraph} = true

Graphs.is_directed(::Type{BipartiteView{H}}) where {H <: AbstractDirectedHypergraph} = true


"""
    SimpleHypergraphs.shortest_path(
	b::BipartiteView{H},
	source::Int,
	target::Int
    ) where {H<:AbstractDirectedHypergraph}

Finds a single shortest path in a graph `b` between vertices `source` and `target`. Note that if
several paths of the same length exist, only one will be returned.

**Arguments**

* `b` : `BipartiteView` of a dihypergraph
* `source` : Vertex (in `BipartiteView`) index
* `target` : Vertex (in `BipartiteView`) index

**Returns**

A path, represented by a sequence of edge indices

"""
function SimpleHypergraphs.shortest_path(b::BipartiteView{H}, source::Int, target::Int) where {H <: AbstractDirectedHypergraph}
    checkbounds(b.h.hg_tail.v2he, source)
    checkbounds(b.h.hg_tail.v2he, target)
    checkbounds(b.h.hg_head.v2he, source)
    checkbounds(b.h.hg_head.v2he, target)

    dj = dijkstra_shortest_paths(b, source)
    return enumerate_paths(dj)[target][1:2:end]
end


"""
    Graphs.SimpleGraphs.fadj(b::BipartiteView{H}) where {H<:AbstractDirectedHypergraph}

Generates a forward adjacency list for the bipartite representation of a dihypergraph.

**Arguments**

* `b` : `BipartiteView` of a dihypergraph

**Returns**

For each vertex `v` in `b`, a `Vector{Int}` with vertex indices reachable via one edge with `v` as
the source

"""
function Graphs.SimpleGraphs.fadj(b::BipartiteView{H}) where {H <: AbstractDirectedHypergraph}
    res = Vector{Vector{Int}}(undef, Graphs.nv(b))

    h_nv = length(b.h.hg_tail.v2he)
    for i in 1:h_nv
        res[i] = h_nv .+ sort!(collect(keys(b.h.hg_tail.v2he[i])))
    end
    for i in 1:length(b.h.hg_head.he2v)
        res[i + h_nv] = sort!(collect(keys(b.h.hg_head.he2v[i])))
    end
    return res
end

"""
    Graphs.SimpleGraphs.badj(b::BipartiteView{H}) where {H<:AbstractDirectedHypergraph}

Generates a backward adjacency list for this view of a directed hypergraph.

**Arguments**

* `b` : `BipartiteView` of a dihypergraph

**Returns**

For each vertex `v` in `b`, a `Vector{Int}` with vertex indices reachable via one edge with `v` as
the destination

"""
function Graphs.SimpleGraphs.badj(b::BipartiteView{H}) where {H <: AbstractDirectedHypergraph}
    res = Vector{Vector{Int}}(undef, Graphs.nv(b))

    h_nv = length(b.h.hg_head.v2he)
    for i in 1:h_nv
        res[i] = h_nv .+ sort!(collect(keys(b.h.hg_head.v2he[i])))
    end
    for i in 1:length(b.h.hg_tail.he2v)
        res[i + h_nv] = sort!(collect(keys(b.h.hg_tail.he2v[i])))
    end
    return res
end

"""
    Graphs.SimpleGraphs.fadj(b::BipartiteView{H}, v::Int) where {H<:AbstractDirectedHypergraph}

Generates a forward adjacency list for one vertex in the bipartite representation of a dihypergraph.

**Arguments**

* `b` : `BipartiteView` of a dihypergraph
* `v` : Vertex index

**Returns**

A `Vector{Int}` with vertex indices reachable via one edge with `v` as the source

"""
Graphs.SimpleGraphs.fadj(b::BipartiteView{H}, v::Int) where {H <: AbstractDirectedHypergraph} = Graphs.outneighbors(b, v)

"""
    Graphs.SimpleGraphs.badj(b::BipartiteView{H}, v::Int) where {H<:AbstractDirectedHypergraph}

Generates a backward adjacency list for one vertex in the bipartite representation of a dihypergraph.

**Arguments**

* `b` : `BipartiteView` of a dihypergraph
* `v` : Vertex index

**Returns**

A `Vector{Int}` with vertex indices reachable via one edge with `v` as the destination

"""
Graphs.SimpleGraphs.badj(b::BipartiteView{H}, v::Integer) where {H <: AbstractDirectedHypergraph} = Graphs.inneighbors(b, v)

# TODO: has_cycles
