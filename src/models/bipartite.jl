"""
  Return the number of edges in a bipartite view `b` of a hypergraph.
"""
Graphs.ne(b::BipartiteView{H}) where {H<:AbstractDirectedHypergraph} = sum(length.(b.h.hg_tail.v2he)) + sum(length.(b.h.hg_head.v2he))


function Graphs.all_neighbors(b::BipartiteView{H}, v::Integer) where {H<:AbstractDirectedHypergraph}
    n1 = nhv(b.h)

    if v <= n1
      t, h = gethyperedges(b.h, v)
      n1 .+ unique([collect(keys(t)); collect(keys(h))])
    else
      t, h = getvertices(b.h, v - n1)
      unique([collect(keys(t)); collect(keys(h))])
    end
end


function Graphs.has_edge(b::BipartiteView{H}, s, d) where {H<:AbstractDirectedHypergraph}
    n1 = nhv(b.h)

    if s <= n1
        d > n1 && haskey(b.h.hg_tail.v2he[s], d - n1)
    else
        d <= n1 && haskey(b.h.hg_head.he2v[s-n1], d)
    end
end


function Graphs.outneighbors(
    b::BipartiteView{H},
    v::Integer
    ) where {H<:AbstractDirectedHypergraph}

    n1 = nhv(b.h)

    if v <= n1
      t, _ = gethyperedges(b.h, v)
      n1 .+ collect(keys(t))
    else
      _, h = getvertices(b.h, v - n1)
      collect(keys(h))
    end

end

function Graphs.inneighbors(
    b::BipartiteView{H},
    v::Integer
    ) where {H<:AbstractDirectedHypergraph}

    n1 = nhv(b.h)

    if v <= n1
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
"""
function Graphs.SimpleGraph(b::BipartiteView{H}) where {H<:AbstractDirectedHypergraph}
    g = SimpleGraph(nv(b))

    n1 = nhv(b.h)
    for v in 1:n1
        t, h = gethyperedges(b.h, v)

        for he in unique([collect(keys(t)); collect(keys(h))])
            add_edge!(g, v, n1 + he)
        end
    end
    g
end


"""
    Graphs.SimpleDiGraph(b::BipartiteView{H}) where {H<:AbstractDirectedHypergraph}

Creates a `Graphs.SimpleDiGraph` representation of a `BipartiteView` b.

This creates a copy of the data. Note that the weights information is not stored
in the created `SimpleDiGraph`.
"""
function Graphs.SimpleDiGraph(b::BipartiteView{H}) where {H<:AbstractDirectedHypergraph}
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

    g
end


Graphs.is_directed(b::BipartiteView{H}) where {H<:AbstractDirectedHypergraph} = true

Graphs.is_directed(::Type{BipartiteView{H}}) where {H<:AbstractDirectedHypergraph} = true


"""
    shortest_path(b::BipartiteView{H}, source::Int, target::Int) where {H<:AbstractDirectedHypergraph}

Finds a single shortest path in a graph `b` between vertices
`source` and `target`.
Note that if several paths of the same length exist, only one
will be returned.

"""
function SimpleHypergraphs.shortest_path(b::BipartiteView{H}, source::Int, target::Int) where {H<:AbstractDirectedHypergraph}
    checkbounds(b.h.hg_tail.v2he, source)
    checkbounds(b.h.hg_tail.v2he, target)
    checkbounds(b.h.hg_head.v2he, source)
    checkbounds(b.h.hg_head.v2he, target)

    dj = dijkstra_shortest_paths(b, source)
    enumerate_paths(dj)[target][1:2:end]
end


"""
    Graphs.SimpleGraphs.fadj(b::BipartiteView{H}) where {H<:AbstractDirectedHypergraph}

Generates a forward adjacency list for this view of a directed hypergraph.
"""
function Graphs.SimpleGraphs.fadj(b::BipartiteView{H}) where {H<:AbstractDirectedHypergraph}
    res = Vector{Vector{Int}}(undef, Graphs.nv(b))

    h_nv = length(b.h.hg_tail.v2he)
    for i in 1:h_nv
       res[i] = h_nv .+ sort!(collect(keys(b.h.hg_tail.v2he[i])))
    end
    for i in 1:length(b.h.hg_head.he2v)
        res[i+h_nv] = sort!(collect(keys(b.h.hg_head.he2v[i])))
    end
    res
end

"""
    Graphs.SimpleGraphs.badj(b::BipartiteView{H}) where {H<:AbstractDirectedHypergraph}

Generates a backward adjacency list for this view of a directed hypergraph.
"""
function Graphs.SimpleGraphs.badj(b::BipartiteView{H}) where {H<:AbstractDirectedHypergraph}
    res = Vector{Vector{Int}}(undef, Graphs.nv(b))

    h_nv = length(b.h.hg_head.v2he)
    for i in 1:h_nv
       res[i] = h_nv .+ sort!(collect(keys(b.h.hg_head.v2he[i])))
    end
    for i in 1:length(b.h.hg_tail.he2v)
        res[i+h_nv] = sort!(collect(keys(b.h.hg_tail.he2v[i])))
    end
    res
end

Graphs.SimpleGraphs.fadj(b::BipartiteView{H}, v::Integer) where {H<:AbstractDirectedHypergraph} = Graphs.outneighbors(b,v)
Graphs.SimpleGraphs.badj(b::BipartiteView{H}, v::Integer) where {H<:AbstractDirectedHypergraph} = Graphs.inneighbors(b,v)

# TODO: has_cycles
