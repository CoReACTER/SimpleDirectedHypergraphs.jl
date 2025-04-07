"""
  Return the number of edges in 2-section view `t` of a directed hypergraph.
"""
function Graphs.ne(t::TwoSectionView{H}) where {H<:AbstractDirectedHypergraph}
    s = 0
    for x in 1:nhe(t.h)
        s += length(t.h.hg_tail.he2v[x]) * length(t.h.hg_head.he2v[x])
    end
    s
end

"""
    Graphs.all_neighbors(
        t::TwoSectionView,
        v::Integer;
        incoming::Bool = false, outgoing::Bool = true
    ) where {H<:AbstractDirectedHypergraph}

Returns N(v) (the vertex v is not included in N(v))

If incoming is true (default true), include nodes connected to v
by directed hyperedges where v is in the head.

If outgoing is true (default true), include nodes connected to v
by directed hyperedges where v is in the tail.

If both incoming and outgoing are false, returns an empty set
"""
function Graphs.all_neighbors(
    t::TwoSectionView{H},
    v::Integer;
    incoming::Bool = true,
    outgoing::Bool = true
) where {H<:AbstractDirectedHypergraph}
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
    collect(neighbors) #returns the corresponding array
end


function Graphs.has_edge(t::TwoSectionView{H}, s, d) where {H<:AbstractDirectedHypergraph}
    s == d && return false
    !isempty(intersect(keys(t.h.hg_tail.v2he[s]), keys(t.h.hg_head.v2he[d])))
end


Graphs.outneighbors(t::TwoSectionView{H}, v::Integer) where {H<:AbstractDirectedHypergraph} =
    Graphs.all_neighbors(t, v, incoming=false)

Graphs.inneighbors(t::TwoSectionView{H}, v::Integer) where {H<:AbstractDirectedHypergraph} =
    Graphs.all_neighbors(t, v, outgoing=false)


"""
    Graphs.SimpleGraph(t::TwoSectionView{H}) where {H<:AbstractDirectedHypergraph}

Creates a `Graphs.SimpleGraph` representation of a `TwoSectionView` t.

This creates a copy of the date. Note that the weights information is not stored
in the created `SimpleGraph`.
"""
function Graphs.SimpleDiGraph(t::TwoSectionView{H}) where {H<:AbstractDirectedHypergraph}
    g = SimpleDiGraph(nv(t))
    for v in Graphs.vertices(t)
        outneighbors_v = Graphs.outneighbors(t, v)

        for neighbor in outneighbors_v
            add_edge!(g, v, neighbor)
        end
    end
    g
end


Graphs.is_directed(t::TwoSectionView{H}) where {H<:AbstractDirectedHypergraph} = true
Graphs.is_directed(::Type{TwoSectionView{H}}) where {H<:AbstractDirectedHypergraph} = true


"""
    shortest_path(t::TwoSectionView{H}, source::Int, target::Int) where {H<:AbstractDirectedHypergraph}

Finds a single shortest path in a graph `b` between vertices
`source` and `target`.
Note that if several paths of the same length exist, only one
will be returned.

"""
function SimpleHypergraphs.shortest_path(t::TwoSectionView{H}, source::Int, target::Int) where {H<:AbstractDirectedHypergraph}
    checkbounds(t.h.hg_tail.v2he, source)
    checkbounds(t.h.hg_head.v2he, target)
    dj = dijkstra_shortest_paths(t, source)
    enumerate_paths(dj)[target]
end


"""
    Graphs.SimpleGraphs.fadj(t::TwoSectionView{H}) where {H<:AbstractDirectedHypergraph}

Generates an adjency list for this view of a directed hypergraph.
"""
function Graphs.SimpleGraphs.fadj(t::TwoSectionView{H}) where {H<:AbstractDirectedHypergraph}
    res = [Vector{Int}() for _ in 1:Graphs.nv(t)]
    
    for he in 1:nhe(t.h)
        vs_tail, vs_head = getvertices(t.h, he)
        for v_tail in keys(vs_tail)
            for v_head in keys(vs_head)
                if v_head != v_tail
                    append!(res[v_tail], v_head)
                end
            end
        end
    end

    sort!.(res)
end


"""
    Graphs.SimpleGraphs.badj(t::TwoSectionView{H}) where {H<:AbstractDirectedHypergraph}

Generates an adjency list for this view of a hypergraph.
"""
function Graphs.SimpleGraphs.badj(t::TwoSectionView{H}) where {H<:AbstractDirectedHypergraph}
    res = [Vector{Int}() for _ in 1:Graphs.nv(t)]
    for he in 1:nhe(t.h)
        vs_tail, vs_head = getvertices(t.h, he)
        for v_tail in keys(vs_tail)
            for v_head in keys(vs_head)
                if v_head != v_tail
                    append!(res[v_head], v_tail)
                end
            end
        end
    end
    sort!.(res)
end


Graphs.SimpleGraphs.fadj(t::TwoSectionView{H}, v::Integer) where {H<:AbstractDirectedHypergraph} = Graphs.outneighbors(t,v)
Graphs.SimpleGraphs.badj(t::TwoSectionView{H}, v::Integer) where {H<:AbstractDirectedHypergraph} = Graphs.inneighbors(t,v)


"""
    get_twosection_adjacency_mx(h::H; count_self_loops::Bool=false,
                                replace_weights::Union{Nothing,Real}=nothing) where {H<:AbstractDirectedHypergraph}

Returns an adjacency matrix for a two section view of a hypergraph `h`.
"""
function SimpleHypergraphs.get_twosection_adjacency_mx(
    h::H;
    count_self_loops::Bool=false,
    replace_weights::Union{Nothing,Real}=nothing
    ) where {T<:Real, H<:AbstractDirectedHypergraph{Tuple{Union{T, Nothing}, Union{T, Nothing}}}}
    mx = zeros(replace_weights === nothing ? Tuple{T,T} : typeof(replace_weights), nhv(h), nhv(h))
    for he in 1:nhe(h)
        for vt in keys(h.hg_tail.he2v[he])
            for vh in keys(h.hg_head.he2v[he])
                vt == vh && !count_self_loops && continue
                mx[vt,vh] += replace_weights === nothing ? (h.hg_tail.he2v[he][vt], h.hg_head.he2v[he][vh]) : replace_weights
            end
        end
    end
    mx
end