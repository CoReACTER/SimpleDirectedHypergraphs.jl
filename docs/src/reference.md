Reference
=========

```@meta
CurrentModule = SimpleDirectedHypergraphs
DocTestSetup = quote
    using SimpleHypergraphs
    using SimpleDirectedHypergraphs
end
```

Abstract types
---------------------

```@docs

AbstractDirectedHypergraph{T} <: AbstractSimpleHypergraph{Tuple{Union{T, Nothing}, Union{T, Nothing}}}
```

Creating a directed hypergraph
---------------------

```@docs

DirectedHypergraph
BasicDirectedHypergraph
SimpleHypergraphs.random_model(::Int, ::Int, ::Type{H}; ::Bool) where {H <: AbstractDirectedHypergraph}
SimpleHypergraphs.random_kuniform_model(::Int, ::Int, ::Int, ::Type{H}) where {H <: AbstractDirectedHypergraph}
SimpleHypergraphs.random_dregular_model(::Int, ::Int, ::Int, ::Type{H}) where {H <: AbstractDirectedHypergraph}
```

Manipulating vertices and hyperedges
------------------------------------
```@docs

SimpleHypergraphs.add_hyperedge!(::DirectedHypergraph{T, V, E, D}; ::D, ::D, ::Union{E,Nothing}, ::Union{E,Nothing}) where {T <: Real, V, E, D <: AbstractDict{Int,T}}
SimpleHypergraphs.add_hyperedge!(::BasicDirectedHypergraph{T, D}; ::D, ::D) where {T <: Real, D <: AbstractDict{Int,T}}

SimpleHypergraphs.add_vertex!(::DirectedHypergraph{T, V, E, D}; ::D, ::D, ::Union{V,Nothing}) where {T <: Real, V, E, D <: AbstractDict{Int,T}}
SimpleHypergraphs.add_vertex!(::BasicDirectedHypergraph{T, D}; ::D, ::D) where {T <: Real, D <: AbstractDict{Int,T}}

SimpleHypergraphs.set_vertex_meta!(DirectedHypergraph{T, V, E, D}, ::Union{V,Nothing}, ::Int) where {T <: Real, V, E, D <: AbstractDict{Int,T}}
SimpleHypergraphs.get_vertex_meta(DirectedHypergraph{T, V, E, D}, ::Int) where {T <: Real, V, E, D <: AbstractDict{Int,T}}
SimpleHypergraphs.set_hyperedge_meta!(::DirectedHypergraph{T, V, E, D}, ::Union{E,Nothing}, ::Union{E,Nothing}, ::Int) where {T <: Real, V, E, D <: AbstractDict{Int,T}}
SimpleHypergraphs.get_hyperedge_meta(::DirectedHypergraph{T, V, E, D}, ::Int) where {T <: Real, V, E, D <: AbstractDict{Int,T}}

SimpleHypergraphs.remove_vertex!(::DirectedHypergraph, ::Int)
SimpleHypergraphs.remove_vertex!(::BasicDirectedHypergraph, ::Int)

SimpleHypergraphs.remove_hyperedge!(::DirectedHypergraph, ::Int)
SimpleHypergraphs.remove_hyperedge!(::BasicDirectedHypergraph, ::Int)
```

Hypergraph array getters and setters
------------------------------------

For a directed hypergraph, an additional index can refer to either the tail (1) or head (2) of a directed hyperedge:
```jldoctest
dh = BasicDirectedHypergraph{Int64}(2,3);
dh[1,1,1] = 1;
dh[2,2,1] = 2;
dh

# output

2×3 BasicDirectedHypergraph{Int64, Dict{Int64, Int64}}:
 (1, nothing)  (nothing, nothing)  (nothing, nothing)
 (nothing, 2)  (nothing, nothing)  (nothing, nothing)
```

Setting with slices is not currently possible, and users must currently use slices on the component undirected hypergraphs instead. However, accessing using slices is allowed:
```jldoctest
dh = BasicDirectedHypergraph{Int64}(4,2);
dh.hg_tail[1:2,1] .= 1;
dh.hg_head[3:4,1] .= 2;
dh[:,1]

# output

4-element Vector{Tuple{Union{Nothing, Int64}, Union{Nothing, Int64}}}:
 (1, nothing)
 (1, nothing)
 (nothing, 2)
 (nothing, 2)
```

```@docs
Base.getindex(::H, ::Vararg{Int,2}) where {H <: AbstractDirectedHypergraph}

Base.setindex!(::H, ::Nothing, ::Vararg{Int,2}) where {H <: AbstractDirectedHypergraph}
Base.setindex!(::H, ::Real, ::Vararg{Int,2}) where {H <: AbstractDirectedHypergraph}
Base.setindex!(::H, ::Tuple{Union{Real, Nothing}, Union{Real, Nothing}}, ::Vararg{Int,2}) where {H <: AbstractDirectedHypergraph}
Base.setindex!(::H, ::Nothing, ::Vararg{Int,3}) where {H <: AbstractDirectedHypergraph}
Base.setindex!(::H, ::Real, ::Vararg{Int,3}) where {H <: AbstractDirectedHypergraph}
```

Hypergraph info
---------------
```@docs
SimpleHypergraphs.nhv(::H) where {H <: AbstractDirectedHypergraph}
SimpleHypergraphs.nhe(::H) where {H <: AbstractDirectedHypergraph}

SimpleHypergraphs.getvertices(::H, ::Int) where {H <: AbstractDirectedHypergraph}
SimpleHypergraphs.gethyperedges(::H, ::Int) where {H <: AbstractDirectedHypergraph}

get_weakly_connected_components(::H) where {H <: AbstractDirectedHypergraph}
get_strongly_connected_components(::H) where {H <: AbstractDirectedHypergraph}

SimpleHypergraphs.get_twosection_adjacency_mx(::H; ::Bool, ::Union{Nothing,Real}) where {H<:AbstractDirectedHypergraph}
SimpleHypergraphs.random_walk(::H, ::Int; ::Function, ::Function, ::Bool) where {H <: AbstractDirectedHypergraph}

SimpleHypergraphs.dual(h::DirectedHypergraph)
SimpleHypergraphs.dual(h::BasicDirectedHypergraph)
```

I/O
---

Directed hypergraphs can be saved as and loaded from JSON- and EHGF-formatted files, where the EHGF format is a close derivative of HGF.

```@docs
SimpleHypergraphs.hg_save
SimpleHypergraphs.hg_load
```
