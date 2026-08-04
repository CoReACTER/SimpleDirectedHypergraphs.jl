"""
    AbstractDirectedHypergraph{T} <: AbstractHypergraph{T}

An abstract directed hypergraph type storing information about vertices and dihyperedges.
"""
abstract type AbstractDirectedHypergraph{T} <: AbstractHypergraph{T} end

@traitimpl SimpleHypergraphs.IsDirected{AbstractDirectedHypergraph}
SimpleHypergraphs.isdirected(::Type{T}) where {T <: AbstractDirectedHypergraph} = true
SimpleHypergraphs.isdirected(X::T) where {T <: AbstractDirectedHypergraph} = true
