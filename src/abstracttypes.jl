"""
    AbstractDirectedHypergraph{T} <: AbstractHypergraph{T}

An abstract directed hypergraph type storing information about vertices and hyperedges.
"""
abstract type AbstractDirectedHypergraph{T} <: AbstractHypergraph{Tuple{Union{T, Nothing}, Union{T, Nothing}}} end

@traitimpl SimpleHypergraphs.IsDirected{AbstractDirectedHypergraph}
SimpleHypergraphs.isdirected(::Type{T}) where {T<:AbstractDirectedHypergraph} = true
SimpleHypergraphs.isdirected(X::T) where {T<:AbstractDirectedHypergraph} = true