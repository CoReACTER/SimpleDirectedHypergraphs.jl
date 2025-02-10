"""
    DirectedHypergraph{T} <: AbstractDirectedHypergraph{T}

A directed hypergraph storing information about vertices and hyperedges.

This implementation is based on guidance from Przemysław Szufel;
    see https://github.com/pszufe/SimpleHypergraphs.jl/issues/45
This allows us to manipulate DirectedHypergraphs using Hypergraph functionality
There is danger of a user manipulating individual `hg_tail` and `hg_head` (undirected) hypergraphs
Is there a smart way to prevent this?
TODO: reconsider this design choice

**Constructors**

    DirectedHypergraph{T}(n::Integer,k::Integer) where {T<:Real}
    DirectedHypergraph{T,V}(n::Integer, k::Integer;
        v_meta=Vector{Union{V,Nothing}}(nothing, n)
        ) where {T<:Real, V}
    DirectedHypergraph{T,E}(n::Integer, k::Integer;
        he_meta_tail=Vector{Union{E,Nothing}}(nothing, k),
        he_meta_head=Vector{Union{E,Nothing}}(nothing, k)
        ) where {T<:Real, E}
    DirectedHypergraph{T,V,E}(n::Integer, k::Integer;
        v_meta=Vector{Union{V,Nothing}}(nothing, n),
        he_meta_tail=Vector{Union{E,Nothing}}(nothing, k),
        he_meta_head=Vector{Union{E,Nothing}}(nothing, k)
        ) where {T<:Real, V, E}
    DirectedHypergraph{T,V,E,D}(n::Integer, k::Integer,
        v_meta=Vector{Union{V,Nothing}}(nothing, n),
        he_meta_tail=Vector{Union{E,Nothing}}(nothing, k),
        he_meta_head=Vector{Union{E,Nothing}}(nothing, k)
        ) where {T<:Real,V,E,D<:AbstractDict{Int,T}}

Construct a hypergraph with a given number of vertices and hyperedges.
Optionally, values of type `V` can be stored at vertices and values of type `E`
can be stored at hyperedges. By default the hypergraph uses a `Dict{Int,T}` for
the internal data storage, however a different dictionary such as `SortedDict`
to ensure result replicability can be used (e.g. when doing stochastic
simulations on hypergraphs).

    DirectedHypergraph(
        m_tail::AbstractMatrix{Union{T, Nothing}},
        m_head::AbstractMatrix{Union{T, Nothing}}
    ) where {T<:Real}    
    DirectedHypergraph{T}(
        m_tail::AbstractMatrix{Union{T, Nothing}},
        m_head::AbstractMatrix{Union{T, Nothing}}
    ) where {T<:Real}
    DirectedHypergraph{T,V}(
        m_tail::AbstractMatrix{Union{T, Nothing}},
        m_head::AbstractMatrix{Union{T, Nothing}};
        v_meta::Vector{Union{Nothing,V}}=Vector{Union{Nothing,V}}(nothing, size(m,1)),
    ) where {T<:Real,V}
    DirectedHypergraph{T,E}(
        m_tail::AbstractMatrix{Union{T, Nothing}},
        m_head::AbstractMatrix{Union{T, Nothing}};
        he_meta_tail::Vector{Union{Nothing,E}}=Vector{Union{Nothing,E}}(nothing, size(m,2)),
        he_meta_head::Vector{Union{Nothing,E}}=Vector{Union{Nothing,E}}(nothing, size(m,2))
    ) where {T<:Real,E}
    DirectedHypergraph{T,V,E}(
        m_tail::AbstractMatrix{Union{T, Nothing}},
        m_head::AbstractMatrix{Union{T, Nothing}};
        v_meta::Vector{Union{Nothing,V}}=Vector{Union{Nothing,V}}(nothing, size(m,1)),
        he_meta_tail::Vector{Union{Nothing,E}}=Vector{Union{Nothing,E}}(nothing, size(m,2)),
        he_meta_head::Vector{Union{Nothing,E}}=Vector{Union{Nothing,E}}(nothing, size(m,2))
    ) where {T<:Real,V,E}
    DirectedHypergraph{T,V,E,D}(
        m_tail::AbstractMatrix{Union{T, Nothing}},
        m_head::AbstractMatrix{Union{T, Nothing}};
        v_meta::Vector{Union{Nothing,V}}=Vector{Union{Nothing,V}}(nothing, size(m,1)),
        he_meta_tail::Vector{Union{Nothing,E}}=Vector{Union{Nothing,E}}(nothing, size(m,2)),
        he_meta_head::Vector{Union{Nothing,E}}=Vector{Union{Nothing,E}}(nothing, size(m,2))
    ) where {T<:Real,V,E,D<:AbstractDict{Int,T}}

Construct a directed hypergraph using its matrix representation.
In the matrix representation rows are vertices and columns are hyperedges.
Optionally, values of type `V` can be stored at vertices and values of type `E`
can be stored at hyperedges. By default the hypergraph uses a `Dict{Int,T}` for
the internal data storage, however a different dictionary such as `SortedDict`
to ensure result replicability can be used (e.g. when doing stochastic
simulations on hypergraphs).

    DirectedHypergraph(g::Graphs.DiGraph)

Constructs a directed hypergraph of degree 2 by making a deep copy of a
Graphs.DiGraph. A `SortedDict` will be used for internal data storage of the
hypergraph.

    DirectedHypergraph{T,V,D}(
        hg_tail::BasicHypergraph{T,D},
        hg_head::BasicHypergraph{T,D};
        v_meta::Vector{Union{Nothing,V}}=Vector{Union{Nothing,V}}(nothing, size(m,1)),
    ) where {T<:Real,V,D<:AbstractDict{Int, T}}
    DirectedHypergraph{T,E,D}(
        hg_tail::BasicHypergraph{T,D},
        hg_head::BasicHypergraph{T,D};
        he_meta_tail::Vector{Union{Nothing,E}}=Vector{Union{Nothing,E}}(nothing, size(m,2)),
        he_meta_head::Vector{Union{Nothing,E}}=Vector{Union{Nothing,E}}(nothing, size(m,2))
    ) where {T<:Real,E,D<:AbstractDict{Int, T}}
    DirectedHypergraph{T,V,E,D}(
        hg_tail::BasicHypergraph{T,D},
        hg_head::BasicHypergraph{T,D};
        v_meta::Vector{Union{Nothing,V}}=Vector{Union{Nothing,V}}(nothing, size(m,1)),
        he_meta_tail::Vector{Union{Nothing,E}}=Vector{Union{Nothing,E}}(nothing, size(m,2)),
        he_meta_head::Vector{Union{Nothing,E}}=Vector{Union{Nothing,E}}(nothing, size(m,2))
    ) where {T<:Real,V,E,D<:AbstractDict{Int, T}}

Constructs a directed hypergraph from two undirected basic hypergraphs, one with hyperedges
containing "tail" vertices and one with hyperedges containing "head"
verticies.

    DirectedHypergraph{T,V,E,D}(
        hg_tail::Hypergraph{T,V,E,D},
        hg_head::Hypergraph{T,V,E,D}
    ) where {T<:Real,V,E,D<:AbstractDict{Int, T}}

Constructs a directed hypergraph from two hypergraphs potentially containing metadata. Throws
an error if the vertex metadata of the two hypergraphs is not element-for-element identical.

**Arguments**

* `T` : type of weight values stored in the hypergraph's adjacency matrix
* `V` : type of values stored in the vertices of the hypergraph
* `E` : type of values stored in the edges of the hypergraph
* `D` : dictionary for storing values the default is `Dict{Int, T}`
* `n` : number of vertices
* `k` : number of hyperedges
* `m` : a matrix representation rows are vertices and columns are hyperedges
* `g` : a (directed) graph representation of the hypergraph
* `hg_tail`: an undirected hypergraph representing the tail half of
    the directed hypergraph
* `hg_head`: an undirected hypergraph representing the head half of
    the directed hypergraph
"""

struct DirectedHypergraph{T<:Real,V,E,D<:AbstractDict{Int, T}} <: AbstractDirectedHypergraph{T}
    hg_tail::BasicHypergraph{T,D}
    hg_head::BasicHypergraph{T,D}

    v_meta::Vector{Union{V,Nothing}}
    he_meta_tail::Vector{Union{E,Nothing}}
    he_meta_head::Vector{Union{E,Nothing}}

    DirectedHypergraph{T,V,E,D}(
        n::Integer, k::Integer,
        v_meta=Vector{Union{V, Nothing}}(nothing, n),
        he_meta_tail=Vector{Union{E, Nothing}}(nothing, k),
        he_meta_head=Vector{Union{E, Nothing}}(nothing, k)
        ) where {T<:Real,V,E,D<:AbstractDict{Int, T}} = 
        new{T,V,E,D}(
            BasicHypergraph{T,D}(n, k),
            BasicHypergraph{T,D}(n, k),
            v_meta, he_meta_tail, he_meta_head
        )

    function DirectedHypergraph{T,V,E,D}(
        hg_tail::BasicHypergraph{T,D},
        hg_head::BasicHypergraph{T,D};
        v_meta::Vector{Union{Nothing,V}}=Vector{Union{Nothing,V}}(nothing, size(hg_tail,1)),
        he_meta_tail::Vector{Union{Nothing,E}}=Vector{Union{Nothing,E}}(nothing, size(hg_tail,2)),
        he_meta_head::Vector{Union{Nothing,E}}=Vector{Union{Nothing,E}}(nothing, size(hg_tail,2))
        ) where {T<:Real,V,E,D<:AbstractDict{Int, T}}
        @assert size(hg_tail) == size(hg_head)

        @assert length(v_meta) == size(hg_tail,1)
        @assert length(he_meta_tail) == size(hg_tail,2)
        @assert length(he_meta_head) == size(hg_head,2)

        new{T,V,E,D}(
            hg_tail,
            hg_head,
            v_meta,
            he_meta_tail,
            he_meta_head
        )
    end
end

DirectedHypergraph{T,V,E}(n::Integer, k::Integer) where {T<:Real, V, E} = DirectedHypergraph{T,V,E,Dict{Int,T}}(n, k)

DirectedHypergraph{T,V}(n::Integer, k::Integer) where {T<:Real, V} = DirectedHypergraph{T,V,Nothing,Dict{Int,T}}(n, k)

DirectedHypergraph{T,D}(n::Integer, k::Integer) where {T<:Real, D<:AbstractDict{Int, T}} = DirectedHypergraph{T,Nothing,Nothing,D}(n, k)

DirectedHypergraph{T}(n::Integer, k::Integer) where {T<:Real} = DirectedHypergraph{T,Nothing,Nothing,Dict{Int,T}}(n, k)

DirectedHypergraph(n::Integer, k::Integer) = DirectedHypergraph{Bool,Nothing,Nothing,Dict{Int,Bool}}(n, k)


function DirectedHypergraph{T,D}(
    hg_tail::BasicHypergraph{T,D},
    hg_head::BasicHypergraph{T,D}
    ) where {T<:Real,D<:AbstractDict{Int, T}}

    DirectedHypergraph{T,Nothing,Nothing,D}(
        hg_tail,
        hg_head
    )
end

function DirectedHypergraph{T,V,D}(
    hg_tail::BasicHypergraph{T,D},
    hg_head::BasicHypergraph{T,D}
    ) where {T<:Real,V,D<:AbstractDict{Int, T}}

    DirectedHypergraph{T,V,Nothing,D}(
        hg_tail,
        hg_head
    )
end

function DirectedHypergraph{T,V,E,D}(
    hg_tail::Hypergraph{T,V,E,D},
    hg_head::Hypergraph{T,V,E,D}
) where {T<:Real,V,E,D<:AbstractDict{Int, T}}
    @assert size(hg_tail) == size(hg_head)

    n, k = size(hg_tail)
    shg_tail = BasicHypergraph(n, k)
    shg_head = BasicHypergraph(n, k)

    # TODO: test behavior on this
    shg_tail .= hg_tail
    shg_head .= hg_head

    if all(hg_tail.v_meta .== hg_head.v_meta)
        DirectedHypergraph{T,V,E,D}(shg_tail, shg_head; v_meta=hg_tail.v_meta, he_meta_tail=hg_tail.he_meta, he_meta_head=hg_head.he_meta)
    else
        @warn "Vertex metadata for tail and head hypergraphs not identical; discarding vertex metadata."
        DirectedHypergraph{T,V,E,D}(shg_tail, shg_head; he_meta_tail=hg_tail.he_meta, he_meta_head=hg_head.he_meta)
    end
end


function DirectedHypergraph{T,V,E,D}(
        m_tail::AbstractMatrix{Union{T, Nothing}},
        m_head::AbstractMatrix{Union{T, Nothing}};
        v_meta::Vector{Union{Nothing,V}}=Vector{Union{Nothing,V}}(nothing, size(m_tail,1)),
        he_meta_tail::Vector{Union{Nothing,E}}=Vector{Union{Nothing,E}}(nothing, size(m_tail,2)),
        he_meta_head::Vector{Union{Nothing,E}}=Vector{Union{Nothing,E}}(nothing, size(m_tail,2))
    ) where {T<:Real,V,E,D<:AbstractDict{Int,T}}

    # TODO: ensure that sizes of matrices identical?

    # Arbitrary, since sizes are identical
    n, k = size(m_tail)

    hg_tail = BasicHypergraph{T,D}(n, k)
    hg_tail .= m_tail
    
    hg_head = BasicHypergraph{T,D}(n, k)
    hg_head .= m_head

    DirectedHypergraph{T,V,E,D}(
        hg_tail, hg_head;
        v_meta=v_meta,
        he_meta_tail=he_meta_tail,
        he_meta_head=he_meta_head
    )
end

function DirectedHypergraph{T,V,E}(
    m_tail::AbstractMatrix{Union{T, Nothing}},
    m_head::AbstractMatrix{Union{T, Nothing}};
    v_meta::Vector{Union{Nothing,V}}=Vector{Union{Nothing,V}}(nothing, size(m_tail,1)),
    he_meta_tail::Vector{Union{Nothing,E}}=Vector{Union{Nothing,E}}(nothing, size(m_tail,2)),
    he_meta_head::Vector{Union{Nothing,E}}=Vector{Union{Nothing,E}}(nothing, size(m_tail,2))
) where {T<:Real,V,E}

    # Arbitrary, since sizes are identical
    n, k = size(m_tail)

    hg_tail = BasicHypergraph{T,Dict{Int,T}}(n, k)
    hg_tail .= m_tail

    hg_head = BasicHypergraph{T,Dict{Int,T}}(n, k)
    hg_head .= m_head

    DirectedHypergraph{T,V,E,Dict{Int,T}}(
        hg_tail, hg_head;
        v_meta=v_meta,
        he_meta_tail=he_meta_tail,
        he_meta_head=he_meta_head
    )
end

function DirectedHypergraph{T,V}(
    m_tail::AbstractMatrix{Union{T, Nothing}},
    m_head::AbstractMatrix{Union{T, Nothing}};
    v_meta::Vector{Union{Nothing,V}}=Vector{Union{Nothing,V}}(nothing, size(m_tail,1)),
) where {T<:Real,V}

    # Arbitrary, since sizes are identical
    n, k = size(m_tail)

    hg_tail = BasicHypergraph{T,Dict{Int,T}}(n, k)
    hg_tail .= m_tail

    hg_head = BasicHypergraph{T,Dict{Int,T}}(n, k)
    hg_head .= m_head

    DirectedHypergraph{T,V,Nothing,Dict{Int,T}}(
        hg_tail,
        hg_head;
        v_meta=v_meta
    )
end

function DirectedHypergraph{T}(
    m_tail::AbstractMatrix{Union{T, Nothing}},
    m_head::AbstractMatrix{Union{T, Nothing}}
) where {T<:Real}

    # Arbitrary, since sizes are identical
    n, k = size(m_tail)

    hg_tail = BasicHypergraph{T,Dict{Int,T}}(n, k)
    hg_tail .= m_tail

    hg_head = BasicHypergraph{T,Dict{Int,T}}(n, k)
    hg_head .= m_head

    DirectedHypergraph{T,Nothing,Nothing,Dict{Int,T}}(hg_tail, hg_head)
end

function DirectedHypergraph(
    m_tail::AbstractMatrix{Union{T, Nothing}},
    m_head::AbstractMatrix{Union{T, Nothing}}
) where {T<:Real}

    # Arbitrary, since sizes are identical
    n, k = size(m_tail)

    hg_tail = BasicHypergraph{T,Dict{Int,T}}(n, k)
    hg_tail .= m_tail

    hg_head = BasicHypergraph{T,Dict{Int,T}}(n, k)
    hg_head .= m_head

    DirectedHypergraph{T,Nothing,Nothing,Dict{Int,T}}(
        hg_tail,
        hg_head
    )
end


function DirectedHypergraph(g::Graphs.DiGraph)
    h = DirectedHypergraph{Bool,Nothing,Nothing,SortedDict{Int,Bool}}(maximum(vertices(g)), ne(g))
    e = 0
    for edge in edges(g)
        e+=1
        h[1,edge.src,e] = true
        h[2,edge.dst,e] = true
    end
    h
end


"""
    BasicDirectedHypergraph{T} <: AbstractDirectedHypergraph{T}

A directed hypergraph storing vertices and hyperedges, with no additional
information other than weight.

This implementation is based on guidance from Przemysław Szufel;
    see https://github.com/pszufe/SimpleHypergraphs.jl/issues/45
This allows us to manipulate DirectedHypergraphs using Hypergraph functionality
There is danger of a user manipulating individual `hg_tail` and `hg_head` (undirected) hypergraphs
Is there a smart way to prevent this?
TODO: reconsider this design choice

**Constructors**

    BasicDirectedHypergraph{T}(n::Integer,k::Integer) where {T<:Real}
    BasicDirectedHypergraph{T,D}(n::Integer, k::Integer) where {T<:Real,D<:AbstractDict{Int,T}}

Construct a hypergraph with a given number of vertices and hyperedges.
By default the hypergraph uses a `Dict{Int,T}` for
the internal data storage, however a different dictionary such as `SortedDict`
to ensure result replicability can be used (e.g. when doing stochastic
simulations on hypergraphs).

    BasicDirectedHypergraph(
        m_tail::AbstractMatrix{Union{T, Nothing}},
        m_head::AbstractMatrix{Union{T, Nothing}}
    ) where {T<:Real}
    BasicDirectedHypergraph{T}(
        m_tail::AbstractMatrix{Union{T, Nothing}},
        m_head::AbstractMatrix{Union{T, Nothing}}
    ) where {T<:Real}
    BasicDirectedHypergraph{T,D}(
        m_tail::AbstractMatrix{Union{T, Nothing}},
        m_head::AbstractMatrix{Union{T, Nothing}}
    ) where {T<:Real,D<:AbstractDict{Int,T}}

Construct a directed hypergraph using its matrix representation.
In the matrix representation rows are vertices and columns are hyperedges.
By default the hypergraph uses a `Dict{Int,T}` for the internal data storage,
however a different dictionary such as `SortedDict` to ensure result
replicability can be used (e.g. when doing stochastic simulations on hypergraphs).

    DirectedHypergraph(g::Graphs.DiGraph)

Constructs a directed hypergraph of degree 2 by making a deep copy of a
Graphs.DiGraph. A `SortedDict` will be used for internal data storage of the
hypergraph.

    DirectedHypergraph{T,D}(
        hg_tail::BasicHypergraph{T,D},
        hg_head::BasicHypergraph{T,D}
    ) where {T<:Real,D<:AbstractDict{Int, T}}

Constructs a directed hypergraph from two undirected basic hypergraphs, one with hyperedges
containing "tail" vertices and one with hyperedges containing "head"
verticies.

**Arguments**

* `T` : type of weight values stored in the hypergraph's adjacency matrix
* `V` : type of values stored in the vertices of the hypergraph
* `E` : type of values stored in the edges of the hypergraph
* `D` : dictionary for storing values the default is `Dict{Int, T}`
* `n` : number of vertices
* `k` : number of hyperedges
* `m` : a matrix representation rows are vertices and columns are hyperedges
* `g` : a (directed) graph representation of the hypergraph
* `hg_tail`: an undirected hypergraph representing the tail half of
    the directed hypergraph
* `hg_head`: an undirected hypergraph representing the head half of
    the directed hypergraph
"""
struct BasicDirectedHypergraph{T<:Real,D<:AbstractDict{Int, T}} <: AbstractDirectedHypergraph{T}
    hg_tail::BasicHypergraph{T,D}
    hg_head::BasicHypergraph{T,D}

    BasicDirectedHypergraph{T,D}(
        n::Integer, k::Integer,
        ) where {T<:Real,D<:AbstractDict{Int, T}} =
        new{T,D}(
            BasicHypergraph{T,D}(n, k),
            BasicHypergraph{T,D}(n, k)
        )

    function BasicDirectedHypergraph{T,D}(
        hg_tail::BasicHypergraph{T,D},
        hg_head::BasicHypergraph{T,D}
        ) where {T<:Real,D<:AbstractDict{Int, T}}
        
        @assert size(hg_tail) == size(hg_head)

        new{T,D}(hg_tail, hg_head)
    end
end


BasicDirectedHypergraph{T}(n::Integer, k::Integer) where {T<:Real} = BasicDirectedHypergraph{T,Dict{Int,T}}(n, k)

BasicDirectedHypergraph(n::Integer, k::Integer) = BasicDirectedHypergraph{Bool,Dict{Int,Bool}}(n, k)


function BasicDirectedHypergraph{T, D}(
    m_tail::AbstractMatrix{Union{T, Nothing}},
    m_head::AbstractMatrix{Union{T, Nothing}}
) where {T<:Real, D<:AbstractDict{Int, T}}

    # Arbitrary, since sizes are identical
    n, k = size(m_tail)

    hg_tail = BasicHypergraph{T,D}(n, k)
    hg_tail .= m_tail

    hg_head = BasicHypergraph{T,D}(n, k)
    hg_head .= m_head

    BasicDirectedHypergraph{T,D}(hg_tail, hg_head)
end


function BasicDirectedHypergraph{T}(
    m_tail::AbstractMatrix{Union{T, Nothing}},
    m_head::AbstractMatrix{Union{T, Nothing}}
) where {T<:Real}

    # Arbitrary, since sizes are identical
    n, k = size(m_tail)

    hg_tail = BasicHypergraph{T,Dict{Int,T}}(n, k)
    hg_tail .= m_tail

    hg_head = BasicHypergraph{T,Dict{Int,T}}(n, k)
    hg_head .= m_head

    BasicDirectedHypergraph{T,Dict{Int,T}}(hg_tail, hg_head)
end

function BasicDirectedHypergraph(
    m_tail::AbstractMatrix{Union{T, Nothing}},
    m_head::AbstractMatrix{Union{T, Nothing}}
) where {T<:Real}

    # Arbitrary, since sizes are identical
    n, k = size(m_tail)

    hg_tail = BasicHypergraph{T,Dict{Int,T}}(n, k)
    hg_tail .= m_tail

    hg_head = BasicHypergraph{T,Dict{Int,T}}(n, k)
    hg_head .= m_head

    BasicDirectedHypergraph{T,Dict{Int,T}}(
        hg_tail,
        hg_head
    )
end


function BasicDirectedHypergraph(g::Graphs.DiGraph)
    h = BasicDirectedHypergraph{Bool,SortedDict{Int,Bool}}(maximum(vertices(g)), ne(g))
    e = 0
    for edge in edges(g)
        e+=1
        h[1,edge.src,e] = true
        h[2,edge.dst,e] = true
    end
    h
end


# implementing traits on types
@traitimpl SimpleHypergraphs.HasMeta{DirectedHypergraph}
SimpleHypergraphs.hasmeta(::Type{DirectedHypergraph}) = true


# TODO: this is awkward...
const DIRECTED_HYPERGRAPH_VALID_FIRST_INDICES = [1,2]

# TODO: can this entirely replace the above? Index setting seems problematic...
@enum HyperedgeDirection begin
    tail = 1
    head = 2
end


# AbstractArray interface functions

"""
    Base.getindex(h::H, idx::Vararg{Int,2}) where {H <: AbstractDirectedHypergraph}

Returns a value for a given vertex-hyperedge pair `idx` for a directed hypergraph `h`.
If a vertex does not belong to a hyperedge `nothing` is returned.

"""
@inline function Base.getindex(h::H, idx::Vararg{Int,2}) where {H <: AbstractDirectedHypergraph}
    @boundscheck checkbounds(h.hg_tail, idx...)
    @boundscheck checkbounds(h.hg_head, idx...)

    tail_value = get(h.hg_tail.v2he[idx[1]], idx[2], nothing)
    head_value = get(h.hg_head.v2he[idx[1]], idx[2], nothing)

    (tail_value, head_value)
end

"""
    Base.setindex!(h::H, ::Nothing, idx::Vararg{Int,2}) where {H <: AbstractDirectedHypergraph}

Removes a vertex from a given hyperedge for a directed hypergraph `h` and a given vertex-hyperedge pair `idx`.
Note that trying to remove a vertex from a hyperedge when it is not present will not throw an error.

"""
@inline function Base.setindex!(h::H, ::Nothing, idx::Vararg{Int,2}) where {H <: AbstractDirectedHypergraph}
    @boundscheck checkbounds(h.hg_tail, idx...)
    @boundscheck checkbounds(h.hg_head, idx...)
    setindex!(h.hg_tail, nothing, idx...)
    setindex!(h.hg_head, nothing, idx...)
    h
end


"""
    Base.setindex!(h::H, v::Real, idx::Vararg{Int,2}) where {H <: AbstractDirectedHypergraph}

Adds a vertex to a hyperedge (represented by indices `idx`) and assigns value
`v` to be stored with that assignment.

"""
@inline function Base.setindex!(h::H, v::Real, idx::Vararg{Int,2}) where {H <: AbstractDirectedHypergraph}
    @boundscheck checkbounds(h.hg_tail, idx...)
    @boundscheck checkbounds(h.hg_head, idx...)

    setindex!(h.hg_tail, v, idx...)
    setindex!(h.hg_head, v, idx...)
    h
end


"""
    Base.setindex!(h::H, v::Tuple{Union{Real, Nothing}, Union{Real, Nothing}}, idx::Vararg{Int,2}) where {H <: AbstractDirectedHypergraph}

Manipulates a hyperedge (represented by indices `idx`), either adding a vertex to the 
ingoing and/or head sides of the hyperedge and assigning a value associated with that assignment,
or else removing a vertex from the ingoing/head sides of the hyperedge.

Here, `v` is a 2-tuple where the first element is the value that will be assigned to the ingoing part of the hyperedge
and the second element is the value that will be assigned to the head part. A value of `nothing` means that the
vertex will be removed from that side of the hyperedge.

"""
@inline function Base.setindex!(h::H, v::Tuple{Union{Real, Nothing}, Union{Real, Nothing}}, idx::Vararg{Int,2}) where {H <: AbstractDirectedHypergraph}
    @boundscheck checkbounds(h.hg_tail, idx...)
    @boundscheck checkbounds(h.hg_head, idx...)
    
    setindex!(h.hg_tail, v[1], idx...)
    setindex!(h.hg_head, v[2], idx...)

    h
end


"""
    Base.setindex!(h::H, ::Nothing, idx::Vararg{Int,3}) where {H <: AbstractDirectedHypergraph}

Removes a vertex from a given hyperedge for a directed hypergraph `h` and a given side-vertex-hyperedge pair `idx`.
If the first index of `idx` is 1, then the vertex will be removed from the tail hyperedge; if `idx` is 2, then
the vertex will be removed from the head hyperedge. 
Note that trying to remove a vertex from a hyperedge when it is not present will not throw an error.

"""
@inline function Base.setindex!(h::H, ::Nothing, idx::Vararg{Int,3}) where {H <: AbstractDirectedHypergraph}
    @boundscheck checkbounds(DIRECTED_HYPERGRAPH_VALID_FIRST_INDICES, idx[1])

    if idx[1] == 1
        side = h.hg_tail
    else
        side = h.hg_head
    end
    
    @boundscheck checkbounds(side, idx[2:end]...)
    
    setindex!(side, nothing, idx[2], idx[3])

    h
end


"""
    Base.setindex!(h::H, v::Real, idx::Vararg{Int,3}) where {H <: AbstractDirectedHypergraph}

Adds a vertex to a hyperedge (represented by indices `idx`, where the first index must be either
1 - referring to an tail hyperedge - or 2 - referring to an head hyperedge) and assigns value
`v` to be stored with that assignment.

"""
@inline function Base.setindex!(h::H, v::Real, idx::Vararg{Int,3}) where {H <: AbstractDirectedHypergraph}
    @boundscheck checkbounds(DIRECTED_HYPERGRAPH_VALID_FIRST_INDICES, idx[1])

    if idx[1] == 1
        side = h.hg_tail
    else
        side = h.hg_head
    end
    
    @boundscheck checkbounds(side, idx[2:end]...)
    
    setindex!(side, v, idx[2]..., idx[3]...)

    h
end


"""
    getvertices(h::H, he_id::Int) where {H <: AbstractDirectedHypergraph}

Returns vertices from a directed hypergraph `a` for a given hyperedge `he_id`.
Vertex indices are given in a tuple `(in, out)`, where `in` are the tail vertices
and `out` are the head vertices

"""
@inline SimpleHypergraphs.getvertices(h::H, he_id::Int) where {H <: AbstractDirectedHypergraph} = (h.hg_tail.he2v[he_id], h.hg_head.he2v[he_id])


"""
    gethyperedges(h::H, v_id::Int) where {H <: AbstractDirectedHypergraph}

Returns hyperedges for a given vertex `v_id` in a directed hypergraph `h`.
Hyperedge indices are given in a tuple `(tail, head)`, where `tail` are the hyperedges where
vertex `v_ind` is on the tail side and `head` are the hyperedges where `v_ind` is on
the head side.

"""
@inline SimpleHypergraphs.gethyperedges(h::H, v_id::Int) where {H <: AbstractDirectedHypergraph} = (h.hg_tail.v2he[v_id], h.hg_head.v2he[v_id])

"""
    to_undirected(h::DirectedHypergraph)

Converts a directed hypergraph into an undirected hypergraph.
Tail and head hyperedges are combined; that is, for all hyperedges he_orig in
the directed hypergraph h, all vertices in the head or tail are added to a
corresponding undirected hyperedge he_new in the undirected hypergraph h'.

Metadata is combined into tuples; i.e., if there was originally tail metadata
t_meta and head metadata h_meta for a given directed hyperedge, the new
undirected hyperedge will have metadata (t_meta, h_meta).

Because vertex-hyperedge weights are restricted to real numbers, we cannot
combine the weights, so we simply set the values to 1.0 if a given vertex
is in a given hyperedge 

"""
function to_undirected(h::DirectedHypergraph{T,V,E,D}) where {T <: Real, V, E, D <: AbstractDict{Int, T}}

    incidence = Matrix{Union{T, Nothing}}(nothing, nhv(h), nhe(h))

    this_nhe = nhe(h)

    for row in 1:nhv(h)
        for column in 1:this_nhe
            tail_val, head_val = h[row, column]
            if tail_val === nothing && head_val === nothing
                incidence[row, column] = nothing
            else
                incidence[row, column] = convert(T, 1.0)
            end
        end
    end

    combined_he_meta = Vector{Union{Tuple{Union{E, Nothing}, Union{E, Nothing}}, Nothing}}(undef, this_nhe)
    fill!(combined_he_meta, nothing)
    for he_index in 1:this_nhe
        tail_meta = h.he_meta_tail[he_index]
        head_meta = h.he_meta_head[he_index]

        if tail_meta !== nothing || head_meta !== nothing
            combined_he_meta[he_index] = (tail_meta, head_meta)
        end
    end

    Hypergraph{T, V, Tuple{Union{E, Nothing},Union{E, Nothing}}, D}(
        incidence,
        v_meta=h.v_meta,
        he_meta=combined_he_meta
    )

end

"""
    to_undirected(h::BasicDirectedHypergraph)

Converts a directed hypergraph into an undirected hypergraph.
Tail and head hyperedges are combined; that is, for all hyperedges he_orig in
the directed hypergraph h, all vertices in the head or tail are added to a
corresponding undirected hyperedge he_new in the undirected hypergraph h'.

Because vertex-hyperedge weights are restricted to real numbers, we cannot
combine the weights, so we simply set the values to 1.0 if a given vertex
is in a given hyperedge 
"""
function to_undirected(h::BasicDirectedHypergraph{T, D}) where {T <: Real, D <: AbstractDict{Int,T}}
    incidence = Matrix{Union{T, Nothing}}(undef, nhv(h), nhe(h))
    fill!(incidence, nothing)

    this_nhe = nhe(h)

    for row in 1:nhv(h)
        for column in 1:this_nhe
            tail_val, head_val = h[row, column]
            if tail_val === nothing && head_val === nothing
                incidence[row, column] = nothing
            else
                incidence[row, column] = convert(T, 1.0)
            end
        end
    end

    BasicHypergraph{T, D}(incidence)

end


"""
    add_vertex!(h::DirectedHypergraph{T, V, E, D};
                hyperedges_tail::D = D(), hyperedges_head::D = D(), v_meta::Union{V,Nothing} = nothing
                ) where {T <: Real, V, E, D <: AbstractDict{Int,T}}

Adds a vertex to a given directed hypergraph `h`. Optionally, the vertex can be added
to existing hyperedges. The `hyperedges_tail` parameter presents a dictionary
of hyperedge identifiers and values stored at the ingoing side of hyperedges, and
the `hyperedges_head` parameter presents a dictionary of hyperedge identifiers and
values stored at the head side of hyperedges.
Additionally, a value can be stored with the vertex using the `v_meta` keyword
parameter.

"""
function SimpleHypergraphs.add_vertex!(h::DirectedHypergraph{T, V, E, D};
                     hyperedges_tail::D = D(), hyperedges_head::D = D(), v_meta::Union{V,Nothing} = nothing
                    ) where {T <: Real, V, E, D <: AbstractDict{Int,T}}
    @boundscheck (checkbounds(h.hg_tail,1,k) for k in keys(hyperedges_tail))
    @boundscheck (checkbounds(h.hg_head,1,k) for k in keys(hyperedges_head))

    push!(h.hg_tail.v2he,hyperedges_tail)
    push!(h.hg_head.v2he,hyperedges_head)

    # Should always be identical to h.hg_head.v2he
    ix = length(h.hg_tail.v2he)

    for k in keys(hyperedges_tail)
        h[1,ix,k]=hyperedges_tail[k]
    end

    for k in keys(hyperedges_head)
        h[2,ix,k]=hyperedges_head[k]
    end

    push!(h.v_meta, v_meta)
    ix
end

"""
    add_vertex!(h::BasicDirectedHypergraph{T, D};
                hyperedges_tail::D = D(), hyperedges_head::D = D()
                ) where {T <: Real, D <: AbstractDict{Int,T}}

Adds a vertex to a given directed hypergraph `h`. Optionally, the vertex can be added
to existing hyperedges. The `hyperedges_tail` parameter presents a dictionary
of hyperedge identifiers and values stored at the ingoing side of hyperedges, and
the `hyperedges_head` parameter presents a dictionary of hyperedge identifiers and
values stored at the head side of hyperedges.

"""
function SimpleHypergraphs.add_vertex!(h::BasicDirectedHypergraph{T, D};
                     hyperedges_tail::D = D(), hyperedges_head::D = D()
                    ) where {T <: Real, D <: AbstractDict{Int,T}}
    @boundscheck (checkbounds(h.hg_tail,1,k) for k in keys(hyperedges_tail))
    @boundscheck (checkbounds(h.hg_head,1,k) for k in keys(hyperedges_head))

    push!(h.hg_tail.v2he,hyperedges_tail)
    push!(h.hg_head.v2he,hyperedges_head)
    
    ix = length(h.hg_tail.v2he)

    for k in keys(hyperedges_tail)
        h[1,ix,k]=hyperedges_tail[k]
    end

    for k in keys(hyperedges_head)
        h[2,ix,k]=hyperedges_head[k]
    end

    ix
end


"""
    remove_vertex!(h::DirectedHypergraph, v::Int)

Removes the vertex `v` from a given directed hypergraph `h`.
Note that running this function will cause reordering of vertices in the
hypergraph; the vertex `v` will replaced by the last vertex of the hypergraph
and the list of vertices will be shrunk.
"""
function SimpleHypergraphs.remove_vertex!(h::DirectedHypergraph, v::Int)
    n = nhv(h)
    if v < n
        h.v_meta[v] = h.v_meta[n]
    end
    
    remove_vertex!(h.hg_tail, v)
    remove_vertex!(h.hg_head, v)

    resize!(h.v_meta, length(h.v_meta) - 1)

    h
end


"""
    remove_vertex!(h::BasicDirectedHypergraph, v::Int)

Removes the vertex `v` from a given directed hypergraph `h`.
Note that running this function will cause reordering of vertices in the
hypergraph; the vertex `v` will replaced by the last vertex of the hypergraph
and the list of vertices will be shrunk.
"""
function SimpleHypergraphs.remove_vertex!(h::BasicDirectedHypergraph, v::Int)
    remove_vertex!(h.hg_tail, v)
    remove_vertex!(h.hg_head, v)

    h
end


"""
    add_hyperedge!(h::DirectedHypergraph{T, V, E, D};
                   vertices_tail::D = D(), vertices_head::D = D(),
                   he_meta_tail::Union{E,Nothing}=nothing, he_meta_head::Union{E,Nothing}=nothing
                   ) where {T <: Real, V, E, D <: AbstractDict{Int,T}}

Adds a hyperedge to a given directed hypergraph `h`.
Optionally, existing vertices can be added to the created hyperedge in the
tail or head directions.
The paramater `vertices_tail` represents a dictionary of vertex identifiers and
values stored at the tail hyperedge; `vertices_head` represented the vertex
identifiers and values stored at the outcoming side of the hyperedge. Additionally, 
a value can be stored with the hyperedge using the `he_meta_tail` and `he_meta_head`
keyword parameters.

"""
function SimpleHypergraphs.add_hyperedge!(h::DirectedHypergraph{T, V, E, D};
                        vertices_tail::D = D(), vertices_head::D = D(),
                        he_meta_tail::Union{E,Nothing}=nothing, he_meta_head::Union{E,Nothing}=nothing
                        ) where {T <: Real, V, E, D <: AbstractDict{Int,T}}
    @boundscheck (checkbounds(h.hg_tail,k,1) for k in keys(vertices_tail))
    @boundscheck (checkbounds(h.hg_head,1,k) for k in keys(vertices_head))
    
    push!(h.hg_tail.he2v,vertices_tail)
    push!(h.hg_head.he2v, vertices_head)

    # Should always be identical to length(h.hg_head.he2v)
    ix = length(h.hg_tail.he2v)
    for k in keys(vertices_tail)
        h[1,k,ix]=vertices_tail[k]
    end
    for k in keys(vertices_head)
        h[2,k,ix]=vertices_head[k]
    end
    push!(h.he_meta_tail, he_meta_tail)
    push!(h.he_meta_head, he_meta_head)
    ix
end


"""
    add_hyperedge!(h::BasicDirectedHypergraph{T, D};
                   vertices::D = D()
                   ) where {T <: Real, D <: AbstractDict{Int,T}}

Adds a hyperedge to a given directed hypergraph `h`.
Optionally, existing vertices can be added to the created hyperedge in the
tail or head directions.
The paramater `vertices_tail` represents a dictionary of vertex identifiers and
values stored at the tail hyperedge; `vertices_head` represented the vertex
identifiers and values stored at the outcoming side of the hyperedge.

"""
function SimpleHypergraphs.add_hyperedge!(h::BasicDirectedHypergraph{T, D};
                        vertices_tail::D = D(), vertices_head::D = D()
                        ) where {T <: Real, D <: AbstractDict{Int,T}}
    @boundscheck (checkbounds(h.hg_tail,k,1) for k in keys(vertices_tail))
    @boundscheck (checkbounds(h.hg_head,1,k) for k in keys(vertices_head))
    
    push!(h.hg_tail.he2v, vertices_tail)
    push!(h.hg_head.he2v, vertices_head)

    # Should always be identical to length(h.hg_head.he2v)
    ix = length(h.hg_tail.he2v)
    for k in keys(vertices_tail)
        h[1,k,ix]=vertices_tail[k]
    end
    for k in keys(vertices_head)
        h[2,k,ix]=vertices_head[k]
    end

    ix
end


"""
    remove_hyperedge!(h::DirectedHypergraph, e::Int)
Removes the hyperedge `e` from a given directed hypergraph `h`.
Note that running this function will cause reordering of hyperedges in the
hypergraph: the hyperedge `e` will replaced by the last hyperedge of the hypergraph
and the list of hyperedges (and hyperedge metadata) will be shrunk.
"""
function SimpleHypergraphs.remove_hyperedge!(h::DirectedHypergraph, e::Int)
    ne = nhe(h)
	@assert(e <= ne)
	if e < ne
	    h.he_meta_tail[e] = h.he_meta_tail[ne]
        h.he_meta_head[e] = h.he_meta_head[ne]
	end

    remove_hyperedge!(h.hg_tail, e)
    remove_hyperedge!(h.hg_head, e)

    resize!(h.he_meta_tail, length(h.he_meta_tail) - 1)
    resize!(h.he_meta_head, length(h.he_meta_head) - 1)

    h
end


"""
    remove_hyperedge!(h::BasicDirectedHypergraph, e::Int)
Removes the hyperedge `e` from a given directed hypergraph `h`.
Note that running this function will cause reordering of hyperedges in the
hypergraph: the hyperedge `e` will replaced by the last hyperedge of the hypergraph
and the list of hyperedges (and hyperedge metadata) will be shrunk.
"""
function SimpleHypergraphs.remove_hyperedge!(h::BasicDirectedHypergraph, e::Int)
    remove_hyperedge!(h.hg_tail, e)
    remove_hyperedge!(h.hg_head, e)

    h
end


"""
    prune_hypergraph!(h::H) where {H <: AbstractDirectedHypergraph}

Remove all vertices with degree 0 and all hyperedges of size 0.

"""
function SimpleHypergraphs.prune_hypergraph!(h::H) where {H <: AbstractDirectedHypergraph}
	for e in reverse(1:nhe(h))
        length(h.hg_tail.he2v[e]) == 0 && length(h.hg_head.he2v[e]) == 0 && SimpleHypergraphs.remove_hyperedge!(h,e)
    end
	for v in reverse(1:nhv(h))
    	length(h.hg_tail.v2he[v]) == 0 && length(h.hg_tail.v2he[v]) == 0 && SimpleHypergraphs.remove_vertex!(h,v)
    end
	h
end


"""
    set_vertex_meta!(h::DirectedHypergraph{T, V, E, D}, new_value::Union{V,Nothing},
        id::Int) where {T <: Real, V, E, D <: AbstractDict{Int,T}}

Sets a new meta value `new_value` for the vertex `id` in the hypergraph `h`.

"""
function SimpleHypergraphs.set_vertex_meta!(h::DirectedHypergraph{T, V, E, D},
        new_value::Union{V,Nothing}, id::Int
        ) where {T <: Real, V, E, D <: AbstractDict{Int,T}}
    checkbounds(h.v_meta, id)
    h.v_meta[id] = new_value
    h.v_meta
end


SimpleHypergraphs.set_vertex_meta!(::BasicDirectedHypergraph, ::Any, ::Int) = throw("Not implemented!")


"""
    get_vertex_meta(h::DirectedHypergraph{T, V, E, D}, id::Int
                    ) where {T <: Real, V, E, D <: AbstractDict{Int,T}}

Returns a meta value stored at the vertex `id` in the directed hypergraph `h`.

"""
function SimpleHypergraphs.get_vertex_meta(h::DirectedHypergraph{T, V, E, D}, id::Int
                         ) where {T <: Real, V, E, D <: AbstractDict{Int,T}}
    checkbounds(h.v_meta, id)
    h.v_meta[id]
end


SimpleHypergraphs.get_vertex_meta(::BasicDirectedHypergraph, ::Int) = throw("Not implemented!")


"""
    set_hyperedge_meta!(h::DirectedHypergraph{T, V, E, D},
        new_value_tail::Union{E,Nothing}, new_value_head::Union{E,Nothing}, id::Int
        ) where {T <: Real, V, E, D <: AbstractDict{Int,T}}

Sets a new meta value `new_value` for the hyperedge `id` in the directed hypergraph `h`.

"""
function SimpleHypergraphs.set_hyperedge_meta!(h::DirectedHypergraph{T, V, E, D},
                             new_value_tail::Union{E,Nothing}, new_value_head::Union{E,Nothing}, id::Int
                             ) where {T <: Real, V, E, D <: AbstractDict{Int,T}}
    checkbounds(h.he_meta_tail, id)
    checkbounds(h.he_meta_head, id)

    h.he_meta_tail[id] = new_value_tail
    h.he_meta_head[id] = new_value_head

    (h.he_meta_tail, h.he_meta_head)
end


"""
    set_hyperedge_meta!(h::DirectedHypergraph{T, V, E, D},
        new_value::Union{E,Nothing}, id::Int, side::HyperedgeDirection
        ) where {T <: Real, V, E, D <: AbstractDict{Int,T}}

Sets a new meta value `new_value` for the hyperedge `id` in the direction `side`
    for a directed hypergraph `h`.

"""
function SimpleHypergraphs.set_hyperedge_meta!(h::DirectedHypergraph{T, V, E, D},
                             new_value::Union{E,Nothing}, id::Int, side::HyperedgeDirection
                             ) where {T <: Real, V, E, D <: AbstractDict{Int,T}}
    
    if side == tail
        checkbounds(h.he_meta_tail, id)
        h.he_meta_tail[id] = new_value
        h.he_meta_tail
    else
        checkbounds(h.he_meta_head, id)
        h.he_meta_head[id] = new_value
        h.he_meta_head
    end

end


SimpleHypergraphs.set_hyperedge_meta!(::BasicDirectedHypergraph, ::Any, ::Any, ::Int) = throw("Not implemented!")
SimpleHypergraphs.set_hyperedge_meta!(::BasicDirectedHypergraph, ::Any, ::Int, ::HyperedgeDirection) = throw("Not implemented!")


"""
    get_hyperedge_meta(h::DirectedHypergraph{T, V, E, D}, id::Int)
        where {T <: Real, V, E, D <: AbstractDict{Int,T}}
Returns a meta value stored at the hyperedge `id` in the directed hypergraph `h`.

"""
function SimpleHypergraphs.get_hyperedge_meta(h::DirectedHypergraph{T, V, E, D}, id::Int
                            ) where {T <: Real, V, E, D <: AbstractDict{Int,T}}
    checkbounds(h.he_meta_tail, id)
    checkbounds(h.he_meta_head, id)

    (h.he_meta_tail[id], h.he_meta_head[id])
end


"""
    get_hyperedge_meta(h::DirectedHypergraph{T, V, E, D}, id::Int, side::HyperedgeDirection)
        where {T <: Real, V, E, D <: AbstractDict{Int,T}}
Returns a meta value stored at the hyperedge `id` in the directed hypergraph `h`.

"""
function SimpleHypergraphs.get_hyperedge_meta(h::DirectedHypergraph{T, V, E, D}, id::Int, side::HyperedgeDirection
                            ) where {T <: Real, V, E, D <: AbstractDict{Int,T}}

    if side == tail
        checkbounds(h.he_meta_tail, id)
        h.he_meta_tail[id]
    else
        checkbounds(h.he_meta_head, id)
        h.he_meta_head[id]
    end
end


SimpleHypergraphs.get_hyperedge_meta(::BasicDirectedHypergraph, ::Int) = throw("Not implemented!")
SimpleHypergraphs.get_hyperedge_meta(::BasicDirectedHypergraph, ::Int, ::HyperedgeDirection) = throw("Not implemented!")


"""
    nhe(h::H) where {H <: AbstractDirectedHypergraph}

Return the number of hyperedges in the directed hypergraph `h`.
"""
function SimpleHypergraphs.nhe(h::H) where {H <: AbstractDirectedHypergraph}
    (length(h.hg_tail.he2v) == length(h.hg_head.he2v)) ? length(h.hg_tail.he2v) : throw("Tail and head sides of hypergraph have different numbers of hyperedges!")
end


"""
    nhv(h::H) where {H <: AbstractDirectedHypergraph}

Return the number of vertices in the directed hypergraph `h`.
"""
function SimpleHypergraphs.nhv(h::H) where {H <: AbstractDirectedHypergraph}
    (length(h.hg_tail.v2he) == length(h.hg_head.v2he)) ? length(h.hg_tail.v2he) : throw("Tail and head sides of hypergraph have different numbers of hyperedges!")
end


function _default_heselect(h::H, v::Int; reverse::Bool=false) where {H <: AbstractDirectedHypergraph}
    he_tail, he_head = gethyperedges(h, v)

    if reverse
        hes = he_head
    else
        hes = he_tail
    end

    sort!(collect(keys(hes))), ones(length(hes))
end


function _default_vselect(h::H, he::Int; reverse::Bool=false) where {H <: AbstractDirectedHypergraph}
    vs_tail, vs_head = getvertices(h, he)

    if reverse
        vs = vs_tail
    else
        vs = vs_head
    end

    sort!(collect(keys(vs))), ones(length(vs))

end


"""
    random_walk(
        h::H,
        start::Int;
        heselect::Function,
        vselect::Function,
        reverse::bool
    ) where {H <: AbstractDirectedHypergraph}

Return a next vertex visited in assuming a random walk starting from vertex `start`.
First a hyperedge is sampled with weights proportional to `heselect` function
(by default each hyperedge is sampled with the same probability).
Next a vertex within hyperedge is with weights proportional to `vselect` function
(by default each vertex, including the source, is sampled with the same probability).

`heselect` and `vselect` functions take two arguments a `Hypergraph` and respectively
a vertex identifier or a hyperedge identifier. The return values of both functions
should be respectively a list of hyperedges or vertices and their weights.
"""
function SimpleHypergraphs.random_walk(h::H, start::Int;
                     heselect::Function=_default_heselect,
                     vselect::Function=_default_vselect,
                     reverse::Bool=false) where {H <: AbstractDirectedHypergraph}
    1 <= start <= nhv(h) || throw(ArgumentError("invalid start vertex index"))
    hes, hew = heselect(h, start, reverse=reverse)
    he = sample(hes, Weights(hew))
    ves, vw = vselect(h, he, reverse=reverse)
    return sample(ves, Weights(vw))
end


"""
    get_weakly_connected_components(h::H) where {H <: AbstractDirectedHypergraph}

Return an array of weakly connected components in the directed hypergraph `h`
(array of vectors of vertices) by first converting the directed hypergraph
into an undirected hypergraph and then obtaining the conected components of
that hypergraph.
"""
function get_weakly_connected_components(h::H) where {H <: AbstractDirectedHypergraph}
    undirected = to_undirected(h)
    get_connected_components(undirected)
end


"""
    _visit!(h::H, v::Int) where {H <: AbstractDirectedHypergraph}

Determines the B-connected component of a vertex `v` in directed hypergraph `h`.
This is an auxiliary function for `get_strongly_connected_components`, which
determines the strongly connected components of a directed hypergraph.
"""
function _visit(
    h::H,
    v::Int
) where {H <: AbstractDirectedHypergraph}
    visited = zeros(Bool, nhv(h))
    visited_tail_nodes = zeros(Int, nhe(h))

    q = Queue{Int}()
    bcc = Set{Int}()
    enqueue!(q, v)

    visited[v] = true

    while length(q) > 0
        u = dequeue!(q)
        push!(bcc, u)

        tail_hes = gethyperedges(h, u)[1]

        for tail_he in keys(tail_hes)
            visited_tail_nodes[tail_he] += 1

            tail_vs, head_vs = getvertices(h, tail_he)

            if visited_tail_nodes[tail_he] == length(tail_vs)
                for head_v in keys(head_vs)
                    if !visited[head_v]
                        visited[head_v] = true
                        enqueue!(q, head_v)
                    end
                end
            end
        end
    end

    bcc
end


"""
    get_strongly_connected_components(h::H) where {H <: AbstractDirectedHypergraph}

Return an array of strongly connected components in the directed hypergraph `h`
(array of vectors of vertices), based on the "naive" algorithm of
Francisco José Martín-Recuerda Moyano (PhD dissertation, 2016).

"""
function get_strongly_connected_components(h::H) where {H <: AbstractDirectedHypergraph}

    T = Dict{Vector{Int}, Set{Int}}()

    for v in 1:nhv(h)
        bcc_v = _visit(h, v)
        bcc_sorted = sort(collect(bcc_v))
        for i in 1:length(bcc_sorted)
            if !haskey(T, bcc_sorted[1:i])
                T[bcc_sorted[1:i]] = Set{Int}()
            end
        end
        push!(T[bcc_sorted], v)
    end

    [sort!(collect(v)) for (k, v) in T if length(v) != 0]
end
