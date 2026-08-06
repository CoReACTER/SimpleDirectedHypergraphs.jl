"""
    DirectedHypergraph{T} <: AbstractDirectedHypergraph{Tuple{Union{T, Nothing}, Union{T, Nothing}}}

A directed hypergraph (dihypergraph) storing information about vertices and directed hyperedges
(dihyperedges).

This implementation is based on guidance from Przemysław Szufel;
    see https://github.com/pszufe/SimpleHypergraphs.jl/issues/45
This allows us to manipulate `DirectedHypergraphs` using `Hypergraph` functionality
There is danger of a user manipulating individual `hg_tail` and `hg_head` (undirected) hypergraphs

Is there a smart way to prevent this?
TODO: reconsider this design choice

**Constructors**

    DirectedHypergraph{T,V,E,D}(
        n::Integer, k::Integer,
        v_meta=Vector{Union{V, Nothing}}(nothing, n),
        he_meta_tail=Vector{Union{E, Nothing}}(nothing, k),
        he_meta_head=Vector{Union{E, Nothing}}(nothing, k)
    ) where {T<:Real,V,E,D<:AbstractDict{Int, T}}
    DirectedHypergraph{T,V,E}(n::Integer, k::Integer) where {T<:Real, V, E}
    DirectedHypergraph{T,V}(n::Integer, k::Integer) where {T<:Real, V}
    DirectedHypergraph{T}(n::Integer, k::Integer) where {T<:Real}
    DirectedHypergraph(n::Integer, k::Integer)

Construct a dihypergraph with a given number of vertices and hyperedges.
Optionally, values of type `V` can be stored at vertices and values of type `E`
can be stored at hyperedges. By default the dihypergraph uses a `Dict{Int,T}` for
the internal data storage; however, a different dictionary such as `SortedDict`
to ensure result replicability can be used (e.g., when doing stochastic
simulations on dihypergraphs).

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

Construct a dihypergraph using its matrix representation. In the matrix representation rows are
vertices and columns are hyperedges. Optionally, values of type `V` can be stored at vertices and
values of type `E` can be stored at hyperedges. By default the hypergraph uses a `Dict{Int,T}` for
the internal data storage, however a different dictionary such as `SortedDict` to ensure result
replicability can be used (e.g. when doing stochastic simulations on dihypergraphs).

    DirectedHypergraph(g::Graphs.DiGraph)

Constructs a dihypergraph of degree 2 by making a deep copy of a `Graphs.DiGraph`. A `SortedDict`
will be used for internal data storage of the dihypergraph.

    DirectedHypergraph(
        hg_tail::Hypergraph{T},
        hg_head::Hypergraph{T}
    ) where {T<:Real}
    DirectedHypergraph{T}(
        hg_tail::Hypergraph{T},
        hg_head::Hypergraph{T}
    ) where {T<:Real}
    DirectedHypergraph{T,V}(
        hg_tail::Hypergraph{T},
        hg_head::Hypergraph{T};
        v_meta::Union{Nothing, Vector{Union{V, Nothing}}}=nothing,
    ) where {T<:Real,V}
     DirectedHypergraph{T,V,E}(
        hg_tail::Hypergraph{T},
        hg_head::Hypergraph{T};
        v_meta::Union{Nothing, Vector{Union{V, Nothing}}}=nothing,
        he_meta_tail::Union{Nothing, Vector{Union{E, Nothing}}}=nothing,
        he_meta_head::Union{Nothing, Vector{Union{E, Nothing}}}=nothing
    ) where {T<:Real,V,E}
     DirectedHypergraph{T,V,E,D}(
        hg_tail::Hypergraph{T,Nothing,Nothing,D},
        hg_head::Hypergraph{T,Nothing,Nothing,D};
        v_meta::Vector{Union{Nothing,V}}=Vector{Union{Nothing,V}}(nothing, size(hg_tail,1)),
        he_meta_tail::Vector{Union{Nothing,E}}=Vector{Union{Nothing,E}}(nothing, size(hg_tail,2)),
        he_meta_head::Vector{Union{Nothing,E}}=Vector{Union{Nothing,E}}(nothing, size(hg_tail,2))
    ) where {T<:Real,V,E,D<:AbstractDict{Int, T}}

Constructs a dihypergraph from two undirected hypergraphs, one with hyperedges containing "tail"
vertices and one with hyperedges containing "head" verticies.

    DirectedHypergraph{T,V,E,D}(
        hg_tail::Hypergraph{T,V,E,D},
        hg_head::Hypergraph{T,V,E,D}
    ) where {T<:Real,V,E,D<:AbstractDict{Int, T}}

Constructs a dihypergraph from two hypergraphs potentially containing metadata. Throws an error if
the vertex metadata of the two hypergraphs is not element-for-element identical.

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
struct DirectedHypergraph{T <: Real, V, E, D <: AbstractDict{Int, T}} <: AbstractDirectedHypergraph{Tuple{Union{T, Nothing}, Union{T, Nothing}}}
    hg_tail::Hypergraph{T, Nothing, Nothing, D}
    hg_head::Hypergraph{T, Nothing, Nothing, D}

    v_meta::Vector{Union{V, Nothing}}
    he_meta_tail::Vector{Union{E, Nothing}}
    he_meta_head::Vector{Union{E, Nothing}}

    DirectedHypergraph{T, V, E, D}(
        n::Integer, k::Integer,
        v_meta = Vector{Union{V, Nothing}}(nothing, n),
        he_meta_tail = Vector{Union{E, Nothing}}(nothing, k),
        he_meta_head = Vector{Union{E, Nothing}}(nothing, k)
    ) where {T <: Real, V, E, D <: AbstractDict{Int, T}} =
        new{T, V, E, D}(
        Hypergraph{T, Nothing, Nothing, D}(n, k),
        Hypergraph{T, Nothing, Nothing, D}(n, k),
        v_meta, he_meta_tail, he_meta_head
    )

    function DirectedHypergraph{T, V, E, D}(
            hg_tail::Hypergraph{T, Nothing, Nothing, D},
            hg_head::Hypergraph{T, Nothing, Nothing, D};
            v_meta::Vector{Union{Nothing, V}} = Vector{Union{Nothing, V}}(nothing, size(hg_tail, 1)),
            he_meta_tail::Vector{Union{Nothing, E}} = Vector{Union{Nothing, E}}(nothing, size(hg_tail, 2)),
            he_meta_head::Vector{Union{Nothing, E}} = Vector{Union{Nothing, E}}(nothing, size(hg_tail, 2))
        ) where {T <: Real, V, E, D <: AbstractDict{Int, T}}
        @assert size(hg_tail) == size(hg_head)

        @assert length(v_meta) == size(hg_tail, 1)
        @assert length(he_meta_tail) == size(hg_tail, 2)
        @assert length(he_meta_head) == size(hg_head, 2)

        return new{T, V, E, D}(
            hg_tail,
            hg_head,
            v_meta,
            he_meta_tail,
            he_meta_head
        )
    end
end


DirectedHypergraph{T, V, E}(n::Integer, k::Integer) where {T <: Real, V, E} = DirectedHypergraph{T, V, E, Dict{Int, T}}(n, k)

DirectedHypergraph{T, V}(n::Integer, k::Integer) where {T <: Real, V} = DirectedHypergraph{T, V, Nothing, Dict{Int, T}}(n, k)

DirectedHypergraph{T}(n::Integer, k::Integer) where {T <: Real} = DirectedHypergraph{T, Nothing, Nothing, Dict{Int, T}}(n, k)

DirectedHypergraph(n::Integer, k::Integer) = DirectedHypergraph{Bool, Nothing, Nothing, Dict{Int, Bool}}(n, k)


function DirectedHypergraph{T, V, E}(
        hg_tail::Hypergraph{T},
        hg_head::Hypergraph{T};
        v_meta::Union{Nothing, Vector{Union{V, Nothing}}} = nothing,
        he_meta_tail::Union{Nothing, Vector{Union{E, Nothing}}} = nothing,
        he_meta_head::Union{Nothing, Vector{Union{E, Nothing}}} = nothing
    ) where {T <: Real, V, E}
    @assert size(hg_tail) == size(hg_head)

    n, k = size(hg_tail)
    shg_tail = Hypergraph{T}(n, k)
    shg_head = Hypergraph{T}(n, k)

    # TODO: test behavior on this
    shg_tail .= hg_tail
    shg_head .= hg_head

    if E === Nothing
        he_meta_tail = fill(nothing, size(hg_tail, 2))
        he_meta_head = fill(nothing, size(hg_tail, 2))
    else
        if isnothing(he_meta_tail)
            he_meta_tail = hg_tail.he_meta
        end
        if isnothing(he_meta_head)
            he_meta_head = hg_head.he_meta
        end
    end

    if V === Nothing
        v_meta = fill(nothing, size(hg_tail, 1))
    end


    return if isnothing(v_meta)
        if !all(hg_tail.v_meta .== hg_head.v_meta)
            @warn "Vertex metadata for tail and head hypergraphs not identical; discarding vertex metadata."
            DirectedHypergraph{T, V, E, Dict{Int, T}}(
                shg_tail,
                shg_head;
                he_meta_tail = he_meta_tail,
                he_meta_head = he_meta_head
            )
        else
            DirectedHypergraph{T, V, E, Dict{Int, T}}(
                shg_tail,
                shg_head;
                v_meta = hg_tail.v_meta,
                he_meta_tail = he_meta_tail,
                he_meta_head = he_meta_head
            )
        end
    else
        DirectedHypergraph{T, V, E, Dict{Int, T}}(
            shg_tail,
            shg_head;
            v_meta = v_meta,
            he_meta_tail = he_meta_tail,
            he_meta_head = he_meta_head
        )
    end
end

function DirectedHypergraph{T, V}(
        hg_tail::Hypergraph{T},
        hg_head::Hypergraph{T};
        v_meta::Union{Nothing, Vector{Union{V, Nothing}}} = nothing,
    ) where {T <: Real, V}

    return DirectedHypergraph{T, V, Nothing}(
        hg_tail,
        hg_head;
        v_meta = v_meta
    )
end

function DirectedHypergraph{T}(
        hg_tail::Hypergraph{T},
        hg_head::Hypergraph{T}
    ) where {T <: Real}

    return DirectedHypergraph{T, Nothing, Nothing}(
        hg_tail,
        hg_head
    )
end

function DirectedHypergraph(
        hg_tail::Hypergraph{T},
        hg_head::Hypergraph{T}
    ) where {T <: Real}

    return DirectedHypergraph{T}(
        hg_tail,
        hg_head
    )
end


function DirectedHypergraph{T, V, E, D}(
        m_tail::AbstractMatrix{Union{T, Nothing}},
        m_head::AbstractMatrix{Union{T, Nothing}};
        v_meta::Vector{Union{Nothing, V}} = Vector{Union{Nothing, V}}(nothing, size(m_tail, 1)),
        he_meta_tail::Vector{Union{Nothing, E}} = Vector{Union{Nothing, E}}(nothing, size(m_tail, 2)),
        he_meta_head::Vector{Union{Nothing, E}} = Vector{Union{Nothing, E}}(nothing, size(m_tail, 2))
    ) where {T <: Real, V, E, D <: AbstractDict{Int, T}}

    @assert size(m_tail) == size(m_head)

    # Arbitrary, since sizes are identical
    n, k = size(m_tail)

    hg_tail = Hypergraph{T, Nothing, Nothing, D}(n, k)
    hg_tail .= m_tail

    hg_head = Hypergraph{T, Nothing, Nothing, D}(n, k)
    hg_head .= m_head

    return DirectedHypergraph{T, V, E, D}(
        hg_tail, hg_head;
        v_meta = v_meta,
        he_meta_tail = he_meta_tail,
        he_meta_head = he_meta_head
    )
end

function DirectedHypergraph{T, V, E}(
        m_tail::AbstractMatrix{Union{T, Nothing}},
        m_head::AbstractMatrix{Union{T, Nothing}};
        v_meta::Vector{Union{Nothing, V}} = Vector{Union{Nothing, V}}(nothing, size(m_tail, 1)),
        he_meta_tail::Vector{Union{Nothing, E}} = Vector{Union{Nothing, E}}(nothing, size(m_tail, 2)),
        he_meta_head::Vector{Union{Nothing, E}} = Vector{Union{Nothing, E}}(nothing, size(m_tail, 2))
    ) where {T <: Real, V, E}

    # Arbitrary, since sizes are identical
    n, k = size(m_tail)

    hg_tail = Hypergraph{T}(n, k)
    hg_tail .= m_tail

    hg_head = Hypergraph{T}(n, k)
    hg_head .= m_head

    return DirectedHypergraph{T, V, E}(
        hg_tail, hg_head;
        v_meta = v_meta,
        he_meta_tail = he_meta_tail,
        he_meta_head = he_meta_head
    )
end

function DirectedHypergraph{T, V}(
        m_tail::AbstractMatrix{Union{T, Nothing}},
        m_head::AbstractMatrix{Union{T, Nothing}};
        v_meta::Vector{Union{Nothing, V}} = Vector{Union{Nothing, V}}(nothing, size(m_tail, 1)),
    ) where {T <: Real, V}

    # Arbitrary, since sizes are identical
    n, k = size(m_tail)

    hg_tail = Hypergraph{T}(n, k)
    hg_tail .= m_tail

    hg_head = Hypergraph{T}(n, k)
    hg_head .= m_head

    return DirectedHypergraph{T, V}(
        hg_tail,
        hg_head;
        v_meta = v_meta
    )
end

function DirectedHypergraph{T}(
        m_tail::AbstractMatrix{Union{T, Nothing}},
        m_head::AbstractMatrix{Union{T, Nothing}}
    ) where {T <: Real}

    # Arbitrary, since sizes are identical
    n, k = size(m_tail)

    hg_tail = Hypergraph{T}(n, k)
    hg_tail .= m_tail

    hg_head = Hypergraph{T}(n, k)
    hg_head .= m_head

    return DirectedHypergraph{T}(hg_tail, hg_head)
end

function DirectedHypergraph(
        m_tail::AbstractMatrix{Union{T, Nothing}},
        m_head::AbstractMatrix{Union{T, Nothing}}
    ) where {T <: Real}

    # Arbitrary, since sizes are identical
    n, k = size(m_tail)

    hg_tail = Hypergraph{T}(n, k)
    hg_tail .= m_tail

    hg_head = Hypergraph{T}(n, k)
    hg_head .= m_head

    return DirectedHypergraph{T}(hg_tail, hg_head)
end


function DirectedHypergraph(g::Graphs.DiGraph)
    h = DirectedHypergraph{Bool, Nothing, Nothing, SortedDict{Int, Bool}}(maximum(vertices(g)), ne(g))
    e = 0
    for edge in edges(g)
        e += 1
        h[1, edge.src, e] = true
        h[2, edge.dst, e] = true
    end
    return h
end


# TODO: this is awkward...
const DIRECTED_HYPERGRAPH_VALID_FIRST_INDICES = [1, 2]

# TODO: can this entirely replace the above? Index setting seems problematic...
@enum HyperedgeDirection begin
    tail = 1
    head = 2
end


# AbstractArray interface functions

"""
    Base.getindex(h::H, idx::Vararg{Int,2}) where {H <: AbstractDirectedHypergraph}

Returns a value for a given vertex-dihyperedge pair `idx` for a dihypergraph `h`.
If a vertex does not belong to a dihyperedge `nothing` is returned.

**Arguments**

* `h` : Dihypergraph
* `idx` : Pair of integers `i` and `j`, where `i` is a vertex index and `j` is a dihyperedge index

**Returns**

`(tail_value, head_value)`, where `tail_value` and `head_value` are either the weight of vertex `i`
in the tail and head of dihyperedge `j`, respectively, or `nothing` if `i` is not in the tail/head.

"""
@inline function Base.getindex(h::H, idx::Vararg{Int, 2}) where {H <: AbstractDirectedHypergraph}
    @boundscheck checkbounds(h.hg_tail, idx...)
    @boundscheck checkbounds(h.hg_head, idx...)

    tail_value = get(h.hg_tail.v2he[idx[1]], idx[2], nothing)
    head_value = get(h.hg_head.v2he[idx[1]], idx[2], nothing)

    return (tail_value, head_value)
end

"""
    Base.setindex!(h::H, ::Nothing, idx::Vararg{Int,2}) where {H <: AbstractDirectedHypergraph}

Removes a vertex from a given hyperedge for a dihypergraph `h` and a given vertex-dihyperedge pair
`idx`. Note that trying to remove a vertex from a dihyperedge when it is not present will not throw
an error.

**Arguments**

* `h` : Dihypergraph
* `nothing`
* `idx` : Pair of integers `i` and `j`, where `i` is a vertex index and `j` is a dihyperedge index

**Returns**

`h`

"""
@inline function Base.setindex!(h::H, ::Nothing, idx::Vararg{Int, 2}) where {H <: AbstractDirectedHypergraph}
    @boundscheck checkbounds(h.hg_tail, idx...)
    @boundscheck checkbounds(h.hg_head, idx...)
    setindex!(h.hg_tail, nothing, idx...)
    setindex!(h.hg_head, nothing, idx...)
    return h
end


"""
    Base.setindex!(h::H, v::Real, idx::Vararg{Int,2}) where {H <: AbstractDirectedHypergraph}

Adds a vertex to a hyperedge (represented by indices `idx`) and assigns value
`v` to be stored with that assignment.

**Arguments**

* `h` : Dihypergraph
* `v` : A real-valued weight which will be applied to a vertex in both the tail and the head of a
    dihyperedge
* `idx` : Pair of integers `i` and `j`, where `i` is a vertex index and `j` is a dihyperedge index

**Returns**

`h`

"""
@inline function Base.setindex!(h::H, v::Real, idx::Vararg{Int, 2}) where {H <: AbstractDirectedHypergraph}
    @boundscheck checkbounds(h.hg_tail, idx...)
    @boundscheck checkbounds(h.hg_head, idx...)

    setindex!(h.hg_tail, v, idx...)
    setindex!(h.hg_head, v, idx...)
    return h
end


"""
    Base.setindex!(h::H, v::Tuple{Union{Real, Nothing}, Union{Real, Nothing}}, idx::Vararg{Int,2}) where {H <: AbstractDirectedHypergraph}

Manipulates a dihyperedge (represented by indices `idx`), either adding a vertex to the tail and/or
head sides of the hyperedge and assigning a value associated with that assignment, or else removing
a vertex from the tail/head sides of the hyperedge.

Here, `v` is a 2-tuple where the first element is the value that will be assigned to the tail part
of the hyperedge and the second element is the value that will be assigned to the head part. A
value of `nothing` means that the vertex will be removed from that side of the hyperedge.

**Arguments**

* `h` : Dihypergraph
* `v` : A 2-tuple `(tail_value, head_value)`. `tail_value` and `head_value` are either `nothing`
    (in which case the vertex is not in the tail or head, respectively, of the dihyperedge) or
    real-valued weights which will be applied to a vertex in the tail and head of a dihyperedge,
    respectively
* `idx` : Pair of integers `i` and `j`, where `i` is a vertex index and `j` is a dihyperedge index

**Returns**

`h`

"""
@inline function Base.setindex!(h::H, v::Tuple{Union{Real, Nothing}, Union{Real, Nothing}}, idx::Vararg{Int, 2}) where {H <: AbstractDirectedHypergraph}
    @boundscheck checkbounds(h.hg_tail, idx...)
    @boundscheck checkbounds(h.hg_head, idx...)

    setindex!(h.hg_tail, v[1], idx...)
    setindex!(h.hg_head, v[2], idx...)

    return h
end


"""
    Base.setindex!(h::H, ::Nothing, idx::Vararg{Int,3}) where {H <: AbstractDirectedHypergraph}

Removes a vertex from a given hyperedge for a dihypergraph `h` and a given side-vertex-dihyperedge
triple `idx`. If the first index of `idx` is 1, then the vertex will be removed from the tail of
the dihyperedge; if `idx` is 2, then the vertex will be removed from the head. Note that trying to
remove a vertex from a hyperedge when it is not present will not throw an error.

**Arguments**

* `h` : Dihypergraph
* `nothing`
* `idx` : Three-tuple of integers `side`, `i, and `j`, where `side` is either `1` (meaning "tail")
    or `2` (meaning "head"), `i` is a vertex index, and `j` is a dihyperedge index

**Returns**

`h`

"""
@inline function Base.setindex!(h::H, ::Nothing, idx::Vararg{Int, 3}) where {H <: AbstractDirectedHypergraph}
    @boundscheck checkbounds(DIRECTED_HYPERGRAPH_VALID_FIRST_INDICES, idx[1])

    if idx[1] == 1
        side = h.hg_tail
    else
        side = h.hg_head
    end

    @boundscheck checkbounds(side, idx[2:end]...)

    setindex!(side, nothing, idx[2], idx[3])

    return h
end


"""
    Base.setindex!(h::H, v::Real, idx::Vararg{Int,3}) where {H <: AbstractDirectedHypergraph}

Adds a vertex to a dihyperedge (represented by indices `idx`, where the first index must be either
`1` - referring to the tail of the dihyperedge - or `2` - referring to the head) and assigns value
`v` to be stored with that assignment.

**Arguments**

* `h` : Dihypergraph
* `v` : A real-valued weight
* `idx` : 3-tuple of integers `side`, `i`, and `j`, where `side` is either `1` ("tail") or `2`
    ("head"), `i` is a vertex index, and `j` is a dihyperedge index

**Returns**

`h`

"""
@inline function Base.setindex!(h::H, v::Real, idx::Vararg{Int, 3}) where {H <: AbstractDirectedHypergraph}
    @boundscheck checkbounds(DIRECTED_HYPERGRAPH_VALID_FIRST_INDICES, idx[1])

    if idx[1] == 1
        side = h.hg_tail
    else
        side = h.hg_head
    end

    @boundscheck checkbounds(side, idx[2:end]...)

    setindex!(side, v, idx[2]..., idx[3]...)

    return h
end


"""
    getvertices(h::H, he_id::Int) where {H <: AbstractDirectedHypergraph}

Returns vertices (with associated weights) from a directed hypergraph `h` for a given hyperedge
with index `he_id`.

**Arguments**

* `h` : Dihypergraph
* `he_id` : Dihyperedge index

**Returns**

2-tuple `(in, out)`, where `in` contains information about the dihyperedge's tail and `out`
contains information about the head. `in` and `out` are both `AbstractDict`s (following the dict
(`D`) type in `h`), with the form `D(v => w)`, where `v` is a vertex index and `w` is a real-valued
weight.

"""
@inline SimpleHypergraphs.getvertices(h::H, he_id::Int) where {H <: AbstractDirectedHypergraph} = (h.hg_tail.he2v[he_id], h.hg_head.he2v[he_id])


"""
    gethyperedges(h::H, v_id::Int) where {H <: AbstractDirectedHypergraph}

Returns dihyperedges (with associated weights) for a given vertex `v_id` in a dihypergraph `h`.

**Arguments**

* `h` : Dihypergraph
* `v_id` : Vertex index

**Returns**

2-tuple `(in, out)`, where `in` contains information about the dihyperedges where vertex `v_id` is
in the tail and `out` contains information about the dihyperedges where `v_id` is in the head. `in`
and `out` are both `AbstractDict`s (following the dict (`D`) type in `h`), with the form
`D(he => w)`, where `he` is a dihyperedge index and `w` is a real-valued weight.

"""
@inline SimpleHypergraphs.gethyperedges(h::H, v_id::Int) where {H <: AbstractDirectedHypergraph} = (h.hg_tail.v2he[v_id], h.hg_head.v2he[v_id])

"""
    to_undirected(
	h::DirectedHypergraph{T, V, E, D}
    ) where {T <: Real, V, E, D <: AbstractDict{Int, T}}

Converts a dihypergraph into an undirected hypergraph. Dihyperedge tails and heads are combined;
that is, for all dihyperedges `he_orig = (tail, head)` in the dihypergraph `h`, all vertices in
`tail` or `head` are added to a corresponding undirected hyperedge `he_new` in the undirected
hypergraph `h'`.

Metadata is combined into tuples; i.e., if there was originally tail metadata `t_meta` and head
metadata `h_meta` for a given dihyperedge `he_orig`, the new hyperedge `he_new` will have metadata
`(t_meta, h_meta)`.

Because vertex-hyperedge weights are restricted to real numbers, we cannot combine the weights, so
we simply set the values to 1.0 if a given vertex is in a given hyperedge.

**Arguments**

* `h` : Dihypergraph

**Returns**

An undirected hypergraph `h'`

"""
function to_undirected(h::DirectedHypergraph{T, V, E, D}) where {T <: Real, V, E, D <: AbstractDict{Int, T}}

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

    return Hypergraph{T, V, Tuple{Union{E, Nothing}, Union{E, Nothing}}, D}(
        incidence,
        v_meta = h.v_meta,
        he_meta = combined_he_meta
    )

end


"""
    add_vertex!(
	h::DirectedHypergraph{T, V, E, D};
        dihyperedges_tail::D = D(),
	dihyperedges_head::D = D(),
	v_meta::Union{V,Nothing} = nothing
    ) where {T <: Real, V, E, D <: AbstractDict{Int,T}}

Adds a vertex to a given dihypergraph `h`. Optionally, the vertex can be added to existing
hyperedges. The `dihyperedges_tail` parameter presents a dictionary of hyperedge indices and weights
stored at tail side of dihyperedges, and the `dihyperedges_head` parameter presents a dictionary of
hyperedge identifiers and weights stored at the head side of dihyperedges. Additionally, a value
can be stored with the vertex using the `v_meta` keyword parameter.

**Arguments**

* `h` : Dihypergraph
* `dihyperedges_tail` : Dihyperedge index-weight dictionary with information about the dihyperedges
    where the new vertex will be in the tail; default is an empty dict of type `D`
* `dihyperedges_head` : Dihyperedge index-weight dictionary with information about the dihyperedges
    where the new vertex will be in the head; default is an empty dict of type `D`
* `v_meta` : (Optional) metadata for the new vertex; default is `nothing`

**Returns**

`ix`, the index of the new vertex

"""
function SimpleHypergraphs.add_vertex!(
        h::DirectedHypergraph{T, V, E, D};
        dihyperedges_tail::D = D(),
        dihyperedges_head::D = D(),
        v_meta::Union{V, Nothing} = nothing
    ) where {T <: Real, V, E, D <: AbstractDict{Int, T}}
    @boundscheck (checkbounds(h.hg_tail, 1, k) for k in keys(dihyperedges_tail))
    @boundscheck (checkbounds(h.hg_head, 1, k) for k in keys(dihyperedges_head))

    push!(h.hg_tail.v2he, dihyperedges_tail)
    push!(h.hg_head.v2he, dihyperedges_head)

    # Should always be identical to h.hg_head.v2he
    ix = length(h.hg_tail.v2he)

    for k in keys(dihyperedges_tail)
        h[1, ix, k] = dihyperedges_tail[k]
    end

    for k in keys(dihyperedges_head)
        h[2, ix, k] = dihyperedges_head[k]
    end

    push!(h.v_meta, v_meta)
    return ix
end

"""
    remove_vertex!(h::DirectedHypergraph, v::Int)

Removes the vertex `v` from a given dihypergraph `h`. Note that running this function will cause
reordering of vertices in the dihypergraph; the vertex `v` will replaced by the last vertex of the
dihypergraph and the list of vertices will be shrunk.

**Arguments**

* `h` : Dihypergraph
* `v` : Vertex index

**Returns**

`h`

"""
function SimpleHypergraphs.remove_vertex!(h::DirectedHypergraph, v::Int)
    n = nhv(h)
    if v < n
        h.v_meta[v] = h.v_meta[n]

        h.hg_tail.v2he[v] = h.hg_tail.v2he[n]
        h.hg_head.v2he[v] = h.hg_head.v2he[n]
    end

    for hv in h.hg_tail.he2v
        if v < n && haskey(hv, n)
            hv[v] = hv[n]
            delete!(hv, n)
        else
            delete!(hv, v)
        end
    end

    for hv in h.hg_head.he2v
        if v < n && haskey(hv, n)
            hv[v] = hv[n]
            delete!(hv, n)
        else
            delete!(hv, v)
        end
    end

    resize!(h.v_meta, n - 1)

    resize!(h.hg_tail.v2he, n - 1)
    resize!(h.hg_tail.v_meta, n - 1)
    resize!(h.hg_head.v2he, n - 1)
    resize!(h.hg_head.v_meta, n - 1)

    return h
end


"""
    add_hyperedge!(
	h::DirectedHypergraph{T, V, E, D};
        vertices_tail::D = D(),
	vertices_head::D = D(),
        he_meta_tail::Union{E,Nothing}=nothing,
	he_meta_head::Union{E,Nothing}=nothing
    ) where {T <: Real, V, E, D <: AbstractDict{Int,T}}

Adds a dihyperedge to a given dihypergraph `h`. Optionally, existing vertices can be added to the
created dihyperedge in the tail and/or head sides. The paramater `vertices_tail` represents a
dictionary of vertex indices and weights for the tail of the dihyperedge; `vertices_head`
represents the vertex indices and values for the head of the dihyperedge. Additionally, metadata
can be stored with the hyperedge using the `he_meta_tail` and `he_meta_head` keyword arguments.

**Arguments**

* `h` : Dihypergraph
* `vertices_tail` : Vertex index-weight dictionary with information about the vertices in the new
    dihyperedge's tail; default is an empty dict of type `D`
* `vertices_head` : Vertex index-weight dictionary with information about the vertices in the new
    dihyperedge's head; default is an empty dict of type `D`
* `he_meta_tail` : (Optional) metadata for the new dihyperedge's tail; default is `nothing`
* `he_meta_head` : (Optional) metadata for the new dihyperedge's head; default is `nothing`

**Returns**

`ix`, the index of the new dihyperedge

"""
function SimpleHypergraphs.add_hyperedge!(
        h::DirectedHypergraph{T, V, E, D};
        vertices_tail::D = D(),
        vertices_head::D = D(),
        he_meta_tail::Union{E, Nothing} = nothing,
        he_meta_head::Union{E, Nothing} = nothing
    ) where {T <: Real, V, E, D <: AbstractDict{Int, T}}
    @boundscheck (checkbounds(h.hg_tail, k, 1) for k in keys(vertices_tail))
    @boundscheck (checkbounds(h.hg_head, 1, k) for k in keys(vertices_head))

    push!(h.hg_tail.he2v, vertices_tail)
    push!(h.hg_head.he2v, vertices_head)

    # Should always be identical to length(h.hg_head.he2v)
    ix = length(h.hg_tail.he2v)
    for k in keys(vertices_tail)
        h[1, k, ix] = vertices_tail[k]
    end
    for k in keys(vertices_head)
        h[2, k, ix] = vertices_head[k]
    end
    push!(h.he_meta_tail, he_meta_tail)
    push!(h.he_meta_head, he_meta_head)
    return ix
end


"""
    remove_hyperedge!(h::DirectedHypergraph, e::Int)

Removes the dihyperedge `e` from a given dihypergraph `h`. Note that running this function will
cause reordering of hyperedges in the dihypergraph: the dihyperedge `e` will replaced by the last
dihyperedge of the dihypergraph and the list of dihyperedges (and dihyperedge metadata) will be shrunk.

**Arguments**

* `h` : Dihypergraph
* `e` : Dihyperedge index

**Returns**

`h`

"""
function SimpleHypergraphs.remove_hyperedge!(h::DirectedHypergraph, e::Int)
    ne = nhe(h)
    @assert(e <= ne)
    if e < ne
        h.he_meta_tail[e] = h.he_meta_tail[ne]
        h.he_meta_head[e] = h.he_meta_head[ne]

        h.hg_tail.he2v[e] = h.hg_tail.he2v[ne]
        h.hg_head.he2v[e] = h.hg_head.he2v[ne]
    end

    for he in h.hg_tail.v2he
        if e < ne && haskey(he, ne)
            he[e] = he[ne]
            delete!(he, ne)
        else
            delete!(he, e)
        end
    end

    for he in h.hg_head.v2he
        if e < ne && haskey(he, ne)
            he[e] = he[ne]
            delete!(he, ne)
        else
            delete!(he, e)
        end
    end

    resize!(h.he_meta_tail, ne - 1)
    resize!(h.he_meta_head, ne - 1)

    resize!(h.hg_tail.he2v, ne - 1)
    resize!(h.hg_head.he2v, ne - 1)
    resize!(h.hg_tail.he_meta, ne - 1)
    resize!(h.hg_head.he_meta, ne - 1)

    return h
end


"""
    prune_hypergraph!(h::H) where {H <: AbstractDirectedHypergraph}

Remove all vertices with degree 0 and all hyperedges of size 0.

**Arguments**

* `h` : Dihypergraph

**Returns**

`h`

"""
function SimpleHypergraphs.prune_hypergraph!(h::H) where {H <: AbstractDirectedHypergraph}
    for e in reverse(1:nhe(h))
        length(h.hg_tail.he2v[e]) == 0 && length(h.hg_head.he2v[e]) == 0 && SimpleHypergraphs.remove_hyperedge!(h, e)
    end
    for v in reverse(1:nhv(h))
        length(h.hg_tail.v2he[v]) == 0 && length(h.hg_tail.v2he[v]) == 0 && SimpleHypergraphs.remove_vertex!(h, v)
    end
    return h
end

"""
    prune_hypergraph(h::H) where {H <: AbstractDirectedHypergraph}

Remove all vertices with degree 0 and all hyperedges of size 0.

**Arguments**

* `h` : Dihypergraph

**Returns**

A deepcopy of `h` with all isolated vertices and empty dihyperedges removed

"""
function SimpleHypergraphs.prune_hypergraph(h::H) where {H <: AbstractDirectedHypergraph}
    return prune_hypergraph!(deepcopy(h))
end

"""
    set_vertex_meta!(
	h::DirectedHypergraph{T, V, E, D},
	new_value::Union{V,Nothing},
        v::Int
    ) where {T <: Real, V, E, D <: AbstractDict{Int,T}}

Sets a new metadata value `new_value` for the vertex `v` in the dihypergraph `h`.

**Arguments**

* `h` : Dihypergraph
* `new_value` : Metadata to replace existing metadata
* `v` : Vertex index 

**Returns**

`h.v_meta`

"""
function SimpleHypergraphs.set_vertex_meta!(
        h::DirectedHypergraph{T, V, E, D},
        new_value::Union{V, Nothing},
        v::Int
    ) where {T <: Real, V, E, D <: AbstractDict{Int, T}}
    checkbounds(h.v_meta, v)
    h.v_meta[v] = new_value
    return h.v_meta
end


"""
    get_vertex_meta(
	h::DirectedHypergraph{T, V, E, D},
	v::Int
    ) where {T <: Real, V, E, D <: AbstractDict{Int,T}}

Returns a metadata value associated with the vertex with index `v` in the dihypergraph `h`.

**Arguments**

* `h` : Dihypergraph
* `v` : Vertex index

**Returns**

Vertex metadata at index `v`

"""
function SimpleHypergraphs.get_vertex_meta(
        h::DirectedHypergraph{T, V, E, D},
        v::Int
    ) where {T <: Real, V, E, D <: AbstractDict{Int, T}}
    @boundscheck checkbounds(h.v_meta, v)
    return h.v_meta[v]
end


"""
    set_hyperedge_meta!(
	h::DirectedHypergraph{T, V, E, D},
        new_value_tail::Union{E,Nothing},
	new_value_head::Union{E,Nothing},
	e::Int
    ) where {T <: Real, V, E, D <: AbstractDict{Int,T}}

Sets new metadata for the hyperedge with index `e` in the dihypergraph `h`.

**Arguments**

* `h` : Dihypergraph
* `new_value_tail` : Metadata value for the dihyperedge's tail (or `nothing`)
* `new_value_head` : Metadata value for the dihyperedge's head (or `nothing`)
* `e` : Index for the dihyperedge

**Returns**

`(h.he_meta_tail, h.he_meta_head)`

"""
function SimpleHypergraphs.set_hyperedge_meta!(
        h::DirectedHypergraph{T, V, E, D},
        new_value_tail::Union{E, Nothing},
        new_value_head::Union{E, Nothing},
        e::Int
    ) where {T <: Real, V, E, D <: AbstractDict{Int, T}}
    @boundscheck checkbounds(h.he_meta_tail, e)
    @boundscheck checkbounds(h.he_meta_head, e)

    h.he_meta_tail[e] = new_value_tail
    h.he_meta_head[e] = new_value_head

    return (h.he_meta_tail, h.he_meta_head)
end


"""
    set_hyperedge_meta!(
	h::DirectedHypergraph{T, V, E, D},
        new_value::Union{E,Nothing},
	e::Int,
	side::HyperedgeDirection
    ) where {T <: Real, V, E, D <: AbstractDict{Int,T}}

Sets a new meta value `new_value` for the hyperedge with index `e` in the direction `side` for a
dihypergraph `h`.

**Arguments**

* `h` : Dihypergraph
* `new_value` : Metadata value for this dihyperedge (or `nothing`)
* `e` : Dihyperedge index
* `side` : a `HyperedgeDirection` (either 1 or 2)

**Returns**

`h.he_meta_tail` (if `side == tail`) or `h.he_meta_head` (if `side == head`)

"""
function SimpleHypergraphs.set_hyperedge_meta!(
        h::DirectedHypergraph{T, V, E, D},
        new_value::Union{E, Nothing},
        e::Int,
        side::HyperedgeDirection
    ) where {T <: Real, V, E, D <: AbstractDict{Int, T}}

    return if side == tail
        checkbounds(h.he_meta_tail, id)
        h.he_meta_tail[id] = new_value
        h.he_meta_tail
    else
        checkbounds(h.he_meta_head, id)
        h.he_meta_head[id] = new_value
        h.he_meta_head
    end

end


"""
    get_hyperedge_meta(
	h::DirectedHypergraph{T, V, E, D},
	e::Int
    ) where {T <: Real, V, E, D <: AbstractDict{Int,T}}

Returns the metadata for the dihyperedge with index `e` in the dihypergraph `h`.

**Arguments**

* `h` : Dihypergraph
* `e` : Dihyperedge index

**Returns**

A 2-tuple `(tail_meta, head_meta)`, where `tail_meta` is the metadata associated with dihyperedge
`e` on the tail side, and `head_meta` is the metadata associated with dihyperedge `e` on the head
side.

"""
function SimpleHypergraphs.get_hyperedge_meta(
        h::DirectedHypergraph{T, V, E, D}, id::Int
    ) where {T <: Real, V, E, D <: AbstractDict{Int, T}}
    checkbounds(h.he_meta_tail, id)
    checkbounds(h.he_meta_head, id)

    return (h.he_meta_tail[id], h.he_meta_head[id])
end


"""
    get_hyperedge_meta(
	h::DirectedHypergraph{T, V, E, D},
	e::Int,
	side::HyperedgeDirection
    ) where {T <: Real, V, E, D <: AbstractDict{Int,T}}

Returns a metadata value for one side of the dihyperedge with index `e` in the dihypergraph `h`.

**Arguments**

* `h` : Dihypergraph
* `e` : Dihyperedge index
* `side` : a `HyperedgeDirection` (either 1 or 2)

**Returns**

`meta`, where `meta` is the metadata associated with the tail (if `side == tail`) or the head
(if `side == head`) of the dihyperedge with index `e`

"""
function SimpleHypergraphs.get_hyperedge_meta(
        h::DirectedHypergraph{T, V, E, D},
        e::Int,
        side::HyperedgeDirection
    ) where {T <: Real, V, E, D <: AbstractDict{Int, T}}

    return if side == tail
        checkbounds(h.he_meta_tail, id)
        h.he_meta_tail[id]
    else
        checkbounds(h.he_meta_head, id)
        h.he_meta_head[id]
    end
end


"""
    nhe(h::H) where {H <: AbstractDirectedHypergraph}

Return the number of dihyperedges in the dihypergraph `h`.

**Arguments**

* `h` : Dihypergraph

**Returns**

The number of dihyperedges in `h`. Note that if, for some reason, the `he2v` lists in `h.hg_tail`
and `h.hg_head` have different lengths (implying that there are different numbers of dihyperedge
tails than heads), an error is raised.

"""
function SimpleHypergraphs.nhe(h::H) where {H <: AbstractDirectedHypergraph}
    return (length(h.hg_tail.he2v) == length(h.hg_head.he2v)) ? length(h.hg_tail.he2v) : throw("Tail and head sides of hypergraph have different numbers of hyperedges!")
end


"""
    nhv(h::H) where {H <: AbstractDirectedHypergraph}

Return the number of vertices in the dihypergraph `h`.

**Arguments**

* `h` : Dihypergraph

**Returns**

The number of vertices in `h`
"""
function SimpleHypergraphs.nhv(h::H) where {H <: AbstractDirectedHypergraph}
    return (length(h.hg_tail.v2he) == length(h.hg_head.v2he)) ? length(h.hg_tail.v2he) : throw("Tail and head sides of hypergraph have different numbers of hyperedges!")
end


"""
    is_loopless(h::H; strict::Bool = true) where {H <: AbstractDirectedHypergraph}

We define a *loop* in two ways: "strict" and "loose".

In the "strict" definition, a dihyperedge is a loop if its tail and head contain exactly the same
vertices (weights can be different). In the "loose" definition, a dihyperedge is a loop if there is
any vertex that appears in both the tail and the head.

A dihypergraph is "loopless" if there are no dihyperedges that are loops.

**Arguments**

* `h` : Dihypergraph
* `strict` : Bool

**Returns**

`loopless`, a `Bool` indicating if `h` is loopless (`true`) or not (`false`)

"""
function is_loopless(h::H; strict::Bool = true) where {H <: AbstractDirectedHypergraph}
    loopless = true
    for e in 1:nhe(h)
        if strict
            if keys(h.hg_tail.he2v[e]) == keys(h.hg_head.he2v[e])
                loopless = false
                break
            end
        else
            if length(intersect(keys(h.hg_tail.he2v[e]), keys(h.hg_head.he2v[e]))) > 0
                loopless = false
                break
            end
        end
    end

    return loopless
end

"""
    is_simple(h::H) where {H <: AbstractDirectedHypergraph}

A dihypergraph is *simple* if it is loopless (using the strict definition; see `is_loopless`)
and if it is without repeated dihyperedge (meaning that there are no two hyperedges in the
dihypergraph with identical tails and heads, ignoring the vertex-hyperedge weights).

**Arguments**

* `h` : Dihypergraph

**Returns**

A `Bool` indicating if `h` is simple (`true`) or not (`false`)

"""
function is_simple(h::H) where {H <: AbstractDirectedHypergraph}
    # Is the dihypergraph without loop?
    if !is_loopless(h)
        return false
    end

    # Is the dihypergraph without repeated hyperedge?
    he_vertices = [(keys(h.hg_tail.he2v[e]), keys(h.hg_head.he2v[e])) for e in 1:nhe(h)]
    if any(x -> x > 1, values(countmap(he_vertices)))
        return false
    end

    return true
end


"""
    is_b_hypergraph(h::H) where {H <: AbstractDirectedHypergraph}

A *B-edge* is a dihyperedge which may have multiple vertices in the tail but which has exactly one
vertex in the head. A *B-hypergraph* is a dihypergraph where all dihyperedges are B-edges.

**Arguments**

* `h` : Dihypergraph

**Returns**

A `Bool` indicating if `h` is a B-hypergraph (`true`) or not (`false`)

"""
function is_b_hypergraph(h::H) where {H <: AbstractDirectedHypergraph}
    return all(x -> length(h.hg_head.he2v[x]) == 1, 1:nhe(h))
end


"""
    is_f_hypergraph(h::H) where {H <: AbstractDirectedHypergraph}

An *F-edge* is a dihyperedge which may have multiple vertices in the head but which has exactly one
vertex in the tail. An *F-hypergraph* is a dihypergraph where all dihyperedges are F-edges.

**Arguments**

* `h` : Dihypergraph

**Returns**

A `Bool` indicating if `h` is an F-hypergraph (`true`) or not (`false`)

"""
function is_f_hypergraph(h::H) where {H <: AbstractDirectedHypergraph}
    return all(x -> length(h.hg_tail.he2v[x]) == 1, 1:nhe(h))
end


"""
    is_bf_hypergraph(h::H) where {H <: AbstractDirectedHypergraph}

A *BF-hypergraph* is a dihypergraph where all dihyperedges are either B-edges (see
`is_b_hypergraph`), meaning that they have exactly one vertex in the head, or F-edges (see
`is_f_hypergraph`), meaning that they have exactly one vertex in the tail.

**Arguments**

* `h` : Dihypergraph

**Returns**

A `Bool` indicating if `h` is a BF-hypergraph (`true`) or not (`false`)

"""
function is_bf_hypergraph(h::H) where {H <: AbstractDirectedHypergraph}
    return all(x -> length(h.hg_tail.he2v[x]) == 1 || length(h.hg_head.he2v[x]) == 1, 1:nhe(h))
end
