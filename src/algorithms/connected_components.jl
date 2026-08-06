"""
    _default_heselect(h::H, v::Int; reverse::Bool = false) where {H <: AbstractDirectedHypergraph}

Default hyperedge selection algorithm (for `random_walk`). All hyperedges are given equal weights.

**Arguments**

* `h` : Dihypergraph
* `v` : Vertex index
* `reverse` : If `true` (default `false`), consider hyperedges where `v` is in the head, rather
    than in the tail

**Returns**

* A sorted list of the possible hyperedge indices
* A weight vector with the same length as the list of indices (all entries are `1`)

"""
function _default_heselect(h::H, v::Int; reverse::Bool = false) where {H <: AbstractDirectedHypergraph}
    he_tail, he_head = gethyperedges(h, v)

    if reverse
        hes = he_head
    else
        hes = he_tail
    end

    return sort!(collect(keys(hes))), ones(length(hes))
end

"""
    _default_vselect(h::H, e::Int; reverse::Bool = false) where {H <: AbstractDirectedHypergraph}

Default vertex selection algorithm (for `random_walk`). All vertices are given equal weights.

**Arguments**

* `h` : Dihypergraph
* `e` : Vertex index
* `reverse` : If `true` (default `false`), consider vertices in the tail of dihyperedge `e`, rather
    than in the head

**Returns**

* A sorted list of the possible vertex indices
* A weight vector with the same length as the list of indices (all entries are `1`)

"""
function _default_vselect(h::H, he::Int; reverse::Bool = false) where {H <: AbstractDirectedHypergraph}
    vs_tail, vs_head = getvertices(h, he)

    if reverse
        vs = vs_tail
    else
        vs = vs_head
    end

    return sort!(collect(keys(vs))), ones(length(vs))

end


"""
    random_walk(
        h::H,
        start::Int;
        heselect::Function,
        vselect::Function,
        reverse::bool
    ) where {H <: AbstractDirectedHypergraph}

Return the next vertex visited in a random walk starting from vertex `start`. First, a hyperedge
is sampled with weights based on the output of the `heselect` function (by default, each hyperedge
is sampled with the same probability). Next, a vertex within that dihyperedge is selected, with
weights based on the output of the `vselect` function (by default, each vertex, including the
source, is sampled with the same probability).

`heselect` and `vselect` functions take two arguments: a dihypergraph and, respectively, a vertex
index or a hyperedge index. The return values of both functions should be, respectively, a list of
dihyperedge or vertex indices and their weights.

**Arguments**

* `h` : Dihypergraph
* `start` : An initial vertex index
* `heselect` : A function (see above) that assigns weights to valid dihyperedges; default assigns
    uniform weights to all dihyperedges
* `vselect` : A function (see above) that assigns weights to valid vertices; default assigns
    uniform weights to all vertices
* `reverse` : If `true` (default `false`), the walk will go along dihyperedges from head to tail,
    rather than from tail to head

**Returns**

A selected vertex index

"""
function SimpleHypergraphs.random_walk(
        h::H, start::Int;
        heselect::Function = _default_heselect,
        vselect::Function = _default_vselect,
        reverse::Bool = false
    ) where {H <: AbstractDirectedHypergraph}
    1 <= start <= nhv(h) || throw(ArgumentError("invalid start vertex index"))
    hes, hew = heselect(h, start, reverse = reverse)
    he = sample(hes, Weights(hew))
    ves, vw = vselect(h, he, reverse = reverse)
    return sample(ves, Weights(vw))
end


"""
    get_weakly_connected_components(h::H) where {H <: AbstractDirectedHypergraph}

Return an array of weakly connected components in the dihypergraph `h` by first converting the
dihypergraph into an undirected hypergraph and then obtaining the conected components of
that hypergraph.

**Arguments**

* `h` : Dihypergraph

**Returns**

A vector where each member is itself a vector consisting of the indices of one weakly connected
component

"""
function get_weakly_connected_components(h::H) where {H <: AbstractDirectedHypergraph}
    undirected = to_undirected(h)
    return get_connected_components(undirected)
end


"""
    _visit!(h::H, v::Int) where {H <: AbstractDirectedHypergraph}

Determines the B-connected component of a vertex `v` in directed hypergraph `h`.
This is an auxiliary function for `get_strongly_connected_components`, which
determines the strongly connected components of a directed hypergraph.

For details, see Francisco José Martín-Recuerda Moyano, "Strong Connectivity in Directed
Hypergraphs and its Application to the Atomic Decomposition of Ontologies", PhD dissertation, 2016.

**Arguments**

* `h` : Dihypergraph
* `v` : Vertex index

**Returns**

The B-connected component containing `v` (a `Set{Int}` containing vertex indices)

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

    return bcc
end


"""
    get_strongly_connected_components(h::H) where {H <: AbstractDirectedHypergraph}

Return an array of strongly connected components in the dihypergraph `h`, based on the "naive"
algorithm of Francisco José Martín-Recuerda Moyano (PhD dissertation, 2016).

**Arguments**

* `h` : Dihypergraph

**Returns**

A vector, where each member of the vector is itself a (sorted) vector containing the vertex indices
of one strongly connected component of `h`

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

    return [sort!(collect(v)) for (_, v) in T if length(v) != 0]
end


"""
    is_connected(h::H) where {H <: AbstractDirectedHypergraph}

A directed hypergraph is *connected* if it has exactly one connected component. See
`get_weakly_connected_components`.

**Arguments**

* `h` : Dihypergraph

**Returns**

A `Bool`; `true` if `h` is (weakly) connected, and `false` otherwise

"""
function is_connected(h::H) where {H <: AbstractDirectedHypergraph}
    return length(get_weakly_connected_components(h)) == 1
end

"""
    is_strongly_connected(h::H) where {H <: AbsctractDirectedHypergraphs}

A directed hypergraph is *strongly connected* if it has exactly one strongly connected component.
See `get_strongly_connected_components`.

**Arguments**

* `h` : Dihypergraph

**Returns**

A `Bool`; `true` if `h` is strongly connected, and `false` otherwise

"""
function is_strongly_connected(h::H) where {H <: AbstractDirectedHypergraph}
    return length(get_strongly_connected_components(h)) == 1
end
