"""
    _inc(h::H, v::Int, e::Int, side::Int) where {H <: AbstractDirectedHypergraph}    

Calculate a binary incidence relation. If vertex `v` is included in side `side` (`1` for "tail",
`2` for "head") of dihyperedge `e` in dihypergraph `h`, then this returns 1; otherwise, it returns
`0`.

**Arguments**

* `h` : Dihypergraph
* `v` : Vertex index
* `e` : Dihyperedge index
* `side` : Either `1` ("tail") or `2` ("head")

**Returns**

`0` for non-incidence, or else `1` for incidence

"""
_inc(h::H, v::Int, e::Int, side::Int) where {H <: AbstractDirectedHypergraph} = isnothing(h[v, e][side]) ? 0 : 1

"""
    _num_quads(h::H, v::Int) where {H <: AbstractDirectedHypergraph}

Calculate the number of quads involving vertex `v` in dihypergraph `h`

**Arguments**

* `h` : Dihypergraph
* `v` : Vertex index

**Returns**

`quads`, the number of quads involving vertex `v`

"""
function _num_quads(h::H, v::Int) where {H <: AbstractDirectedHypergraph}
    quads = 0

    nv = nhv(h)

    rel_hes = sort!(collect(union(keys(h.hg_tail.v2he[v]), keys(h.hg_head.v2he[v]))))

    for α in 1:length(rel_hes)
        a = rel_hes[α]
        for β in (α + 1):length(rel_hes)
            b = rel_hes[β]

            for j in 1:nv
                if v == j
                    continue
                end

                if !(
                        (isnothing(h[j, a][1]) && isnothing(h[j, a][2]))
                            || (isnothing(h[j, b][1]) && isnothing(h[j, b][2]))
                    )
                    quads += 1
                end
            end
        end
    end

    return quads

end

"""
    _max_num_quads(
        h::H,
        v::Int,
        tail_deg_exc::AbstractVector{Int},
        head_deg_exc::AbstractVector{Int}
    ) where {H <: AbstractDirectedHypergraph}

Calculate the maximum number of quads vertex `v` could participate in in dihypergraph `h`.

**Arguments**

* `h` : Dihypergraph
* `v` : Vertex index
* `tail_deg_exc` : A vector containing the tail degrees of every dihyperedge in `h`, excluding `v`
* `head_deg_exc` : A vector containing the head degrees of every dihyperedge in `h`, excluding `v`

**Returns**

`qmax`, the maximum possible number of quads involving `v` in `h`

"""
function _max_num_quads(
        h::H,
        v::Int,
        tail_deg_exc::AbstractVector{Int},
        head_deg_exc::AbstractVector{Int}
    ) where {H <: AbstractDirectedHypergraph}

    ne = nhe(h)

    # TODO: there must be a better implementation
    qmax = 0
    for α in 1:ne
        for β in (α + 1):ne
            inc_ab = (_inc(h, v, α, 1) + _inc(h, v, α, 2)) * (_inc(h, v, β, 1) + _inc(h, v, β, 2))

            if inc_ab == 0
                continue
            end

            if tail_deg_exc[α] == head_deg_exc[α] && tail_deg_exc[β] == head_deg_exc[β]
                W = 4 * min(tail_deg_exc[α], tail_deg_exc[β])
            elseif tail_deg_exc[α] == head_deg_exc[α] && tail_deg_exc[β] != head_deg_exc[β]
                degs = [tail_deg_exc[α], tail_deg_exc[β], head_deg_exc[β]]
                min_deg = minimum(degs)
                max_deg = maximum(degs)
                med_deg = Statistics.median(degs)

                if min_deg == tail_deg_exc[α]
                    W = 4 * min_deg
                elseif length(Set(degs)) == 3
                    W = 2 * min_deg + 2 * med_deg
                else
                    W = 2 * min_deg + 2 * max_deg
                end
            elseif tail_deg_exc[α] != head_deg_exc[α] && tail_deg_exc[β] == head_deg_exc[β]
                degs = [tail_deg_exc[α], head_deg_exc[α], tail_deg_exc[β]]
                min_deg = minimum(degs)
                max_deg = maximum(degs)
                med_deg = Statistics.median(degs)

                if min_deg == tail_deg_exc[β]
                    W = 4 * min_deg
                elseif length(Set(degs)) == 3
                    W = 2 * min_deg + 2 * med_deg
                else
                    W = 2 * min_deg + 2 * max_deg
                end
            elseif length(Set([tail_deg_exc[α], tail_deg_exc[β], head_deg_exc[α], head_deg_exc[β]])) != 4
                degs = [tail_deg_exc[α], tail_deg_exc[β], head_deg_exc[α], head_deg_exc[β]]
                setdeg = Set(degs)
                min_deg = minimum(degs)
                max_deg = maximum(degs)
                med_deg = Statistics.median(setdeg)

                if length(setdeg) == 2
                    W = 3 * min_deg + max_deg
                elseif min(degs[1], degs[3]) == min(degs[2], degs[4])
                    W = 3 * min_deg + med_deg
                elseif max(degs[1], degs[3]) == min(degs[2], degs[4]) || min(degs[1], degs[3]) == max(degs[2], degs[4])
                    W = 2 * min_deg + 2 * med_deg
                else
                    W = 2 * min_deg + med_deg + max_deg
                end
            else
                degs = [tail_deg_exc[α], tail_deg_exc[β], head_deg_exc[α], head_deg_exc[β]]
                degs_sorted = sort(degs)

                if max(degs[1], degs[3]) < min(degs[2], degs[4]) || min(degs[1], degs[3]) > max(degs[2], degs[4])
                    W = 2 * degs_sorted[1] + 2 * degs_sorted[2]
                else
                    W = 2 * degs_sorted[1] + degs_sorted[2] + degs_sorted[3]
                end
            end
            qmax += inc_ab * W
        end
    end
    return qmax
end

"""
    quad_clustering_coefficient(hg::H, v::Int) where {H <: AbstractDirectedHypergraph} 

    quad_clustering_coefficient(hg::H) where {H <: AbstractDirectedHypergraph}

Implements the "quad clustering coefficient" (QCC) for directed hypergraphs, as described in:
Ha, Neri, and Annibale, Chaos 34, 043102 (2024), DOI: 10.1063/5.0188246

A *quad* is the shortest simple cycle in a hypergraph, consisting of two vertices `i` and `j` that
are both incident on the same two hyperedges `α` and `β`. The QCC is a density, describing the
fraction of all possible "quads" a particular vertex `i` participates in. It is always true that
`0 <= QCC(hg, i) <= 1`.

**Arguments**

* `h` : Dihypergraph
* `v` : Vertex index

**Returns**

If a vertex index is given, then this function returns the QCC of that vertex. If no index is
given, then instead a vector with the QCCs of all vertices (in index order) is returned.

"""
function SimpleHypergraphs.quad_clustering_coefficient(hg::H, i::Int) where {H <: AbstractDirectedHypergraph}
    # Degrees of hyperedge tails and heads, not including vertex `i`
    tail_deg_exc = zeros(Int, nhe(hg))
    head_deg_exc = zeros(Int, nhe(hg))
    for he in 1:nhe(hg)
        tail_deg_exc[he] = length([k for k in keys(hg.hg_tail.he2v[he]) if k != i])
        head_deg_exc[he] = length([k for k in keys(hg.hg_head.he2v[he]) if k != i])
    end

    # Can this vertex participate in any quads, based on its connectivity?
    degree_thresh = 0
    for he in union(Set(keys(hg.hg_tail.v2he[i])), Set(keys(hg.hg_head.v2he[i])))
        degree_thresh += tail_deg_exc[he] + head_deg_exc[he]
    end

    if degree_thresh < 2
        return 0.0
    end

    q = _num_quads(hg, i)
    qmax = _max_num_quads(hg, i, tail_deg_exc, head_deg_exc)

    if qmax == 0
        return 0.0
    end

    return q / qmax
end

function SimpleHypergraphs.quad_clustering_coefficient(hg::H) where {H <: AbstractDirectedHypergraph}
    return [quad_clustering_coefficient(hg, v) for v in 1:nhv(hg)]
end
