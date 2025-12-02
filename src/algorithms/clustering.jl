_inc(hg::H, i::Int, α::Int, side::Int) where {H <: AbstractDirectedHypergraph} = isnothing(hg[i, α][side]) ? 0 : 1

function _num_quads(hg::H, i::Int) where {H <: AbstractDirectedHypergraph}
    quads = 0
    nv = nhv(hg)
    ne = nhe(hg)

    # TODO: there must be a better implementation
    for α in 1:ne
        for β in α+1:ne
            for j in 1:nv
                if i == j
                    continue
                end

                quads += (
                    (_inc(hg, i, α, 1) + _inc(hg, i, α, 2))
                    * (_inc(hg, i, β, 1) + _inc(hg, i, β, 2))
                    * (_inc(hg, j, α, 1) + _inc(hg, j, α, 2))
                    * (_inc(hg, j, β, 1) + _inc(hg, j, β, 2))
                )
            end
        end
    end
    
    quads
end

function _max_num_quads(
    hg::H,
    i::Int,
    tail_deg_exc::AbstractVector{Int},
    head_deg_exc::AbstractVector{Int}
) where {H <: AbstractDirectedHypergraph}
    ne = nhe(hg)

    # TODO: there must be a better implementation
    qmax = 0
    for α in 1:ne
        for β in α+1:ne
            inc_ab = (_inc(hg, i, α, 1) + _inc(hg, i, α, 2)) * (_inc(hg, i, β, 1) + _inc(hg, i, β, 2))

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
                    W = 2 * degs_sorted[1] +  2 * degs_sorted[2]
                else
                    W =  2 * degs_sorted[1] + degs_sorted[2] + degs_sorted[3]
                end
            end
            qmax += inc_ab * W
        end
    end
    qmax
end

"""
    SimpleHypergraphs.quad_clustering_coefficient(hg::H, i::Int) where {H <: AbstractDirectedHypergraph} 

    SimpleHypergraphs.quad_clustering_coefficient(hg::H) where {H <: AbstractDirectedHypergraph}

    Implements the "quad clustering coefficient" (QCC) for directed hypergraphs, as described in:
    Ha, Neri, and Annibale, Chaos 34, 043102 (2024), DOI: 10.1063/5.0188246

    A *quad* is the shortest simple cycle in a hypergraph, consisting of two vertices `i` and `j` that are both
    incident on the same two hyperedges `α` and `β`. The QCC is a density, describing the fraction of all possible
    "quads" a particular vertex `i` participates in. It is always true that `0 <= QCC(hg, i) <= 1`.
"""
function SimpleHypergraphs.quad_clustering_coefficient(hg::H, i::Int) where {H <: AbstractDirectedHypergraph}
    # Degrees of hyperedge tails and heads, not including vertex `i`
    tail_deg_exc = zeros(Int, nhe(hg))
    head_deg_exc = zeros(Int, nhe(hg))
    for he in 1:nhe(hg)
        tail_deg_exc[he] = length([k for k in keys(hg.hg_tail.he2v[he]) if k != i])
        tail_deg_exc[he] = length([k for k in keys(hg.hg_head.he2v[he]) if k != i])
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

    return q / qmax
end

function SimpleHypergraphs.quad_clustering_coefficient(hg::H) where {H <: AbstractDirectedHypergraph}
    return [SimpleHypergraphs.quad_clustering_coefficient(hg, i) for i in 1:nhv(hg)]
end