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

function _max_num_quads(hg::H, i::Int) where {H <: AbstractDirectedHypergraph}
    ne = nhe(hg)
    tail_degrees = length.(hg.hg_tail.he2v)
    head_degrees = length.(hg.hg_head.he2v)
    # TODO: there must be a better implementation
    qmax = 0
    for α in 1:ne
        for β in α+1:ne
            inc_ab = (_inc(hg, i, α, 1) + _inc(hg, i, α, 2)) * (_inc(hg, i, β, 1) + _inc(hg, i, β, 2))

            if inc_ab == 0
                continue
            end

            if tail_degrees[α] == tail_degrees[β] && head_degrees[α] == head_degrees[β]
                W = 4 * (min(tail_degrees[α], head_degrees[α]) - 1)
            elseif tail_degrees[α] == tail_degrees[β] || head_degrees[α] == head_degrees[β]
                degs = [tail_degrees[α], tail_degrees[β], head_degrees[β]]
                min_deg = minimum(degs)
                max_deg = maximum(degs)
                med_deg = Statistics.median(degs)

                if min_deg == tail_degrees[α]
                    W = 4 * (min_deg - 1)
                elseif min_deg != tail_degrees[α] && length(Set(degs)) == 3
                    W = 2 * (min_deg - 1) + 2 * (med_deg - 1)
                else
                    W = 2 * (min_deg - 1) + 2 * (max_deg - 1)
                end
            elseif length(Set([tail_degrees[α], tail_degrees[β], head_degrees[α], head_degrees[β]])) != 4
                degs = [tail_degrees[α], tail_degrees[β], head_degrees[α], head_degrees[β]]
                setdeg = Set(degs)
                min_deg = minimum(degs)
                max_deg = maximum(degs)
                med_deg = Statistics.median(setdeg)

                if length(setdeg) == 2
                    W = 3 * (min_deg - 1) + (max_deg - 1)
                elseif min(degs[1], degs[3]) == min(degs[2], degs[4])
                    W = 3 * (min_deg - 1) + (med_deg - 1)
                elseif max(degs[1], degs[3]) == min(degs[2], degs[4]) || min(degs[1], degs[3]) == max(degs[2], degs[4])
                    W = 2 * (min_deg - 1) + 2 * (med_deg - 1)
                else
                    W = 2 * (min_deg - 1) + (med_deg - 1) + (max_deg - 1)
                end
            else
                degs = [tail_degrees[α], tail_degrees[β], head_degrees[α], head_degrees[β]]
                degs_sorted = sort(degs)

                if max(degs[1], degs[3]) < min(degs[2], degs[4]) || min(degs[1], degs[3]) > max(degs[2], degs[4])
                    W = 2 * (degs_sorted[1] - 1) +  2 * (degs_sorted[2] - 1)
                else
                    W =  2 * (degs_sorted[1] - 1) + (degs_sorted[2] - 1) + (degs_sorted[3] - 1)
                end
            end
            qmax += inc_ab * W
        end
    end
    qmax
end

"""

"""
function SimpleHypergraphs.quad_clustering_coefficient(hg::H, i::Int) where {H <: AbstractDirectedHypergraph}
end

"""

"""
function SimpleHypergraphs.quad_clustering_coefficient(hg::H) where {H <: AbstractDirectedHypergraph}
end