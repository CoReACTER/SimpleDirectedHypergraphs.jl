# TODO: weights (default vector of 1.0 of length nve(hg))

struct DiHyperPathState{T} where {T<:Real}
    reached_vs::BitVector
    marked_hes::BitVector
    removed_hes::BitVector
    hes_tail_count::Vector{Int}
    he_inedges::Vector{Set{Int}}
    edge_weights::Vector{T}
    edge_costs::Vector{T}
    edge_heap_points::Vector{Union{Nothing, Int}}
end

function initialize_dihyperpath_state(
        hg::H,
        hyperedge_weights::Vector{T}
    ) where {H <: AbstractDirectedHypergraph, T <: Real}
    nv = nhv(hg)
    ne = nhe(hg)

    edge_heap_points = Vector{Union{Nothing, Int}}(undef, ne)
    edge_heap_points .= nothing

    return DiHyperPathState{T}(
        BitVector(zeros(nv)),
        BitVector(zeros(ne)),
        BitVector(zeros(ne)),
        length.(keys.(hg.hg_tail.he2v)),
        [Set{Int}() for _ in 1:ne],
        hyperedge_weights,
        fill(typemax(T), ne),
        edge_heap_points
    )
end

"""
Kind of silly function
"""
function validate_weights(w::AbstractVector{T}) where {T <: Real}
    @assert all(w .>= 0)
end

"""

"""
function forward_reachable(hg::H, source::Int, state::DiHyperPathState{T}) where {H <: AbstractDirectedHypergraph, T <: Real}
    # Priority queue of reached vertices
    Q = Queue{Int}()
    enqueue!(Q, source)
    
    state.reached_vs[source] = true

    # Which vertices/hyperedges have been reached?
    reached_vs = Set{Int}()
    reached_es = Set{Int}()

    while length(Q) > 0
        v = popfirst!(Q)
        push!(reached_vs, v)

        for out_e in keys(hg.hg_tail.v2he[v])
            # Following pseudocode exactly. This feels awkward; how slow would it be to just query reached_vs?
            state.hes_tail_count[out_e] -= 1

            if state.hes_tail_count[out_e] == 0
                push!(reached_es, out_e)

                for w in keys(hg.hg_head.he2v[out_e])
                    if !state.reached_vs[w]
                        enqueue!(Q, w)
                        state.reached_vs[w] = true
                    end
                end
            end
        end
    end

    # Restore state
    for v in reached_vs
        state.reached_vs[v] = false
        for out_e in hg.hg_tail.v2he[v]
            state.hes_tail_count[out_e] += 1
        end
    end

    # Return reached vertices and hyperedges
    return reached_vs, reached_es
end

"""

"""
function backward_traceable(hg::H, target::Int, state::DiHyperPathState{T}) where {H <: AbstractDirectedHypergraph, T <: Real}
    # Priority queue of reached vertices
    Q = Queue{Int}()
    enqueue!(Q, target)
    
    state.reached_vs[target] = true

    # Which vertices/hyperedges have been reached?
    reached_vs = Set{Int}()
    reached_es = Set{Int}()

    while length(Q) > 0
        v = popfirst!(Q)
        push!(reached_vs, v)

        for in_e in keys(hg.hg_head.v2he[v])
            if !state.marked_hes[in_e]
                push!(reached_es, in_e)
                state.marked_hes[in_e] = true

                for w in keys(hg.hg_tail.he2v[in_e])
                    if !state.reached_vs[w]
                        enqueue!(Q, w)
                        state.reached_vs[w] = true
                    end
                end
            end
        end
    end

    # Restore state
    for v in reached_vs
        state.reached_vs[v] = false
    end
    for he in reached_es
        state.marked_hes[he] = false
    end

    # Return reached vertices and hyperedges
    return reached_vs, reached_es
end


"""

"""
function shortest_hyperpath_kk_heuristic(
    hg::H,
    source::Int,
    target::Int,
    cost_function::Function,
    hyperedge_weights::Vector{T}
) where {H <: AbstractDirectedHypergraph, T <: Real}
    state = initialize_dihyperpath_state(hg, hyperedge_weights)

    # Doubly reachable hyperedges
    dr_hes = intersect(
        forward_reachable(hg, source, state)[2],
        backward_traceable(hg, target, state)[2]
    )

    # Eliminate non-doubly reachable hyperedges
    hg_copy = deepcopy(hg)
    hg_copy[:, Not(dr_hes)] .= nothing

    # Min-heap for hyperedges
    Hmin = MutableBinaryMinHeap{Tuple{Int, T}}()
    for out_e in keys(hg.hg_tail.v2he[source])
        # If only the source is needed for this hyperedge
        if length(hg.hg_tail.he2v[out_e]) == 1
            # TODO: what do cost functions need?
            state.edge_heap_points[out_e] = push!(Hmin, cost_function(hg, state, out_e))
        end
    end

    state.reached_vs[source] = true

    while length(Hmin) > 0
        e = pop!(Hmin)
        state.removed_hes[e] = true

        path = short_hyperpath_vhe(hg, source, e)
        # TODO: nature of cost function for hyperedge vs. cost function for path
        state.edge_costs[e] = cost_function(hg, state, path)

        out_edges = Set{Int}()
        for v in keys(hg.hg_head.he2v[e])
            for f in keys(hg.hg_tail.v2he[v])
                if !state.reached_vs[v]
                    state.hes_tail_count[f] -= 1
                end

                if !state.marked_hes[f]
                    push!(out_edges, f)
                    state.marked_hes[f] = true
                end
            end
            state.reached_vs[v] = true
        end

        for f in out_edges
            push!(state.he_inedges[f], e)
            if !isnothing(state.edge_heap_points[f]) && !state.removed_hes[f]
                update!(Hmin, state.edge_heap_points[f], cost_function(hg, state, short_hyperpath_vhe(hg, source, f)))
            elseif isnothing(state.edge_heap_points[f]) && state.hes_tail_count[f] == 0
                state.edge_heap_points[f] = push!(Hmin, cost_function(hg, state, short_hyperpath_vhe(hg, source, f)))
            end
        end
    end

    path = Set{Int}()
    cost = typemax(T)
    for in_e in keys(hg.hg_head.v2he[target])
        if !isnothing(state.edge_heap_points[in_e])
            p = short_hyperpath_vhe(hg, source, in_e)
            cost_p = cost_function(hg, state, p)
            if cost_p < cost
                path = p
                cost = cost_p
            end
        end
    end

    return path
end

function shortest_hyperpath_kk_heuristic(
    hg::H,
    source::Int,
    targets::Set{Int},
    cost_function::Function
) where {H <: AbstractDirectedHypergraph}

end

function shortest_hyperpath_kk_heuristic(
    hg::H,
    sources::Set{Int},
    target::Int,
    cost_function::Function
) where {H <: AbstractDirectedHypergraph}

end

function shortest_hyperpath_kk_heuristic(
    hg::H,
    sources::Set{Int},
    targets::Set{Int},
    cost_function::Function
) where {H <: AbstractDirectedHypergraph}

end

"""

"""
function short_hyperpath_vhe(hg::H, v::Int, he::Int) where {H <: AbstractDirectedHypergraph}
    Q = Queue{Int}()
    for e in state.he_inedges[he]
        enqueue!(Q, e)
        state.marked_hes[e] = true
    end

    superpath = Set{Int}(he)
    path = Set{Int}(he)

    # Construct (likely redundant) superpath by backtracking from target
    while length(Q) > 0
        e = popfirst!(Q)
        push!(superpath, e)

        for f in state.he_inedges[e]
            if !state.marked_hes[f]
                enqueue!(Q, f)
                state.marked_hes[f] = true
            end
        end
    end

    for e in superpath
        state.marked_hes[e] = false
    end

    # TODO: try to be more clever about this
    superpath = sort(collect(superpath), by=x -> state.edge_costs[x], rev=true)
    hg_copy = deepcopy(hg)
    # Eliminate all hyperedges not on superpath
    hg_copy[:, Not(superpath)] .= nothing

    # Remove target from superpath; does not make sense to remove target in the next stage
    filter!(x -> x != he, superpath)

    # Try to minimize the size of the path by eliminating unnecessary hyperedges
    for e in superpath
        hg_copy[:, e] .= nothing

        # Only if hyperedge is essential for reaching target,
        if !is_reachable(hg_copy, s, he)
            # Restore hyperedge to hypergraph copy
            hg_copy[:, e] .= hg[:, e]
            push!(path, e)
        end
    end

    return path
end

"""

"""
function is_reachable(hg::H, v::Int, he::Int) where {H <: AbstractDirectedHypergraph}
    # TODO
end

"""

"""
function all_hyperpaths(hg::H, source::Int, target::Int) where {H <: AbstractDirectedHypergraph}

end

function all_hyperpaths(hg::H, source::Int, targets::Set{Int}) where {H <: AbstractDirectedHypergraph}

end

function all_hyperpaths(hg::H, sources::Set{Int}, target::Int) where {H <: AbstractDirectedHypergraph}

end

function all_hyperpaths(hg::H, sources::Set{Int}, targets::Set{Int}) where {H <: AbstractDirectedHypergraph}

end