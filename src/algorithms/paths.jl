"""
    struct DiHyperPathState{T}(
        reached_vs::BitVector,
        marked_hes::BitVector,
        removed_hes::BitVector,
        hes_tail_count::Vector{Int},
        he_inedges::Vector{Set{Int}},
        edge_weights::Vector{T},
        edge_costs::Vector{T},
        edge_heap_points::Vector{Union{Nothing, Int}}
    ) where {T<:Real}

    Stores the state of a pathfinding operation.

    `reached_vs` and `marked_hes` track which vertices and hyperedges, respectively, have been visited and/or need to
    be searched. `removed_hes` tracks which hyperedges have been removed from the edge heap during the main heuristic
    pathfinding operation (see `shortest_hyperpath_kk_heuristic`). `hes_tail_count` is a vector containing the size of
    the tail of each hyperedge. `he_inedges` tracks which hyperedges flow into which other hyperedges (i.e., for a
    hyperedge `f`, the *inedges* are those hyperedges `e` where there exists some vertex `v` that is in both the head
    of `e` and the tail of `f`). `edge_weights` are a set of initially-defined costs associated with each directed
    hyperedge, and `edge_costs` are the (estimated) costs from the source vertex to each hyperedge. Finally,
    `edge_heap_points` tracks the references for each hyperedge in the heap used during heuristic pathfinding.

"""
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

"""
    initialize_dihyperpath_state(
        hg::H,
        hyperedge_weights::Vector{T}
    ) where {H <: AbstractDirectedHypergraph, T <: Real}

    Construct an initial state for directed hypergraph pathfinding on hypergraph `hg` with hyperedge weights (i.e.,
    costs) `hyperedge_weights`.
"""
function initialize_dihyperpath_state(
        hg::H,
        hyperedge_weights::Vector{T}
    ) where {H <: AbstractDirectedHypergraph, T <: Real}
    # Hyperedge weights need to be nonnegative
    @assert all(hyperedge_weights .>= 0)

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
    forward_reachable(
        hg::H,
        source::Int,
        state::DiHyperPathState{T}
    ) where {H <: AbstractDirectedHypergraph, T <: Real}

    Traverses a hypergraph `hg` starting from vertex with index `source` to determine all other vertices and hyperedges
    that are reachable, following hyperedges along their forward direction (i.e., from tail to head).

    `DiHyperPathState` `state` is used to track what vertices/hyperedges have been traversed.
"""
function forward_reachable(
    hg::H,
    source::Int,
    state::DiHyperPathState{T}
) where {H <: AbstractDirectedHypergraph, T <: Real}
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
    backward_traceable(
        hg::H,
        target::Int,
        state::DiHyperPathState{T}
    ) where {H <: AbstractDirectedHypergraph, T <: Real}

    Traverses a hypergraph `hg` starting from vertex with index `target` to determine all other vertices and hyperedges
    that are reachable, following hyperedges along their reverse direction (i.e., from head to tail).

    `DiHyperPathState` `state` is used to track what vertices/hyperedges have been traversed.
"""
function backward_traceable(
    hg::H,
    target::Int,
    state::DiHyperPathState{T}
) where {H <: AbstractDirectedHypergraph, T <: Real}
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
    shortest_hyperpath_kk_heuristic(
        hg::H,
        source::Int,
        target::Int,
        cost_function::Function,
        hyperedge_weights::Vector{T}
    ) where {H <: AbstractDirectedHypergraph, T <: Real}

    shortest_hyperpath_kk_heuristic(
        hg::DirectedHypergraph{T, V, E, D},
        source::Int,
        targets::Set{Int},
        cost_function::Function,
        hyperedge_weights::Vector{T}
    ) where {T <: Real, V, E, D <: AbstractDict{Int,T}}

    shortest_hyperpath_kk_heuristic(
        hg::DirectedHypergraph{T, V, E, D},
        sources::Set{Int},
        target::Int,
        cost_function::Function,
        hyperedge_weights::Vector{T}
    ) where {T <: Real, V, E, D <: AbstractDict{Int,T}}

    shortest_hyperpath_kk_heuristic(
        hg::DirectedHypergraph{T, V, E, D},
        sources::Set{Int},
        targets::Set{Int},
        cost_function::Function,
        hyperedge_weights::Vector{T}
    ) where {T <: Real, V, E, D <: AbstractDict{Int,T}}

    Implements the heuristic directed hypergraph pathfinding algorithm of Krieger & Kececioglu (2022),
    DOI: 10.1186/s13015-022-00217-9. This algorithm is not guaranteed to find the optimal pathway from `source` to
    `target` based on some `cost_function` (which must take in the hypergraph `hg`, the current `state`, and some
    pathway (collection of hyperedge)), but in practice, it produces the optimal pathway approximately 99% of the time.
    
    Note that, ostensibly, this algorithm only works for single-source, single-sink pathfinding (i.e., with a single
    `source` and a single `target`). If the user provides multiple `sources` and/or multiple `targets`, the
    multi-source/multi-sink problem will be reformulated as a single-source, single-sink problem by adding a
    *metasource* vertex (connected to all source vertices by a single, 0-cost hyperedge) and/or *metatarget* vertex
    (connected to all target vertices by a single, 0-cost hyperedge).
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

        path = short_hyperpath_vhe(hg, source, e, state)
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
                update!(Hmin, state.edge_heap_points[f], cost_function(hg, state, short_hyperpath_vhe(hg, source, f, state)))
            elseif isnothing(state.edge_heap_points[f]) && state.hes_tail_count[f] == 0
                state.edge_heap_points[f] = push!(Hmin, cost_function(hg, state, short_hyperpath_vhe(hg, source, f, state)))
            end
        end
    end

    path = Set{Int}()
    cost = typemax(T)
    for in_e in keys(hg.hg_head.v2he[target])
        if !isnothing(state.edge_heap_points[in_e])
            p = short_hyperpath_vhe(hg, source, in_e, state)
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
    hg::DirectedHypergraph{T, V, E, D},
    source::Int,
    targets::Set{Int},
    cost_function::Function,
    hyperedge_weights::Vector{T}
) where {T <: Real, V, E, D <: AbstractDict{Int,T}}
    hg_copy = deepcopy(hg)

    # Add a single "metatarget" vertex to reformulate as single-source, single-sink pathfinding problem
    # The hyperedge from the targets to the metatarget will have a cost of 0 associated with it
    metatarget = add_vertex!(hg_copy)
    meta_he = add_hyperedge!(
        hg_copy;
        vertices_tail=D( x => convert(T, 0) for x in targets),
        vertices_head=D(metatarget, convert(T, 0))
    )

    path = shortest_hyperpath_kk_heuristic(
        hg_copy,
        source,
        metatarget,
        cost_function,
        vcat(hyperedge_weights, convert(T, 0))
    )

    # Remove the fictitious hyperedge from the targets to the metatarget
    setdiff(path, Set{Int}(meta_he))
end

function shortest_hyperpath_kk_heuristic(
    hg::DirectedHypergraph{T, V, E, D},
    sources::Set{Int},
    target::Int,
    cost_function::Function,
    hyperedge_weights::Vector{T}
) where {T <: Real, V, E, D <: AbstractDict{Int,T}}
    hg_copy = deepcopy(hg)

    # Add a single "metasource" vertex to reformulate as single-source, single-sink pathfinding problem
    # The hyperedge from the metasource to the sources will have a cost of 0 associated with it
    metasource = add_vertex!(hg_copy)
    meta_he = add_hyperedge!(
        hg_copy;
        vertices_tail=D(metasource, convert(T, 0)),
        vertices_head=D( x => convert(T, 0) for x in sources)
    )

    path = shortest_hyperpath_kk_heuristic(
        hg_copy,
        metasource,
        target,
        cost_function,
        vcat(hyperedge_weights, convert(T, 0))
    )

    # Remove the fictitious hyperedge from the metasource to the sources
    setdiff(path, Set{Int}(meta_he))
end

function shortest_hyperpath_kk_heuristic(
    hg::DirectedHypergraph{T, V, E, D},
    sources::Set{Int},
    targets::Set{Int},
    cost_function::Function,
    hyperedge_weights::Vector{T}
) where {T <: Real, V, E, D <: AbstractDict{Int,T}}
    hg_copy = deepcopy(hg)

    # Add a single "metasource" vertex to reformulate as single-source, single-sink pathfinding problem
    # The hyperedge from the metasource to the sources will have a cost of 0 associated with it
    metasource = add_vertex!(hg_copy)
    meta_he_source = add_hyperedge!(
        hg_copy;
        vertices_tail=D(metasource, convert(T, 0)),
        vertices_head=D( x => convert(T, 0) for x in sources)
    )

    # Add a single "metatarget" vertex to reformulate as single-source, single-sink pathfinding problem
    # The hyperedge from the targets to the metatarget will have a cost of 0 associated with it
    metatarget = add_vertex!(hg_copy)
    meta_he_target = add_hyperedge!(
        hg_copy;
        vertices_tail=D( x => convert(T, 0) for x in targets),
        vertices_head=D(metatarget, convert(T, 0))
    )

    path = shortest_hyperpath_kk_heuristic(
        hg_copy,
        metasource,
        metatarget,
        cost_function,
        vcat(hyperedge_weights, [convert(T, 0), convert(T, 0)])
    )

    # Remove fictitious hyperedges
    setdiff(path, Set{Int}([meta_he_source, meta_he_target]))
end

"""
    short_hyperpath_vhe(
        hg::H,
        v::Int,
        he::Int,
        state::DiHyperPathState{T}
    ) where {H <: AbstractDirectedHypergraph, T <: Real}

    Obtain a (relatively, but not necessarily optimally) short hyperpath in hypergraph `hg` from a vertex `v` to a
    hyperedge `he`, using `state` to track the hypergraph traversal during pathfinding. `short_hyperpath_vhe` uses a
    greedy algorithm to first select hyperedges for a superpath and then prune unnecessary hyperedges to achieve a
    (generally shorter) hyperpath. 
"""
function short_hyperpath_vhe(
    hg::H,
    v::Int,
    he::Int,
    state::DiHyperPathState{T}
) where {H <: AbstractDirectedHypergraph, T <: Real}
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
        if !is_reachable(hg_copy, v, he, state)
            # Restore hyperedge to hypergraph copy
            hg_copy[:, e] .= hg[:, e]
            push!(path, e)
        end
    end

    return path
end

"""
    is_reachable(
        hg::H,
        source_v::Int,
        target_he::Int,
        state::DiHyperPathState{T}
    ) where {H <: AbstractDirectedHypergraph, T <: Real}

    A short-circuiting version of `forward_reachable` that returns `true` if a hyperpath in hypergraph `hg` exists from
    a source vertex with index `source_v` to a target hyperedge with index `target_he`. A `DiHyperPathState` object
    `state` is used to keep track of what vertices and hyperedges have been visited during a traversal.
"""
function is_reachable(
    hg::H,
    source_v::Int,
    target_he::Int,
    state::DiHyperPathState{T}
) where {H <: AbstractDirectedHypergraph, T <: Real}
    # Priority queue of reached vertices
    Q = Queue{Int}()
    enqueue!(Q, source_v)
    
    state.reached_vs[source_v] = true

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
                if out_e == target_he
                    return true
                end

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

    return false
end

"""
    is_reachable(
        hg::H,
        source_v::Int,
        target_v::Int,
        state::DiHyperPathState{T}
    ) where {H <: AbstractDirectedHypergraph, T <: Real}

    A short-circuiting version of `forward_reachable` that returns `true` if a hyperpath in hypergraph `hg` exists from
    a source vertex with index `source_v` to a target vertex with index `target_v`. A `DiHyperPathState` object
    `state` is used to keep track of what vertices and hyperedges have been visited during a traversal.
"""
function is_reachable(
    hg::H,
    source_v::Int,
    target_v::Int,
    state::DiHyperPathState{T}
) where {H <: AbstractDirectedHypergraph, T <: Real}
    # Priority queue of reached vertices
    Q = Queue{Int}()
    enqueue!(Q, source_v)
    
    state.reached_vs[source_v] = true

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
                    if w == target_v
                        return true
                    end

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

    return false
end

"""
    get_hyperpath(hg::H, source::Int, target::Int, out::Set{Int}) where {H <: AbstractDirectedHypergraph}

    If one exists, obtain a hyperpath in directed hypergraph `hg` from a source vertex with index `source` to a target
    vertex with index `target`. The hyperpath cannot include any hyperedge with index included in the set `out`.
"""
function get_hyperpath(hg::H, source::Int, target::Int, out::Set{Int}) where {H <: AbstractDirectedHypergraph}
    # Remove excluded hyperedges
    hg_copy = deepcopy(hg)
    hg_copy[:, sort(collect(out))] .= nothing

    weights = ones(nhe(hg_copy))
    state = initialize_dihyperpath_state(hg_copy, weights)

    reached_vs, reached_es = forward_reachable(hg_copy, source, state)

    # Path does not exist
    if target ∉ reached_vs
        return Set{Int}()
    end
    
    path = Set{Int}()

    # Try to minimize the size of the path by eliminating unnecessary hyperedges
    for e in reached_es
        hg_copy[:, e] .= nothing

        # Only if hyperedge is essential for reaching target,
        if !isreachable(hg_copy, source, target, state)
            # Restore hyperedge to hypergraph copy
            hg_copy[:, e] .= hg[:, e]
            push!(path, e)
        end
    end

    return path
end

"""
    all_hyperpaths(hg::H, source::Int, target::Int) where {H <: AbstractDirectedHypergraph}

    all_hyperpaths(hg::H, source::Int, targets::Set{Int}) where {H <: AbstractDirectedHypergraph}

    all_hyperpaths(hg::H, sources::Set{Int}, target::Int) where {H <: AbstractDirectedHypergraph}    

    all_hyperpaths(hg::H, sources::Set{Int}, targets::Set{Int}) where {H <: AbstractDirectedHypergraph}

    Exhaustively (but efficiently) generate all hyperpaths in directed hypergraph `hg` from some source(s) to some
    target(s), using the algorithm described by Krieger & Kececioglu (2022), DOI: 10.1186/s13015-022-00217-9. 
    
    Note that, ostensibly, this algorithm only works for single-source, single-sink pathfinding (i.e., with a single
    `source` index and a single `target` index). If the user provides multiple `sources` and/or multiple `targets`, the
    multi-source/multi-sink problem will be reformulated as a single-source, single-sink problem by adding a
    *metasource* vertex (connected to all source vertices by a single hyperedge) and/or *metatarget* vertex
    (connected to all target vertices by a single hyperedge).
"""
function all_hyperpaths(hg::H, source::Int, target::Int) where {H <: AbstractDirectedHypergraph}
    # Queue of subproblems
    # Subproblem is defined as a "out" set of hyperedges (which must not be present in a path) and "keep" hyperedges
    # which must be present in the path
    Q = Queue{Tuple{Set{Int}, Set{Int}}}()
    
    paths = Set{Set{Int}}()

    # Start with no restrictions
    enqueue!(Q, (Set{Int}(), Set{Int}()))

    while length(Q) > 0
        out, keep = popfirst!(Q)

        path = get_hyperpath(hg, source, target, out)

        if length(path) > 0 && path ∉ paths
            push!(paths, path)

            k = keep
            for e in setdiff(path, keep)
                enqueue!(Q, (union(out, e), k))
                push!(k, e)
            end
        end
    end

    paths
end

function all_hyperpaths(hg::H, source::Int, targets::Set{Int}) where {H <: AbstractDirectedHypergraph}
    hg_copy = deepcopy(hg)

    # Add a single "metatarget" vertex to reformulate as single-source, single-sink pathfinding problem
    # The hyperedge from the targets to the metatarget will have a cost of 0 associated with it
    metatarget = add_vertex!(hg_copy)
    meta_he = add_hyperedge!(
        hg_copy;
        vertices_tail=D( x => convert(T, 0) for x in targets),
        vertices_head=D(metatarget, convert(T, 0))
    )

    paths = all_hyperpaths(
        hg_copy,
        source,
        metatarget,
    )

    # Remove the fictitious hyperedge from the targets to the metatarget
    return Set(setdiff(p, Set{Int}(meta_he)) for p in paths)
end

function all_hyperpaths(hg::H, sources::Set{Int}, target::Int) where {H <: AbstractDirectedHypergraph}
    hg_copy = deepcopy(hg)

    # Add a single "metasource" vertex to reformulate as single-source, single-sink pathfinding problem
    # The hyperedge from the metasource to the sources will have a cost of 0 associated with it
    metasource = add_vertex!(hg_copy)
    meta_he = add_hyperedge!(
        hg_copy;
        vertices_tail=D(metasource, convert(T, 0)),
        vertices_head=D( x => convert(T, 0) for x in sources)
    )

    paths = all_hyperpaths(
        hg_copy,
        metasource,
        target,
    )

    # Remove the fictitious hyperedge from the metasource to the sources
    return Set(setdiff(p, Set{Int}(meta_he)) for p in paths)
end

function all_hyperpaths(hg::H, sources::Set{Int}, targets::Set{Int}) where {H <: AbstractDirectedHypergraph}
    hg_copy = deepcopy(hg)

    # Add a single "metasource" vertex to reformulate as single-source, single-sink pathfinding problem
    # The hyperedge from the metasource to the sources will have a cost of 0 associated with it
    metasource = add_vertex!(hg_copy)
    meta_he_source = add_hyperedge!(
        hg_copy;
        vertices_tail=D(metasource, convert(T, 0)),
        vertices_head=D( x => convert(T, 0) for x in sources)
    )

    # Add a single "metatarget" vertex to reformulate as single-source, single-sink pathfinding problem
    # The hyperedge from the targets to the metatarget will have a cost of 0 associated with it
    metatarget = add_vertex!(hg_copy)
    meta_he_target = add_hyperedge!(
        hg_copy;
        vertices_tail=D( x => convert(T, 0) for x in targets),
        vertices_head=D(metatarget, convert(T, 0))
    )

    paths = all_hyperpaths(
        hg_copy,
        metasource,
        metatarget,
    )

    # Remove fictitious hyperedges
    return Set(setdiff(p, Set{Int}([meta_he_source, meta_he_target])) for p in paths)
end

function initialize_ilp_model(
    hg::H,
    source::Int,
    target::Int,
    cost_function::Function,
    hyperedge_weights::Vector{T}) where {H<:AbstractDirectedHypergraph, T<:Real}

    # Initialize state
    state = initialize_dihyperpath_state(hg, hyperedge_weights)
    
    # First, verify that the problem is well-posed
    # That is, can `target` be reached from `source`
    @assert is_reachable(hg, source, target, state)

    # Initialize integer linear programming model
    model = Model(GLPK.Optimizer)
    
    # Create one binary variable for each hyperedge in `hg`
    @variable(model, x[1:nhe(hg)], Bin)
    set_start_value.(x, 1)

    # Verify that all variables are bound to be either 0 (not present in hyperpath) or 1 (present in hyperpath)
    @assert all(is_binary.(x))

    # Define initial constraints
    for i in 1:nhe(hg)
        # Tail-covering inequalities
        for v in keys(hg.hg_tail.he2v[i])
            if v == source
                continue
            end

            in_hes = collect(keys(hg.hg_head.v2he[v]))
            @constraint(model, con_scalar, sum([x[ih] for ih in in_hes]) >= x[i])
        end

        # Head-hitting inequalities
        if target ∉ hg.hg_head.he2v[i]
            hits = Int[]
            for j in 1:nhe(hg)
                if i != j && length(intersect(Set(keys(hg.hg_head.he2v[i])), Set(keys(hg.hg_tail.he2v[j])))) >= 1
                    push!(hits, j)
                end
            end

            @constraint(model, con_scalar, sum([x[j] for j in hits]) >= x[i])
        end
    end

    # Target-production inequality
    @constraint(model, con_scalar, sum([x[e] for e in keys(hg.hg_head.v2he[target])]) >= 1)

    # Distance-based inequalities
    dist_ests = fill(typemax(T), nhv(hg))
    for v in forward_reachable(hg, source, state)[1]
        # TODO: this is inefficient
        # Currently, will repeat a lot of work
        dist_ests[v] = cost_function(
            hg,
            state,
            shortest_hyperpath_kk_heuristic(hg, source, v, cost_function, hyperedge_weights)
        )
    end

    cuts = Set{Int}[]
    crosses = BitVector[]

    for d in unique(values(dist_ests))
        if d > dist_ests[target]
            continue
        end

        # Construct an s,t-cut, where the source is within the cut set of vertices and the target is not in the set
        cut_d = Set(v for (v, dist) in enumerate(dist_ests) if dist < d)
        push!(cut_d, source)
        setdiff!(cut_d, Set(target))

        push!(cuts, cut_d)

        cross = [
            issubset(Set(keys(hg.hg_tail.he2v[i])), cut_d) && !issubset(Set(keys(hg.hg_head.he2v[i])), cut_d)
            for i in 1:nhe(hg)
        ]
        push!(crosses, BitVector(cross))

        @constraint(model, con_scalar, dot(x, cross) >= 1)
    end

    # Define objective function
    @objective(model, Min, dot(x, state.edge_weights))

    return model, x, cuts, crosses
    
end

"""
    shortest_hyperpath_kk_ilp(
        hg::H,
        source::Int,
        target::Int,
        cost_function::Function,
        hyperedge_weights::Vector{T}
    ) where {H<:AbstractDirectedHypergraph, T<:Real}

    shortest_hyperpath_kk_ilp(
        hg::H,
        source::Int,
        targets::Set{Int},
        cost_function::Function,
        hyperedge_weights::Vector{T}
    ) where {H<:AbstractDirectedHypergraph, T<:Real}

    shortest_hyperpath_kk_ilp(
        hg::H,
        sources::Set{Int},
        target::Int,
        cost_function::Function,
        hyperedge_weights::Vector{T}
    ) where {H<:AbstractDirectedHypergraph, T<:Real}

    shortest_hyperpath_kk_ilp(
        hg::H,
        sources::Set{Int},
        targets::Set{Int},
        cost_function::Function,
        hyperedge_weights::Vector{T}
    ) where {H<:AbstractDirectedHypergraph, T<:Real}

    Implements the exact directed hypergraph pathfinding algorithm of Krieger & Kececioglu (2023),
    DOI: 10.1089/cmb.2023.0242. This algorithm is guaranteed to find the optimal pathway from `source` to
    `target` based on some `cost_function` (which must take in the hypergraph `hg`, the current `state`, and some
    pathway (collection of hyperedge)), if and only if such a path exists.

    Note that, ostensibly, this algorithm only works for single-source, single-sink pathfinding (i.e., with a single
    `source` and a single `target`). If the user provides multiple `sources` and/or multiple `targets`, the
    multi-source/multi-sink problem will be reformulated as a single-source, single-sink problem by adding a
    *metasource* vertex (connected to all source vertices by a single, 0-cost hyperedge) and/or *metatarget* vertex
    (connected to all target vertices by a single, 0-cost hyperedge).
"""
function shortest_hyperpath_kk_ilp(
    hg::H,
    source::Int,
    target::Int,
    cost_function::Function,
    hyperedge_weights::Vector{T}
) where {H<:AbstractDirectedHypergraph, T<:Real}

    # TODO: do I need to carry `x` over like this? Not sure about variable scope
    model, x, cuts, crosses = initialize_ilp_model(hg, source, target, cost_function, hyperedge_weights)

    optimize!(model)
    
    # Convert floating-point solution into BitVector
    # TODO: Is this necessary w/ JuMP? Or will the output really be binary? 
    solution = value.(x) .> 0.5

    # Check for s,t-cuts that are not crossed by the current solution
    new_cuts, new_crosses = expand_cuts(hg, source, target, cuts, crosses, solution)

    while length(new_cuts) > 0
        # Add new constraints to the model
        for cross in new_crosses
            @constraint(model, con_scalar, dot(x, cross) >= 1)
        end

        # Re-optimize model with new cut-constraints
        optimize!(model)
        solution = value.(x) .> 0.5

        new_cuts, new_crosses = expand_cuts(hg, source, target, new_cuts, new_crosses, solution)
    end

    return Set(findall(solution))
end

function shortest_hyperpath_kk_ilp(
    hg::H,
    source::Int,
    targets::Set{Int},
    cost_function::Function,
    hyperedge_weights::Vector{T}
) where {H<:AbstractDirectedHypergraph, T<:Real}
    hg_copy = deepcopy(hg)

    # Add a single "metatarget" vertex to reformulate as single-source, single-sink pathfinding problem
    # The hyperedge from the targets to the metatarget will have a cost of 0 associated with it
    metatarget = add_vertex!(hg_copy)
    meta_he = add_hyperedge!(
        hg_copy;
        vertices_tail=D( x => convert(T, 0) for x in targets),
        vertices_head=D(metatarget, convert(T, 0))
    )

    path = shortest_hyperpath_kk_ilp(
        hg_copy,
        source,
        metatarget,
        cost_function,
        hyperedge_weights
    )

    # Remove the fictitious hyperedge from the targets to the metatarget
    return setdiff(path, Set{Int}(meta_he))
end

function shortest_hyperpath_kk_ilp(
    hg::H,
    sources::Set{Int},
    target::Int,
    cost_function::Function,
    hyperedge_weights::Vector{T}
) where {H<:AbstractDirectedHypergraph, T<:Real}
    hg_copy = deepcopy(hg)

    # Add a single "metasource" vertex to reformulate as single-source, single-sink pathfinding problem
    # The hyperedge from the metasource to the sources will have a cost of 0 associated with it
    metasource = add_vertex!(hg_copy)
    meta_he = add_hyperedge!(
        hg_copy;
        vertices_tail=D(metasource, convert(T, 0)),
        vertices_head=D( x => convert(T, 0) for x in sources)
    )

    path = shortest_hyperpath_kk_ilp(
        hg_copy,
        metasource,
        target,
        cost_function,
        vcat(hyperedge_weights, convert(T, 0))
    )

    # Remove the fictitious hyperedge from the metasource to the sources
    setdiff(path, Set{Int}(meta_he))
end

function shortest_hyperpath_kk_ilp(
    hg::H,
    sources::Set{Int},
    targets::Set{Int},
    cost_function::Function,
    hyperedge_weights::Vector{T}
) where {H<:AbstractDirectedHypergraph, T<:Real}
    hg_copy = deepcopy(hg)

    # Add a single "metasource" vertex to reformulate as single-source, single-sink pathfinding problem
    # The hyperedge from the metasource to the sources will have a cost of 0 associated with it
    metasource = add_vertex!(hg_copy)
    meta_he_source = add_hyperedge!(
        hg_copy;
        vertices_tail=D(metasource, convert(T, 0)),
        vertices_head=D( x => convert(T, 0) for x in sources)
    )

    # Add a single "metatarget" vertex to reformulate as single-source, single-sink pathfinding problem
    # The hyperedge from the targets to the metatarget will have a cost of 0 associated with it
    metatarget = add_vertex!(hg_copy)
    meta_he_target = add_hyperedge!(
        hg_copy;
        vertices_tail=D( x => convert(T, 0) for x in targets),
        vertices_head=D(metatarget, convert(T, 0))
    )

    path = shortest_hyperpath_kk_ilp(
        hg_copy,
        metasource,
        metatarget,
        cost_function,
        vcat(hyperedge_weights, [convert(T, 0), convert(T, 0)])
    )

    # Remove fictitious hyperedges
    setdiff(path, Set{Int}([meta_he_source, meta_he_target]))
end

"""
    expand_cuts(
        hg::H,
        source::Int,
        target::Int,
        cuts::Vector{Set{Int}},
        crosses::Vector{BitVector},
        curr_sol::BitVector
    ) where {H <: AbstractDirectedHypergraph}

    A helper function for `shortest_hyperpath_kk_ilp`

    Efficiently identify new `source`-`target` cuts of a directed hypergraph `hg` that are violated by current solution
    of the integer linear programming problem `curr_sol`.
"""
function expand_cuts(
    hg::H,
    source::Int,
    target::Int,
    cuts::Vector{Set{Int}},
    crosses::Vector{BitVector},
    curr_sol::BitVector
) where {H <: AbstractDirectedHypergraph}
    new_cuts = Set{Int}[]
    new_crosses = BitVector[]

    for (old_cut, old_cross) in zip(cuts, crosses)
        # Source augmentation
        new_cut = old_cut
        new_cross = old_cross
        for e in 1:nhe(hg)
            # If `e` is an active hyperedge that crosses the current cut, add the head of `e` to the cut so `e` no
            # longer crosses it
            if curr_sol[e] && new_cross[e]
                new_cut = union(new_cut, Set(keys(hg.hg_head.he2v[e])))
                new_cross = [
                    issubset(Set(keys(hg.hg_tail.he2v[i])), new_cut) && !issubset(Set(keys(hg.hg_head.he2v[i])), new_cut)
                    for i in 1:nhe(hg)
                ]
            end
        end
        
        # If this is a valid s,t-cut that no active hyperedges cross, add it to the new cut list
        if !(any(new_cross) || target ∈ new_cut)
            push!(new_cuts, new_cut)
            push!(new_crosses, new_cross)
        end

        # Sink augmentation
        # TODO: you could (should) combine these two into one loop
        new_cut = old_cut
        new_cross = old_cross
        while any(curr_sol .&& new_cross)
            for e in 1:nhe(hg)
                # If `e` is an active hyperedge that crosses the current cut, remove vertices from the tail of `e` to the
                # cut so `e` no longer crosses it
                if curr_sol[e] && new_cross[e]
                    # Greedily pick the vertex in the tail of `e` that causes the fewest hyperedges to newly cross this cut
                    new_cut_ev = Dict{Int, Set{Int}}()
                    new_cross_ev = Dict{Int, Vector{Bool}}()
                    for v in keys(hg.hg_tail.he2v[e])
                        new_cut_ev[v] = setdiff(new_cut, Set(v))
                        new_cross_ev[v] = [
                            issubset(Set(keys(hg.hg_tail.he2v[i])), new_cut_ev[v]) &&
                            !issubset(Set(keys(hg.hg_head.he2v[i])), new_cut_ev[v])
                            for i in 1:nhe(hg)
                        ]
                    end

                    greedy_v = minimum(q -> length(findall(!old_cross .&& new_cross_ev[q])), keys(new_cross_ev))

                    new_cut = new_cut_ev[greedy_v]
                    new_cross = new_cross_ev[greedy_v]
                end
            end
        end

        # If this is still a valid s,t-cut
        if !any(new_corss) && source ∈ new_cut
            push!(new_cuts, new_cut)
            push!(new_crosses, new_cross)
        end

    end

    return new_cuts, new_crosses
end