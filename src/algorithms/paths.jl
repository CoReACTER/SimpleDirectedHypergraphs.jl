"""
    forward_reachable(
        hg::H,
        source::Int,
    ) where {H <: AbstractDirectedHypergraph}

    Traverses a hypergraph `hg` starting from vertex with index `source` to determine all other vertices and hyperedges
    that are reachable, following hyperedges along their forward direction (i.e., from tail to head).

"""
function forward_reachable(
    hg::H,
    source::Int,
) where {H<:AbstractDirectedHypergraph}
    # Priority queue of reached vertices
    Q = Queue{Int}()
    enqueue!(Q, source)

    reached_vs = BitVector(falses(nhv(hg)))
    reached_vs[source] = true

    hes_tail_count = length.(keys.(hg.hg_tail.he2v))

    # Which vertices/hyperedges have been reached?
    vs = Set{Int}()
    es = Set{Int}()

    while length(Q) > 0
        v = dequeue!(Q)
        push!(vs, v)

        for out_e in keys(hg.hg_tail.v2he[v])
            # Following pseudocode exactly. This feels awkward; how slow would it be to just query reached_vs?
            hes_tail_count[out_e] -= 1

            if hes_tail_count[out_e] == 0
                push!(es, out_e)

                for w in keys(hg.hg_head.he2v[out_e])
                    if !reached_vs[w]
                        enqueue!(Q, w)
                        reached_vs[w] = true
                    end
                end
            end
        end
    end

    # Return reached vertices and hyperedges
    return (vs, es)
end

"""
    backward_traceable(
        hg::H,
        target::Int,
    ) where {H <: AbstractDirectedHypergraph}

    Traverses a hypergraph `hg` starting from vertex with index `target` to determine all other vertices and hyperedges
    that are reachable, following hyperedges along their reverse direction (i.e., from head to tail).
"""
function backward_traceable(
    hg::H,
    target::Int,
) where {H<:AbstractDirectedHypergraph}
    # Priority queue of reached vertices
    Q = Queue{Int}()
    enqueue!(Q, target)

    reached_vs = BitVector(falses(nhv(hg)))
    reached_vs[target] = true

    marked_hes = BitVector(falses(nhe(hg)))

    # Which vertices/hyperedges have been reached?
    vs = Set{Int}()
    es = Set{Int}()

    while length(Q) > 0
        v = dequeue!(Q)
        push!(vs, v)

        for in_e in keys(hg.hg_head.v2he[v])
            if !marked_hes[in_e]
                push!(es, in_e)
                marked_hes[in_e] = true

                for w in keys(hg.hg_tail.he2v[in_e])
                    if !reached_vs[w]
                        enqueue!(Q, w)
                        reached_vs[w] = true
                    end
                end
            end
        end
    end

    # Return reached vertices and hyperedges
    return (vs, es)
end


"""
    shortest_hyperpath_kk_heuristic(
        hg::H,
        source::Int,
        target::Int,
        hyperedge_weights::Vector{T}
    ) where {H <: AbstractDirectedHypergraph, T <: Real}

    shortest_hyperpath_kk_heuristic(
        hg::DirectedHypergraph{T, V, E, D},
        source::Int,
        targets::Set{Int},
        hyperedge_weights::Vector{T}
    ) where {T <: Real, V, E, D <: AbstractDict{Int,T}}

    shortest_hyperpath_kk_heuristic(
        hg::DirectedHypergraph{T, V, E, D},
        sources::Set{Int},
        target::Int,
        hyperedge_weights::Vector{T}
    ) where {T <: Real, V, E, D <: AbstractDict{Int,T}}

    shortest_hyperpath_kk_heuristic(
        hg::DirectedHypergraph{T, V, E, D},
        sources::Set{Int},
        targets::Set{Int},
        hyperedge_weights::Vector{T}
    ) where {T <: Real, V, E, D <: AbstractDict{Int,T}}

    Implements the heuristic directed hypergraph pathfinding algorithm of Krieger & Kececioglu (2022),
    DOI: 10.1186/s13015-022-00217-9. This algorithm is not guaranteed to find the optimal pathway from `source` to
    `target` based on some nonnegative `hyperedge_weights`), but in practice, it produces the optimal pathway
    approximately 99% of the time.
    
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
    hyperedge_weights::Vector{T}
) where {H<:AbstractDirectedHypergraph,T<:Real}

    reached_vs = BitVector(falses(nhv(hg)))
    reached_vs[source] = true

    marked_hes = BitVector(falses(nhe(hg)))
    removed_hes = BitVector(falses(nhe(hg)))

    hyperedge_inedges = [Set{Int}() for _ in 1:nhe(hg)]
    hes_tail_count = length.(keys.(hg.hg_tail.he2v))

    hyperedge_costs = fill(typemax(T), nhe(hg))

    hyperedge_heap_points = Vector{Union{Nothing,Int}}(nothing, nhe(hg))

    # Verify that the target can be reached
    fr = forward_reachable(hg, source)
    @assert target ∈ fr[1]

    # Doubly reachable hyperedges
    dr_hes = sort!(collect(intersect(
        fr[2],
        backward_traceable(hg, target)[2]
    )))

    # Eliminate non-doubly reachable hyperedges
    hg_copy = deepcopy(hg)
    hg_copy.hg_tail[:, InvertedIndices.Not(dr_hes)] .= nothing
    hg_copy.hg_head[:, InvertedIndices.Not(dr_hes)] .= nothing

    # Min-heap for hyperedges
    Hmin = MutableBinaryMinHeap{Tuple{T, Int}}()
    for out_e in keys(hg.hg_tail.v2he[source])
        # If only the source is needed for this hyperedge
        if length(hg.hg_tail.he2v[out_e]) == 1
            hyperedge_heap_points[out_e] = push!(Hmin, (hyperedge_weights[out_e], out_e))
        end
    end

    while length(Hmin) > 0
        e = pop!(Hmin)[2]
        removed_hes[e] = true

        path = short_hyperpath_vhe(hg, source, e, hyperedge_inedges, hyperedge_costs)
        hyperedge_costs[e] = sum(hyperedge_weights[x] for x in path)

        out_edges = Set{Int}()
        for v in keys(hg.hg_head.he2v[e])
            for f in keys(hg.hg_tail.v2he[v])
                if !marked_hes[f]
                    if !reached_vs[v]
                        hes_tail_count[f] -= 1
                    end
                    
                    if hes_tail_count[f] == 0
                        push!(out_edges, f)
                        marked_hes[f] = true
                    end
                end
            end
            reached_vs[v] = true
        end

        for f in out_edges
            marked_hes[f] = false
        end

        for f in out_edges
            push!(hyperedge_inedges[f], e)
            if !isnothing(hyperedge_heap_points[f]) && !removed_hes[f]
                update!(
                    Hmin,
                    hyperedge_heap_points[f],
                    (
                        sum(
                            hyperedge_weights[x]
                            for x in short_hyperpath_vhe(hg, source, f, hyperedge_inedges, hyperedge_costs)
                        ),
                        f
                    )
                )
            elseif isnothing(hyperedge_heap_points[f]) && hes_tail_count[f] == 0
                hyperedge_heap_points[f] = push!(
                    Hmin,
                    (
                        sum(
                            hyperedge_weights[x]
                            for x in short_hyperpath_vhe(hg, source, f, hyperedge_inedges, hyperedge_costs)
                        ),
                        f
                    )
                )
            end
        end
    end

    path = Set{Int}()
    cost = typemax(T)
    for in_e in keys(hg.hg_head.v2he[target])
        if !isnothing(hyperedge_heap_points[in_e])
            p = short_hyperpath_vhe(hg, source, in_e, hyperedge_inedges, hyperedge_costs)
            cost_p = sum(hyperedge_weights[e] for e in p)
            if cost_p < cost
                path = p
                cost = cost_p
            end
        end
    end

    return path
end

function shortest_hyperpath_kk_heuristic(
    hg::DirectedHypergraph{T,V,E,D},
    source::Int,
    targets::Set{Int},
    hyperedge_weights::Vector{S}
) where {S<:Real,T<:Real,V,E,D<:AbstractDict{Int,T}}
    hg_copy = deepcopy(hg)

    # Add a single "metatarget" vertex to reformulate as single-source, single-sink pathfinding problem
    # The hyperedge from the targets to the metatarget will have a cost of 0 associated with it
    metatarget = SimpleHypergraphs.add_vertex!(hg_copy)
    meta_he = SimpleHypergraphs.add_hyperedge!(
        hg_copy;
        vertices_tail=D(x => convert(T, 0) for x in targets),
        vertices_head=D(metatarget => convert(T, 0))
    )

    path = shortest_hyperpath_kk_heuristic(
        hg_copy,
        source,
        metatarget,
        vcat(hyperedge_weights, convert(S, 0))
    )

    # Remove the fictitious hyperedge from the targets to the metatarget
    setdiff(path, Set{Int}(meta_he))
end

function shortest_hyperpath_kk_heuristic(
    hg::DirectedHypergraph{T,V,E,D},
    sources::Set{Int},
    target::Int,
    hyperedge_weights::Vector{S}
) where {S<:Real,T<:Real,V,E,D<:AbstractDict{Int,T}}
    hg_copy = deepcopy(hg)

    # Add a single "metasource" vertex to reformulate as single-source, single-sink pathfinding problem
    # The hyperedge from the metasource to the sources will have a cost of 0 associated with it
    metasource = SimpleHypergraphs.add_vertex!(hg_copy)
    meta_he = SimpleHypergraphs.add_hyperedge!(
        hg_copy;
        vertices_tail=D(metasource => convert(T, 0)),
        vertices_head=D(x => convert(T, 0) for x in sources)
    )

    path = shortest_hyperpath_kk_heuristic(
        hg_copy,
        metasource,
        target,
        vcat(hyperedge_weights, convert(S, 0))
    )

    # Remove the fictitious hyperedge from the metasource to the sources
    setdiff(path, Set{Int}(meta_he))
end

function shortest_hyperpath_kk_heuristic(
    hg::DirectedHypergraph{T,V,E,D},
    sources::Set{Int},
    targets::Set{Int},
    hyperedge_weights::Vector{S}
) where {S<:Real,T<:Real,V,E,D<:AbstractDict{Int,T}}
    hg_copy = deepcopy(hg)

    # Add a single "metasource" vertex to reformulate as single-source, single-sink pathfinding problem
    # The hyperedge from the metasource to the sources will have a cost of 0 associated with it
    metasource = SimpleHypergraphs.add_vertex!(hg_copy)
    meta_he_source = SimpleHypergraphs.add_hyperedge!(
        hg_copy;
        vertices_tail=D(metasource => convert(T, 0)),
        vertices_head=D(x => convert(T, 0) for x in sources)
    )

    # Add a single "metatarget" vertex to reformulate as single-source, single-sink pathfinding problem
    # The hyperedge from the targets to the metatarget will have a cost of 0 associated with it
    metatarget = SimpleHypergraphs.add_vertex!(hg_copy)
    meta_he_target = SimpleHypergraphs.add_hyperedge!(
        hg_copy;
        vertices_tail=D(x => convert(T, 0) for x in targets),
        vertices_head=D(metatarget => convert(T, 0))
    )

    path = shortest_hyperpath_kk_heuristic(
        hg_copy,
        metasource,
        metatarget,
        vcat(hyperedge_weights, [convert(S, 0), convert(S, 0)])
    )

    # Remove fictitious hyperedges
    setdiff(path, Set{Int}([meta_he_source, meta_he_target]))
end

"""
    short_hyperpath_vhe(
        hg::H,
        v::Int,
        he::Int,
    ) where {H <: AbstractDirectedHypergraph, T <: Real}

    Obtain a (relatively, but not necessarily optimally) short hyperpath in hypergraph `hg` from a vertex `v` to a
    hyperedge `he`, `short_hyperpath_vhe` uses a greedy algorithm to first select hyperedges for a superpath and then
    prune unnecessary hyperedges to achieve a (generally shorter) hyperpath. 
"""
function short_hyperpath_vhe(
    hg::H,
    v::Int,
    he::Int,
    he_inedges::Vector{Set{Int}},
    he_costs::Vector{T}
) where {H<:AbstractDirectedHypergraph, T<:Real}
    marked_hes = BitVector(falses(nhe(hg)))

    Q = Queue{Int}()
    for e in he_inedges[he]
        enqueue!(Q, e)
        marked_hes[e] = true
    end

    superpath = Set{Int}(he)
    path = Set{Int}(he)

    # Construct (likely redundant) superpath by backtracking from target
    while length(Q) > 0
        e = dequeue!(Q)
        push!(superpath, e)

        for f in he_inedges[e]
            if !marked_hes[f]
                enqueue!(Q, f)
                marked_hes[f] = true
            end
        end
    end

    # TODO: try to be more clever about this
    superpath = sort(collect(superpath), by=x -> he_costs[x], rev=true)
    hg_copy = deepcopy(hg)
    # Eliminate all hyperedges not on superpath
    hg_copy.hg_tail[:, InvertedIndices.Not(superpath)] .= nothing
    hg_copy.hg_head[:, InvertedIndices.Not(superpath)] .= nothing

    # Remove target from superpath; does not make sense to remove target in the next stage
    filter!(x -> x != he, superpath)

    # Try to minimize the size of the path by eliminating unnecessary hyperedges
    for e in superpath
        hg_copy.hg_tail[:, e] .= nothing
        hg_copy.hg_head[:, e] .= nothing

        # Only if hyperedge is essential for reaching target,
        if !is_reachable(hg_copy, v, he, :hyperedge)
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
        source::Int,
        target::Int,
        target_type::Symbol
    ) where {H<:AbstractDirectedHypergraph}

    Use `forward_reachable` to determine if `target` (either a vertex index, if `target_type === :vertex` or a
    hyperedge index if `target_type === :hyperedge`) can be reached from vertex with index`source` in directed
    hypergraph `hg`.
"""
function is_reachable(
    hg::H,
    source::Int,
    target::Int,
    target_type::Symbol
) where {H<:AbstractDirectedHypergraph}
    @assert target_type ∈ [:vertex, :hyperedge] "`target_type` must be :vertex or :hyperedge"

    fr = forward_reachable(hg, source)

    #TODO: make a short-circuiting version of this. Can't for the life of me figure out why this isn't working...
    if target_type === :vertex
        return target ∈ fr[1]
    else
        return target ∈ fr[2]
    end
end

"""
    get_hyperpath(
        hg::H,
        source::Int,
        target::Int,
        out::Set{Int}
    ) where {H <: AbstractDirectedHypergraph}

    If one exists, obtain a hyperpath in directed hypergraph `hg` from a source vertex with index `source` to a target
    vertex with index `target`. The hyperpath cannot include any hyperedge with index included in the set `out`.
"""
function get_hyperpath(hg::H, source::Int, target::Int, out::Set{Int}) where {H<:AbstractDirectedHypergraph}
    # Remove excluded hyperedges
    hg_copy = deepcopy(hg)
    inds = sort(collect(out))
    hg_copy.hg_tail[:, inds] .= nothing
    hg_copy.hg_head[:, inds] .= nothing

    reached_vs, reached_es = forward_reachable(hg_copy, source)

    # Path does not exist
    if target ∉ reached_vs
        return Set{Int}()
    end

    path = Set{Int}()

    # Try to minimize the size of the path by eliminating unnecessary hyperedges
    for e in reached_es
        hg_copy.hg_tail[:, e] .= nothing
        hg_copy.hg_head[:, e] .= nothing

        # Only retain if hyperedge is essential for reaching target
        if !is_reachable(hg_copy, source, target, :vertex)
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
function all_hyperpaths(hg::H, source::Int, target::Int) where {H<:AbstractDirectedHypergraph}
    # Queue of subproblems
    # Subproblem is defined as a "out" set of hyperedges (which must not be present in a path) and "keep" hyperedges
    # which must be present in the path
    Q = Queue{Tuple{Set{Int},Set{Int}}}()

    paths = Set{Set{Int}}()

    # Start with no restrictions
    enqueue!(Q, (Set{Int}(), Set{Int}()))

    while length(Q) > 0
        out, keep = dequeue!(Q)

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

function all_hyperpaths(
    hg::DirectedHypergraph{T,V,E,D},
    source::Int,
    targets::Set{Int}
) where {T<:Real,V,E,D<:AbstractDict{Int,T}}
    hg_copy = deepcopy(hg)

    # Add a single "metatarget" vertex to reformulate as single-source, single-sink pathfinding problem
    # The hyperedge from the targets to the metatarget will have a cost of 0 associated with it
    metatarget = SimpleHypergraphs.add_vertex!(hg_copy)
    meta_he = SimpleHypergraphs.add_hyperedge!(
        hg_copy;
        vertices_tail=D(x => convert(T, 0) for x in targets),
        vertices_head=D(metatarget => convert(T, 0))
    )

    paths = all_hyperpaths(
        hg_copy,
        source,
        metatarget,
    )

    # Remove the fictitious hyperedge from the targets to the metatarget
    return Set(setdiff(p, Set{Int}(meta_he)) for p in paths)
end

function all_hyperpaths(
    hg::DirectedHypergraph{T,V,E,D},
    sources::Set{Int},
    target::Int
) where {T<:Real,V,E,D<:AbstractDict{Int,T}}
    hg_copy = deepcopy(hg)

    # Add a single "metasource" vertex to reformulate as single-source, single-sink pathfinding problem
    # The hyperedge from the metasource to the sources will have a cost of 0 associated with it
    metasource = SimpleHypergraphs.add_vertex!(hg_copy)
    meta_he = SimpleHypergraphs.add_hyperedge!(
        hg_copy;
        vertices_tail=D(metasource => convert(T, 0)),
        vertices_head=D(x => convert(T, 0) for x in sources)
    )

    paths = all_hyperpaths(
        hg_copy,
        metasource,
        target,
    )

    # Remove the fictitious hyperedge from the metasource to the sources
    return Set(setdiff(p, Set{Int}(meta_he)) for p in paths)
end

function all_hyperpaths(
    hg::DirectedHypergraph{T,V,E,D},
    sources::Set{Int},
    targets::Set{Int}
) where {T<:Real,V,E,D<:AbstractDict{Int,T}}
    hg_copy = deepcopy(hg)

    # Add a single "metasource" vertex to reformulate as single-source, single-sink pathfinding problem
    # The hyperedge from the metasource to the sources will have a cost of 0 associated with it
    metasource = SimpleHypergraphs.add_vertex!(hg_copy)
    meta_he_source = SimpleHypergraphs.add_hyperedge!(
        hg_copy;
        vertices_tail=D(metasource => convert(T, 0)),
        vertices_head=D(x => convert(T, 0) for x in sources)
    )

    # Add a single "metatarget" vertex to reformulate as single-source, single-sink pathfinding problem
    # The hyperedge from the targets to the metatarget will have a cost of 0 associated with it
    metatarget = SimpleHypergraphs.add_vertex!(hg_copy)
    meta_he_target = SimpleHypergraphs.add_hyperedge!(
        hg_copy;
        vertices_tail=D(x => convert(T, 0) for x in targets),
        vertices_head=D(metatarget => convert(T, 0))
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
    hyperedge_weights::Vector{T}
) where {H<:AbstractDirectedHypergraph, T<:Real}

    # First, verify that the problem is well-posed
    # That is, can `target` be reached from `source`
    @assert is_reachable(hg, source, target, :vertex)

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
            @constraint(model, sum([x[ih] for ih in in_hes]) >= x[i])
        end

        # Head-hitting inequalities
        if target ∉ keys(hg.hg_head.he2v[i])
            hits = Int[]
            for j in 1:nhe(hg)
                if i != j && length(intersect(Set(keys(hg.hg_head.he2v[i])), Set(keys(hg.hg_tail.he2v[j])))) >= 1
                    push!(hits, j)
                end
            end

            @constraint(model, sum([x[j] for j in hits]) >= x[i])
        end
    end

    # Target-production inequality
    @constraint(model, sum([x[e] for e in keys(hg.hg_head.v2he[target])]) >= 1)

    # Distance-based inequalities
    dist_ests = fill(typemax(T), nhv(hg))
    for v in forward_reachable(hg, source)[1]
        # TODO: this is inefficient
        # Currently, will repeat a lot of work
        # TODO: you are here; don't yet even understand what the problem is...
        dist_ests[v] = sum(
            [hyperedge_weights[e] for e in shortest_hyperpath_kk_heuristic(hg, source, v, hyperedge_weights)]
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

        @constraint(model, dot(x, cross) >= 1)
    end

    # Define objective function
    @objective(model, Min, dot(x, hyperedge_weights))

    return model, x, cuts, crosses
end

"""
    shortest_hyperpath_kk_ilp(
        hg::H,
        source::Int,
        target::Int,
        hyperedge_weights::Vector{T}
    ) where {H<:AbstractDirectedHypergraph, T<:Real}

    shortest_hyperpath_kk_ilp(
        hg::DirectedHypergraph{T,V,E,D},
        source::Int,
        targets::Set{Int},
        hyperedge_weights::Vector{S}
    ) where {S<:Real,T<:Real,V,E,D<:AbstractDict{Int,T}}

    shortest_hyperpath_kk_ilp(
        hg::DirectedHypergraph{T,V,E,D},
        sources::Set{Int},
        target::Int,
        hyperedge_weights::Vector{S}
    ) where {S<:Real,T<:Real,V,E,D<:AbstractDict{Int,T}}

    shortest_hyperpath_kk_ilp(
        hg::DirectedHypergraph{T,V,E,D},
        sources::Set{Int},
        targets::Set{Int},
        hyperedge_weights::Vector{S}
    ) where {S<:Real,T<:Real,V,E,D<:AbstractDict{Int,T}}

    Implements the exact directed hypergraph pathfinding algorithm of Krieger & Kececioglu (2023),
    DOI: 10.1089/cmb.2023.0242. This algorithm is guaranteed to find the optimal pathway from `source` to
    `target` based on some nonnegative `hyperedge_weights`, if and only if such a path exists.

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
    hyperedge_weights::Vector{T}
) where {H<:AbstractDirectedHypergraph,T<:Real}

    # TODO: do I need to carry `x` over like this? Not sure about variable scope
    model, x, cuts, crosses = initialize_ilp_model(hg, source, target, hyperedge_weights)

    optimize!(model)

    # Convert floating-point solution into BitVector
    # TODO: Is this necessary w/ JuMP? Or will the output really be binary? 
    solution = value.(x) .> 0.5

    # Check for s,t-cuts that are not crossed by the current solution
    new_cuts, new_crosses = expand_cuts(hg, source, target, cuts, crosses, solution)

    while length(new_cuts) > 0
        # Add new constraints to the model
        for cross in new_crosses
            @constraint(model, dot(x, cross) >= 1)
        end

        # Re-optimize model with new cut-constraints
        optimize!(model)
        solution = value.(x) .> 0.5

        new_cuts, new_crosses = expand_cuts(hg, source, target, new_cuts, new_crosses, solution)
    end

    return Set(findall(solution))
end

function shortest_hyperpath_kk_ilp(
    hg::DirectedHypergraph{T,V,E,D},
    source::Int,
    targets::Set{Int},
    hyperedge_weights::Vector{S}
) where {S<:Real,T<:Real,V,E,D<:AbstractDict{Int,T}}
    hg_copy = deepcopy(hg)

    # Add a single "metatarget" vertex to reformulate as single-source, single-sink pathfinding problem
    # The hyperedge from the targets to the metatarget will have a cost of 0 associated with it
    metatarget = SimpleHypergraphs.add_vertex!(hg_copy)
    meta_he = SimpleHypergraphs.add_hyperedge!(
        hg_copy;
        vertices_tail=D(x => convert(T, 0) for x in targets),
        vertices_head=D(metatarget => convert(T, 0))
    )

    path = shortest_hyperpath_kk_ilp(
        hg_copy,
        source,
        metatarget,
        vcat(hyperedge_weights, convert(S, 0))
    )

    # Remove the fictitious hyperedge from the targets to the metatarget
    return setdiff(path, Set{Int}(meta_he))
end

function shortest_hyperpath_kk_ilp(
    hg::DirectedHypergraph{T,V,E,D},
    sources::Set{Int},
    target::Int,
    hyperedge_weights::Vector{S}
) where {S<:Real,T<:Real,V,E,D<:AbstractDict{Int,T}}
    hg_copy = deepcopy(hg)

    # Add a single "metasource" vertex to reformulate as single-source, single-sink pathfinding problem
    # The hyperedge from the metasource to the sources will have a cost of 0 associated with it
    metasource = SimpleHypergraphs.add_vertex!(hg_copy)
    meta_he = SimpleHypergraphs.add_hyperedge!(
        hg_copy;
        vertices_tail=D(metasource => convert(T, 0)),
        vertices_head=D(x => convert(T, 0) for x in sources)
    )

    path = shortest_hyperpath_kk_ilp(
        hg_copy,
        metasource,
        target,
        vcat(hyperedge_weights, convert(S, 0))
    )

    # Remove the fictitious hyperedge from the metasource to the sources
    setdiff(path, Set{Int}(meta_he))
end

function shortest_hyperpath_kk_ilp(
    hg::DirectedHypergraph{T,V,E,D},
    sources::Set{Int},
    targets::Set{Int},
    hyperedge_weights::Vector{S}
) where {S<:Real,T<:Real,V,E,D<:AbstractDict{Int,T}}
    hg_copy = deepcopy(hg)

    # Add a single "metasource" vertex to reformulate as single-source, single-sink pathfinding problem
    # The hyperedge from the metasource to the sources will have a cost of 0 associated with it
    metasource = SimpleHypergraphs.add_vertex!(hg_copy)
    meta_he_source = SimpleHypergraphs.add_hyperedge!(
        hg_copy;
        vertices_tail=D(metasource => convert(T, 0)),
        vertices_head=D(x => convert(T, 0) for x in sources)
    )

    # Add a single "metatarget" vertex to reformulate as single-source, single-sink pathfinding problem
    # The hyperedge from the targets to the metatarget will have a cost of 0 associated with it
    metatarget = SimpleHypergraphs.add_vertex!(hg_copy)
    meta_he_target = SimpleHypergraphs.add_hyperedge!(
        hg_copy;
        vertices_tail=D(x => convert(T, 0) for x in targets),
        vertices_head=D(metatarget => convert(T, 0))
    )

    path = shortest_hyperpath_kk_ilp(
        hg_copy,
        metasource,
        metatarget,
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
    ) where {H<:AbstractDirectedHypergraph}

    A helper function for `shortest_hyperpath_kk_ilp`.

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
) where {H<:AbstractDirectedHypergraph}
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
                    new_cut_ev = Dict{Int,Set{Int}}()
                    new_cross_ev = Dict{Int,Vector{Bool}}()
                    greedy_v = 0
                    min_length = typemax(Int)
                    for v in keys(hg.hg_tail.he2v[e])
                        new_cut_ev[v] = setdiff(new_cut, Set(v))
                        new_cross_ev[v] = [
                            issubset(Set(keys(hg.hg_tail.he2v[i])), new_cut_ev[v]) &&
                            !issubset(Set(keys(hg.hg_head.he2v[i])), new_cut_ev[v])
                            for i in 1:nhe(hg)
                        ]
                        v_length = length(findall(map(!, old_cross) .&& new_cross_ev[v]))
                        if v_length < min_length
                            min_length = v_length
                            greedy_v = v
                        end
                    end

                    new_cut = new_cut_ev[greedy_v]
                    new_cross = new_cross_ev[greedy_v]
                end
            end
        end

        # If this is still a valid s,t-cut
        if !any(new_cross) && source ∈ new_cut
            push!(new_cuts, new_cut)
            push!(new_crosses, new_cross)
        end

    end

    return new_cuts, new_crosses
end