"""
    random_model(nVertices::Int, nEdges::Int, HType::Type{H}) where {H<:AbstractDirectedHypergraph}

Generate a *random* directed hypergraph (in the style of Erdős–Rényi random graphs) without any structural constraints.

# The Algorithm

Given two integer parameters *nVertices* and *nEdges* (the number of nodes and hyperedges, respectively),
the algorithm computes - for each hyperedge *he={1,...,m}* -
two random numbers *s_t ϵ [1, n]* (i.e., the size of the hyperedge tail) and *s_h ϵ [1, n]*
(i.e., the size of the hyperedge head).
Then, the algorithm selects uniformly at random *s_t* vertices from *V* to be added in the tail of *he*
and *s_h* vertices from *V* to be added into the head of *he*.
If *no_self_loops* is true (default false), then vertices will be chosen such that the same vertex *v*
cannot appear in both the head and the tail of the same hyperedge *he*.
"""
function SimpleHypergraphs.random_model(
    nVertices::Int,
    nEdges::Int,
    HType::Type{H};
    no_self_loops::Bool=false
) where {H<:AbstractDirectedHypergraph}
    if no_self_loops && nVertices == 1
        # Impossible; all directed hyperedges with non-empty tails & heads will be self-loops
        error("Impossible to avoid self-loops in non-empty directed hyperedges in a directed hypergraph with one (1) vertex!")
    end

    mx_tail = Matrix{Union{Nothing,Bool}}(nothing, nVertices, nEdges)
    mx_head = Matrix{Union{Nothing,Bool}}(nothing, nVertices, nEdges)
    if nEdges == 0
        return HType(nVertices, nEdges)
    end

    for e in 1:size(mx_tail,2)
        if no_self_loops
            # To avoid self-loops, can't have all 
            nv = rand(1:size(mx_tail,1)-1)
        else
            nv = rand(1:size(mx_tail,1))
        end
        mx_tail[sample(1:size(mx_tail,1), nv;replace=false), e] .= true

        if no_self_loops
            valid_indices = collect(setdiff(Set(1:size(mx_tail, 1)), Set(findall(x->x===true, mx_tail[:, e]))))
            nv = rand(1:length(valid_indices))
            mx_head[sample(valid_indices, nv;replace=false), e] .= true
        else
            nv = rand(1:size(mx_head,1))
            mx_head[sample(1:size(mx_head,1), nv;replace=false), e] .= true
        end
    end

    h = HType(mx_tail, mx_head)
    if all(length.(h.hg_tail.v2he) .> 0) && all(length.(h.hg_head.v2he) .> 0)
        return h
    else
        return random_model(nVertices, nEdges, HType, no_self_loops=no_self_loops)
    end
end


# TODO: overwrite issue
"""
    random_kuniform_model(nVertices::Int, nEdges::Int, k::Int; HType::Type{H}=DirectedHypergraph, no_self_loops::Bool=false) where {H<:AbstractDirectedHypergraph}

Generates a *k*-uniform directed hypergraph, i.e., an hypergraph where each hyperedge has size *k = k_tail + k_head*.
In this implementation, *k_tail* and *k_head* are not necessarily equal.

# The Algorithm

The algorithm proceeds as the *random_model*, forcing the size of each hyperedge equal to *k*.
"""
function SimpleHypergraphs.random_kuniform_model(
    nVertices::Int,
    nEdges::Int,
    k::Int,
    HType::Type{H}=DirectedHypergraph;
    no_self_loops::Bool=false
) where {H<:AbstractDirectedHypergraph}
    mx_tail = Matrix{Union{Nothing,Bool}}(nothing, nVertices,nEdges)
    mx_head = Matrix{Union{Nothing,Bool}}(nothing, nVertices,nEdges)
    
    if no_self_loops && nVertices == 1
        # Impossible; all directed hyperedges with non-empty tails & heads will be self-loops
        error("Impossible to avoid self-loops in non-empty directed hyperedges in a directed hypergraph with one (1) vertex!")
    end

    for e in 1:size(mx_tail,2)
        if no_self_loops
            k_tail = rand(1:k-1)
        else
           k_tail = rand(1:k)
        end
        k_head = k - k_tail

        mx_tail[sample(1:size(mx_tail,1), k_tail;replace=false), e] .= true

        if no_self_loops
            valid_indices = collect(setdiff(Set(1:size(mx_tail, 1)), Set(findall(x->x===true, mx_tail[:, e]))))
            mx_head[sample(valid_indices, k_head;replace=false), e] .= true
        else
            mx_head[sample(1:size(mx_head,1), k_head;replace=false), e] .= true
        end
    end
    
    HType(mx_tail, mx_head)
end


"""
    random_dregular_model(
        nVertices::Int,
        nEdges::Int,
        d::Int,
        HType::Type{H};
        no_self_loops::Bool=false
    ) where {H<:AbstractDirectedHypergraph}
Generates a *d*-regular directed hypergraph, where each node has degree *d*.

# The Algorithm

The algorithm exploits the *k*-uniform approach described for the *random_kuniform_model* method
to build a *d*-regular hypergraph *H* having *nVertices* nodes and *nEdges* edges.
It returns the hypergraph *H^* dual of *H*.
"""
function SimpleHypergraphs.random_dregular_model(
    nVertices::Int,
    nEdges::Int,
    d::Int,
    HType::Type{H};
    no_self_loops::Bool = false
) where {H<:AbstractDirectedHypergraph}
    if no_self_loops && nVertices == 1
        # Impossible; all directed hyperedges with non-empty tails & heads will be self-loops
        error("Impossible to avoid self-loops in non-empty directed hyperedges in a directed hypergraph with one (1) vertex!")
    end

    mx_tail = Matrix{Union{Nothing,Bool}}(nothing, nVertices, nEdges)
    mx_head = Matrix{Union{Nothing,Bool}}(nothing, nVertices, nEdges)
    
    for v in 1:size(mx_tail,1)
        d_tail = rand(1:d)
        d_head = d - d_tail

        if d_tail > 0
            mx_tail[v, sample(1:size(mx_tail,2), d_tail;replace=false)] .= true
        end

        if d_head > 0
            if no_self_loops
                valid_indices = collect(setdiff(Set(1:size(mx_tail, 2)), Set(findall(x->x===true, mx_tail[v, :]))))
                mx_head[v, sample(valid_indices, d_head;replace=false)] .= true
            else
                mx_head[v, sample(1:size(mx_head,2), d_head;replace=false)] .= true
            end
        end
    end

    HType(mx_tail, mx_head)
end

function random_crn_model(
    nVertices::Int,
    nEdges::Int;
    nTrials::Int=10,
    maxTailSize::Int=nVertices,
    maxHeadSize::Int=nVertices,
    allow_catalysts::Bool=true,
    symmetric::Bool=true,
    closed_system::Bool=true,
    # ensure_connected::Symbol=:none,
    rng::AbstractRNG=Random.default_rng()
)
    # No point making a random CRN with 0 vertices or 0 hyperedges
    @assert nVertices > 0
    @assert nEdges > 0

    # Algorithm either does not check for connectivity (:none), checks for weak connectivity (:weak), or checks for
    # strong connectivity (:strong)
    # @assert ensure_connected ∈ [:none, :weak, :strong]
    
    # If hyperedge generation cannot be attempted, algorithm is futile
    @assert nTrials >= 1

    # Cannot have a symmetric CRN with an odd number of reactions
    @assert !(symmetric && isodd(nEdges))

    # In symmetric case, we must restrict the tail/head size
    # If maxTailSize != maxHeadSize initially, then there will be some hyperedges that are valid in the forwards
    # direction but invalid in the reverse direction
    if symmetric && (maxTailSize != maxHeadSize)
        @warn "`maxTailSize != maxHeadSize; changing both to `min(maxTailSize, maxHeadSize)` to ensure symmetry."
        maxSize = min(maxTailSize, maxHeadSize)
        maxTail = maxSize
        maxHead = maxSize
    else
        maxTail = maxTailSize
        maxHead = maxHeadSize
    end

    # Directed hypergraph with `Hij = (t, h)`, where `t` and `h` are the multiplicities of vertex/species `i` in the
    # tail and head of hyperedge/reaction `j`, respectively.
    H = DirectedHypergraph{Int}(nVertices, nEdges)

    # Uniformly select `nEdges` directed hyperedges `e` obeying:
    # - `minTailSize ≤ |tail(e)| ≤ maxTailSize`
    # - `minHeadSize ≤ |head(e)| ≤ maxHeadSize`
    # - (if `!allow_catalysts`) `tail(e) ∩ head(e) = ∅`
    # - Let `mnet = {m_head(e, i) - m_tail(e, i) | 1 ≤ i ≤ nVertices}`, where `m_head(e, i)` is the multiplicity (i.e.,
    #   the stoichiometric coefficient) of species/vertex `i` in the head of reaction/hyperedge `e` and `m_tail(e, i)`
    #   is, accordingly, the multiplicity of `i` in the tail of `e`. Then,
    #       - (if `!allow_import_export`), `∃ x1 ∈ mnet, x1 > 0 ⩓ ∃ x2 ∈ mnet, x2 < 0`
    #       - (in general) `∃ x ∈ mnet, x ≠ 0`
    # After `nTrials` attempts, if `nEdges` valid, unique hyperedges haven't been generated, give up and move on
    i = 0
    n = 0

    if symmetric
        target = nEdges / 2
    else
        target = nEdges
    end

    # If system is closed, no reaction (hyperedge) can have an empty tail or head
    # If the system is open, then reactions can have empty tail or empty heads
    # We still don't allow reactions with non-empty tails and heads where the net reaction has an empty tail or head
    # (see below)
    if closed_system
        minTailHead = 1
    else
        minTailHead = 0
    end

    chosen = Tuple{Accumulator{Int64, Int64}, Accumulator{Int64, Int64}}[]
    
    while n < target && i < nTrials
        for _ in 1:(target - n)
            # Select tail size
            nTail = rand(rng, minTailHead:maxTail)

            # Select tail vertices/species
            thisTail = counter(rand(rng, 1:nVertices, nTail))

            # Select head size
            nHead = rand(rng, minTailHead:maxHead)

            # Select & partition head vertices/species
            if allow_catalysts
                # Select any `nHead` vertices with replacement
                # During validation (below) we will ensure that this is not a futile reaction
                thisHead = counter(rand(rng, 1:nVertices, nHead))
            else
                # If no catalysts are allowed
                allowed_vertices = setdiff(1:nVertices, keys(thisTail))
                thisHead = counter(rand(rng, allowed_vertices, nHead))
            end

            valid = true
            
            # Check for duplication
            if (thisTail, thisHead) ∈ chosen
                valid = false
            end

            # Check for futility
            if allow_catalysts
                # If tail or head is empty, reaction is valid (must be open system)
                if !(length(thisTail) == 0 || length(thisHead) == 0)
                    fordiff = setdiff(thisHead, thisTail)
                    revdiff = setdiff(thisTail, thisHead)

                    # Otherwise, empty tail/head in *net* reaction is disallowed
                    if length(fordiff) == 0 || length(revdiff) == 0
                        valid = false
                    end
                end
            else
                # If catalysts not allowed, the only possible futile reaction is ∅ -> ∅
                if length(thisTail) == 0 && length(thisHead) == 0
                    valid = false
                end
            end

            if valid
                n += 1
                # Add hyperedge to hypergraph
                push!(chosen, (thisTail, thisHead))

                if symmetric
                    # Add reverse hyperedge to hypergraph
                    push!(chosen, (thisHead, thisTail))
                end
            end
        end

        i += 1
    end

    # Use the hyperedges generated to populate the hypergraph
    for (i, he) in enumerate(chosen)
        for (k, v) in he[1]
            H[1, k, i] = v
        end

        for (k, v) in he[2]
            H[2, k, i] = v
        end
    end

    # # Check connectivity
    # # If possible, try to ensure connectivity while maintaining network shape (i.e., w/ exactly `nEdges` hyperedges)
    # if ensure_connected === :weak

    # elseif ensure_connected === :strong

    # else
    #     # If user doesn't care about connectivity, return hypergraph directly
    #     return H
    # end

    return H

end

# TODO: preferential attachment
# TODO: hypercurveball algorithm
# TODO:
# - Erdős–Rényi–Gilbert model for (directed) hypergraphs?
# - Random preferential model for directed hypergraphs
#     - I think some new math might need to be done here
