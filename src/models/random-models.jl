# TODO:
# - Erdős–Rényi–Gilbert model for (directed) hypergraphs?
# - Random preferential model for directed hypergraphs
#     - I think some new math might need to be done here


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
