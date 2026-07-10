# TODO: maybe more fancy file format and correctness checking should be done

struct EHGF_Format <: Abstract_HG_format end


"""
    hg_save(io::IO, h::H, format::EHGF_Format) where {H <: AbstractDirectedHypergraph}

Saves an undirected hypergraph `h` to an output stream `io` in `ehgf` format.

TODO: what to do about metadata?

"""
function SimpleHypergraphs.hg_save(io::IO, h::H, format::EHGF_Format) where {H <: AbstractDirectedHypergraph}
    
    h_size = Base.size(h)
    
    println(io, h_size[1], " ", h_size[2])
    for i in 1:h_size[2]
        tail_keys = sort(collect(keys(h.hg_tail.he2v[i])))
        head_keys = sort(collect(keys(h.hg_head.he2v[i])))
        print(
            io, 
            join(["$k=$(h.hg_tail.he2v[i][k])" for k in tail_keys], ' ')
        )
        print(io, " || ")
        print(
            io,
            join(["$k=$(h.hg_head.he2v[i][k])" for k in head_keys], ' ')
        )
        print(io, "\n")
    end
end


"""
    hg_save(io::IO, h::DirectedHypergraph, format::JSON_Format; pretty::Bool = false)

Saves a directed hypergraph `h` to an output stream `io` in `json` format.

If `h` has `Composite Types` either for vertex metadata or hyperedges metadata,
the user has to explicit tell the JSON package about it, for instance using:

TODO: this

See the (JSON.jl documentation)[https://juliaio.github.io/JSON.jl/stable/] for more details.

The `json` in output contains the following information (keys):

* `n` : number of vertices
* `k` : number of hyperedges
* `tail` : a matrix representation of the tails of `h`, where rows are vertices and columns are hyperedges
* `head` : a matrix representation of the heads of `h`, where rows are vertices and columns are hyperedges
* `v_meta` : vertices metadata
* `he_meta_tail` : metadata for hyperedge tails
* `he_meta_head` : metadata for hyperedge heads

"""
function SimpleHypergraphs.hg_save(io::IO, h::DirectedHypergraph, format::JSON_Format; pretty::Bool = false)
    json_hg = OrderedDict{Symbol, Any}()

    json_hg[:n] = nhv(h)
    json_hg[:k] = nhe(h)

    json_hg[:tail] = Matrix(h.hg_tail)
    json_hg[:head] = Matrix(h.hg_head)
    
    json_hg[:v_meta] = h.v_meta
    json_hg[:he_meta_tail] = h.he_meta_tail
    json_hg[:he_meta_head] = h.he_meta_head

    JSON.json(io, json_hg; pretty)
end

"""
    hg_save(
	io::IO,
	h::DirectedHypergraph{T, V, E, D},
	::HIF_Format;
	pretty::Bool=false
    ) where {T, V, E, D}

Saves a directed hypergraph `h` to an output stream `io` in `HIF` format (see Coll et al.,
DOI: 10.1017/nws.2025.10018)

If `h` has `Composite Types` either for vertex metadata or hyperedges metadata,
the user has to explicit tell the JSON package about it, for instance using:

TODO: this

See the (JSON.jl documentation)[https://juliaio.github.io/JSON.jl/stable/] for more details.
"""
function SimpleHypergraphs.hg_save(
    io::IO,
    h::DirectedHypergraph{T, V, E, D},
    format::HIF_Format;
    pretty::Bool = false
) where {T, V, E, D}
    incidences = Vector{OrderedDict{Symbol, Union{Int, T}}}()

    for i in 1:nhv(h)
	hes = gethyperedges(h, i)
	# Tails
	for j in sort!(collect(keys(hes[1])))
            push!(incidences, OrderedDict{Symbol, Union{Int, T}}(
		:node => i,
		:edge => j,
		:weight => T(h[1, i, j]),
		:direction => "tail"
	    ))
        end

	# Heads
	for j in sort!(collect(keys(hes[2])))
            push!(incidences, OrderedDict{Symbol, Union{Int, T}}(
		:node => i,
		:edge => j,
		:weight => T(h[2, i, j]),
		:direction => "head"
	    ))
        end
    end

    # Decide whether to include metadata for nodes and edges
    # There are two poossible reasons to include metadata:
    #	1. there is at least one metadata entry
    #	2. there is at least one node or edge with no connections (isolated vertex or empty hyperedge)
    node_meta_included = any(x -> !(isnothing(x)), h.v_meta) || any(
	v -> isempty(h.hg_tail.v2he[v]) && isempty(h.hg_head.v2he[v]
    ), 1:nhv(h))
    
    edge_meta_included = (
	any(x -> !(isnothing(x)), h.he_meta_tail) 
	|| any(x -> !(isnothing(x)), h.he_meta_head)
	|| any(e -> isempty(h.hg_tail.he2v[e]) && isempty(h.hg_head.he2v[e]), 1:nhe(h))
    )

    json_node_meta = Vector{OrderedDict{Symbol, Any}}()
    json_edge_meta = Vector{OrderedDict{Symbol, Any}}()
    
    if node_meta_included
        for i in 1:nhv(h)
            node_entry = OrderedDict{Symbol, Union{Int, typeof(h.v_meta[i])}}(:node => i)
            if !(isnothing(h.v_meta[i]))
                node_entry[:attrs] = h.v_meta[i]
            end
            push!(json_node_meta, node_entry)
        end
    end
    if edge_meta_included
        for j in 1:nhe(h)
	    edge_entry = OrderedDict{
		Symbol, 
		Union{
		    Int,
		    Vector{Union{typeof(h.hg_meta_tail[j]), typeof(h.hg_meta_head[j])}}
		}
	    }(:edge => j)
	    if !(isnothing(h.he_meta_tail[j]) && isnothing(h.he_meta_head[j]))
		edge_entry[:attrs] = [h.he_meta_tail[j], h.he_meta_head[j]]
            end
            push!(json_edge_meta, edge_entry)
        end
    end

    json_dhg = OrderedDict{
	Symbol,
	Union{typeof(incidences), typeof(json_node_meta), typeof(json_edge_meta)}
    }()

    json_dhg[:incidences] = incidences
    if length(json_node_meta) > 0
	json_dhg[:nodes] = json_node_meta
    end
    if length(json_edge_meta) > 0
	json_dhg[:edges] = json_edge_meta
    end
    
    JSON.json(io, json_dhg; pretty)
end



"""
    dhg_load(
        io::IO,
        format::EHGF_Format;
        HType::Type{H} = DirectedHypergraph,
        T::Type{U} = Bool,
        D::Type{<:AbstractDict{Int, U}} = Dict{Int, T},
    ) where {U <: Real, H <: AbstractDirectedHypergraph}

Loads a hypergraph from a stream `io` from `ehgf` format.

**Arguments**

* `T` : type of weight values stored in the hypergraph's adjacency matrix
* `D` : dictionary for storing values the default is `Dict{Int, T}`

Skips a single initial comment.

"""
function dhg_load(
    io::IO,
    format::EHGF_Format;
    HType::Type{H} = DirectedHypergraph,
    T::Type{U} = Bool,
    D::Type{<:AbstractDict{Int, U}} = Dict{Int, T},
) where {U <: Real, H <: AbstractDirectedHypergraph}
    line = readline(io)

    if startswith(line, "\"\"\"")
      singleline = true
        while(
            !( (!singleline && endswith(line, "\"\"\"")) ||
            (singleline && endswith(line, "\"\"\"") && length(line)>5)
            ) &&
            !eof(io)
            )
                line = readline(io)
                singleline = false
        end
        if eof(io)
            throw(ArgumentError("malformed input"))
        end
       line = readline(io)
    end

    l = split(strip(line))
    length(l) == 2 || throw(ArgumentError("expected two integers"))
    n, k = parse.(Int, l)
    h = HType{T, D}(n, k)

    for i in 1:k
        lastv = 0

        ht = split(readline(io), " || ")
        length(ht) == 2 || throw(ArgumentError("Expected one head and one tail!"))

        he_tail, he_head = ht
        
        for pos in split.(strip.(he_tail))
            entry = split(pos, '=')
            length(entry) == 2 || throw(ArgumentError("Expected format: vertex=weight"))

            v = parse(Int, entry[1])
            w = parse(T, entry[2])

            if v > lastv
                lastv = v
            else
                throw(ArgumentError("Vertices in hyperedge must be sorted!"))
            end
            h.hg_tail[v, i] = w
        end

        lastv = 0
        for pos in split.(strip.(he_head))
            entry = split(pos, '=')
            length(entry) == 2 || throw(ArgumentError("Expected format: vertex=weight"))

            v = parse(Int, entry[1])
            w = parse(T, entry[2])

            if v > lastv
                lastv = v
            else
                throw(ArgumentError("Vertices in hyperedge must be sorted!"))
            end
            h.hg_head[v, i] = w
        end

    end
    # we ignore lines beyond k+1 in the file
    h
end


"""
    dhg_load(
        io::IO,
        T::Type{H},
        format::JSON_Format;
        T::Type{U} = Bool,
        D::Type{<:AbstractDict{Int, U}} = Dict{Int,U},
        V = Nothing,
        E = Nothing
    ) where {H <: AbstractDirectedHypergraph, U <: Real}

Loads a hypergraph from a stream `io` from `json` format.

**Arguments**

* `T` : type of weight values stored in the hypergraph's adjacency matrix
* `D` : dictionary for storing values the default is `Dict{Int, T}`
* `V` : type of values stored in the vertices of the hypergraph
* `E` : type of values stored in the edges of the hypergraph

"""
function dhg_load(
        io::IO,
        format::JSON_Format;
        HType::Type{H} = DirectedHypergraph,
        T::Type{U} = Bool,
        D::Type{<:AbstractDict{Int, U}} = Dict{Int, T},
        V = Nothing,
        E = Nothing
    ) where {U <: Real, H<:AbstractDirectedHypergraph}
    json_hg = JSON.parse(readline(io))

    m_tail = reshape(JSON.parse(json_hg.tail, Array{Union{T, Nothing}}), json_hg.n, json_hg.k)
    m_head = reshape(JSON.parse(json_hg.head, Array{Union{T, Nothing}}), json_hg.n, json_hg.k)

    if V != Nothing
        v_meta = JSON.parse(json_hg.v_meta, Array{Union{V, Nothing}})
    else
	v_meta = nothing
    end

    if E != Nothing
        he_meta_tail = JSON.parse(json_hg.he_meta_tail, Array{Union{E, Nothing}})
        he_meta_head = JSON.parse(json_hg.he_meta_head, Array{Union{E, Nothing}})
    else
	he_meta_tail = nothing
	he_meta_head = nothing
    end

    HType{T, V, E, D}(m_tail, m_head; v_meta=v_meta, he_meta_tail=he_meta_tail, he_meta_head=he_meta_head)
end

"""
    dhg_load(
	io::IO,
	format::HIF_Format;
	HType::Type{H} = DirectedHypergraph,
	T::Type{U} = Bool,
	D::Type{<:AbstractDict{Int, U}} = Dict{Int, T},
	V::Union{Type, Symbol} = :auto,
	E::Union{Type, Symbol} = :auto,
	sort_by_id::Bool = false,
	add_original_id_to_meta::Union{Symbol, Nothing} = nothing,
	show_warning::Bool = true
    )

Loads a directed hypergraph from an input stream `io` in `HIF` format (Coll et al.,
DOI: 10.1017/nws.2025.10018).

`T` is the type of weight values stored in the directed hypergraph matrix's matrix representation,
and `D` is the dictionary type used to store weights for vertices in directed hyperedges. If the
directed hypergraph has vertex or hyperedge metadata, their types can be specified using `V` and
`E`, respectively.

Vertex and hyperedge indices are regenerated to match the 1-based indexing scheme of
`SimpleDirectedHypergraphs.jl`. The original indices can be preserved by setting
`add_original_id_to_meta` to a Symbol representing the key under which those indices will be stored
in the metadata dictionary.
"""
function dhg_load(
    io::IO,
    format::HIF_Format;
    HType::Type{H} = DirectedHypergraph,
    T::Type{U} = Bool,
    D::Type{<:AbstractDict{Int, U}} = Dict{Int, T},
    V::Union{Type, Symbol} = :auto,
    E::Union{Type, Symbol} = :auto,
    sort_by_id::Bool = false,
    add_original_id_to_meta::Union{Symbol, Nothing} = nothing,
    show_warning::Bool = true,
    validate_schema::Bool = true
) 
    data = JSON.parse(io; dicttype=Dict{Symbol, Any})

    if validate_schema
	schema_url = "https://raw.githubusercontent.com/pszufe/HIF-standard/main/schemas/hif_schema.json"
	schema = String(HTTP.get(schema_url).body)
	validator = Schema(schema)
	@assert JSONSchema.validate(validator, data) "Failed HIF schema validation!"
    else
	# More basic, minimal validation
	haskey(data, :incidences) || throw(ArgumentError("Missing required attribute 'incidences'"))
    end

    if isempty(data[:incidences])
        if isempty(get(data, :edges, [])) && isempty(get(data, :nodes, []))
            return Hypergraph{
                T, 
                V == :auto ? Nothing : V, 
                E == :auto ? Nothing : E,
                D,
            }(0, 0)
        end
    end

    nodesdf = _build_attr_dataframe(data, :nodes, V, add_original_id_to_meta)
    edgesdf = _build_attr_dataframe(data, :edges, E, add_original_id_to_meta)

    attr_nodes_N = nrow(nodesdf)  
    attr_edges_N = nrow(edgesdf)
    if attr_nodes_N == 0 && isnothing(add_original_id_to_meta)
        # no node attributes found so all attrs set to Nothing
        nodesdf.attrs = Nothing[]
    end
    if attr_edges_N == 0 && isnothing(add_original_id_to_meta)
        # no edge attributes found so all attrs set to Nothing
        edgesdf.attrs = Nothing[]
    end

    _add_nodes_and_edges_from_incidences!(data, nodesdf, edgesdf, add_original_id_to_meta)

    # narrow types for attrs if V or E is :auto
    if V == :auto
        _sanitize_types_items!(nodesdf)
    end
    if E == :auto
        _sanitize_types_items!(edgesdf)
    end

    _add_id_sort_column!(nodesdf)
    _add_id_sort_column!(edgesdf)
    # if all nodes or edges were discovered from incidences, sort by id to have consistent ordering
    if attr_nodes_N == 0
        sort!(nodesdf, :id_sort)
    end
    if attr_edges_N == 0
        sort!(edgesdf, :id_sort)
    end

    if sort_by_id
        sort!(nodesdf, :id_sort)
        sort!(edgesdf, :id_sort)
    end

    if show_warning
        if nrow(nodesdf) > 0 && nodesdf.id != 1:nrow(nodesdf)
            @warn "Nodes in the source file were not sorted or not consistent - their order will change"
        end

        if nrow(edgesdf) > 0 && edgesdf.id != 1:nrow(edgesdf)
            @warn "Edges in the source file were not sorted or not consistent - their order will change"
        end
    end

    hg = HType{
        T, 
        eltype(nodesdf.attrs), 
        eltype(edgesdf.attrs), 
        D,
    }(nrow(nodesdf), nrow(edgesdf), nodesdf.attrs, edgesdf.attrs)

    _add_weights_from_incidences!(data, hg, edgesdf, nodesdf)
    hg    
end

"""
    dhg_load(
        fname::AbstractString;
        format::Abstract_HG_format = HGF_Format(),
        HType::Type{H} = DirectedHypergraph,
        T::Type{U} = Bool,
        D::Type{<:AbstractDict{Int, U}} = Dict{Int, T},
        V = Nothing,
        E = Nothing
    ) where {U <: Real, H <: AbstractDirectedHypergraph}

Loads a hypergraph from a file `fname`.
The default saving format is `hgf`.

**Arguments**

* `HType`: type of hypergraph to store data in
* `T` : type of weight values stored in the hypergraph's adjacency matrix
* `D` : dictionary for storing values the default is `Dict{Int, T}`
* `V` : type of values stored in the vertices of the hypergraph
* `E` : type of values stored in the edges of the hypergraph

"""
function dhg_load(
        fname::AbstractString;
        format::Abstract_HG_format = EHGF_Format(),
        HType::Type{H} = DirectedHypergraph,
        T::Type{U} = Bool,
        D::Type{<:AbstractDict{Int, U}} = Dict{Int, T},
        V = Nothing,
        E = Nothing
    ) where {U <: Real, H <: AbstractDirectedHypergraph}

    if format == EHGF_Format()
        if HType == DirectedHypergraph
            open(io -> dhg_load(io, format; HType=HType, T=T, D=D), fname, "r")
        else
            error("EHGF loading only implemented for DirectedHypergraph")
        end
    else
        open(io -> dhg_load(io, format; HType=HType, T=T, D=D, V=V, E=E), fname, "r")
    end

end
