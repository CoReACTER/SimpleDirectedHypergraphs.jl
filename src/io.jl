# TODO: maybe more fancy file format and correctness checking should be done

struct EHGF_Format <: Abstract_HG_format end


"""
    hg_save(io::IO, h::H, format::EHGF_Format; pretty::Bool = false) where {H <: AbstractDirectedHypergraph}

Saves an undirected hypergraph `h` to an output stream `io` in `ehgf` format.

**Arguments**

* `io`: input-output stream object
* `h`: directed hypergraph to be saved
* `format`: file format (here `EHGF_Format`, for the "extended hypergraph format")
* `pretty`: should this be pretty-printed? (NOTE: for EHGF output, this argument does nothing, but
    it's still required by the `SimpleHypergraphs.jl` API)

"""
function SimpleHypergraphs.hg_save(io::IO, h::H, format::EHGF_Format; pretty::Bool = false) where {H <: AbstractDirectedHypergraph}

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
    return
end


"""
    hg_save(
	io::IO,
	h::DirectedHypergraph{T, V, E, D},
	::HIF_Format;
	pretty::Bool=false
    ) where {T, V, E, D}

Saves a directed hypergraph `h` to an output stream `io` in `HIF` format (see Coll et al.,
DOI: 10.1017/nws.2025.10018).

TODO: handling for composite metadata types

**Arguments**

* `io`: input-output stream object
* `h`: directed hypergraph to be saved
* `format`: file format (here `HIF_Format`, for the "hypergraph interchange format")
* `pretty`: should this be pretty-printed?

**Returns**

`JSON.json` output in `HIF` format

"""
function SimpleHypergraphs.hg_save(
        io::IO,
        h::DirectedHypergraph{T, V, E, D},
        format::HIF_Format;
        pretty::Bool = false
    ) where {T, V, E, D}
    incidences = Vector{OrderedDict{Symbol, Union{Int, String, T}}}()

    for i in 1:nhv(h)
        hes = gethyperedges(h, i)
        # Tails
        for j in sort!(collect(keys(hes[1])))
            push!(
                incidences, OrderedDict{Symbol, Union{Int, String, T}}(
                    :node => i,
                    :edge => j,
                    :weight => T(h[i, j][1]),
                    :direction => "tail"
                )
            )
        end

        # Heads
        for j in sort!(collect(keys(hes[2])))
            push!(
                incidences, OrderedDict{Symbol, Union{Int, String, T}}(
                    :node => i,
                    :edge => j,
                    :weight => T(h[i, j][2]),
                    :direction => "head"
                )
            )
        end
    end

    # Decide whether to include metadata for nodes and edges
    # There are two poossible reasons to include metadata:
    #	1. there is at least one metadata entry
    #	2. there is at least one node or edge with no connections (isolated vertex or empty hyperedge)
    node_meta_included = any(x -> !(isnothing(x)), h.v_meta) || any(
        v -> isempty(h.hg_tail.v2he[v]) && isempty(
            h.hg_head.v2he[v]
        ), 1:nhv(h)
    )

    edge_meta_included = (
        any(x -> !(isnothing(x)), h.he_meta_tail)
            || any(x -> !(isnothing(x)), h.he_meta_head)
            || any(e -> isempty(h.hg_tail.he2v[e]) && isempty(h.hg_head.he2v[e]), 1:nhe(h))
    )

    json_node_meta = Vector{OrderedDict{Symbol, Any}}()
    json_edge_meta = Vector{OrderedDict{Symbol, Any}}()

    if node_meta_included
        for i in 1:nhv(h)
            node_entry = OrderedDict{
                Symbol,
                Union{Int, OrderedDict{Symbol, typeof(h.v_meta[i])}},
            }(:node => i)
            if !(isnothing(h.v_meta[i]))
                node_entry[:attrs] = OrderedDict{Symbol, typeof(h.v_meta[i])}(:value => h.v_meta[i])
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
                    OrderedDict{Symbol, Union{typeof(h.he_meta_tail[j]), typeof(h.he_meta_head[j])}},
                },
            }(:edge => j)
            if !(isnothing(h.he_meta_tail[j]) && isnothing(h.he_meta_head[j]))
                edge_entry[:attrs] = OrderedDict{Symbol, Union{typeof(h.he_meta_tail[j]), typeof(h.he_meta_head[j])}}(
                    :tail => h.he_meta_tail[j],
                    :head => h.he_meta_head[j]
                )
            end
            push!(json_edge_meta, edge_entry)
        end
    end

    # Make Dict to be output as JSON
    json_dhg = OrderedDict{
        Symbol,
        Union{typeof(incidences), typeof(json_node_meta), typeof(json_edge_meta)},
    }()

    json_dhg[:incidences] = incidences
    if length(json_node_meta) > 0
        json_dhg[:nodes] = json_node_meta
    end
    if length(json_edge_meta) > 0
        json_dhg[:edges] = json_edge_meta
    end

    return JSON.json(io, json_dhg; pretty)
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

Skips a single initial comment.

**Arguments**

* `io` : input-output stream object form which to load a directed hypergraph
* `format` : file format from which to load (here, `EHGF_Format`, or the "extended hypergraph
    format")
* `HType` : directed hypergraph datatype to store data from file (default is `DirectedHypergraph`)
* `T` : real-valued datatype for incidence matrix entries/weights (default is `Bool`)
* `D` : dictionary datatype for internal storage (default is `Dict{Int, T}`)

**Returns**

A dihypergraph of type `HType{T, Any, Any, D}`

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
        while (
                !(
                    (!singleline && endswith(line, "\"\"\"")) ||
                        (singleline && endswith(line, "\"\"\"") && length(line) > 5)
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
    h = HType{T, Any, Any, D}(n, k)

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
    return h
end


"""
    _add_weights_from_incidences!(
	data::Dict{Symbol, Any},
	hg::AbstractDirectedHypergraph{Tuple{Union{T,Nothing},Union{T,Nothing}}},
	edges::DataFrame,
	nodes::DataFrame
    )

THIS FUNCTION IS INTERNAL AND SHOULD NOT BE CALLED DIRECTLY.
Adds weights to the hypergraph `hg` based on the incidences provided in `data`.
The `edges` and `nodes` DataFrames are used to map edge and node identifiers to
their respective indices in the hypergraph.

This code is heavily inspired by the corresponding function in `SimpleHypergraphs.jl`; see
https://github.com/pszufe/SimpleHypergraphs.jl.

**Arguments**

* `data` : parsed `HIF` (hypergraph interchange format) data in dictionary format
* `hg` : directed hypergraph object to be modified with incidence weights
* `edges` : `DataFrame` with directed hyperedge metadata
* `nodes` : `DataFrame` with vertex metadata

"""
function _add_weights_from_incidences!(
        data::Dict{Symbol, Any},
        hg::AbstractDirectedHypergraph{Tuple{Union{T, Nothing}, Union{T, Nothing}}},
        edges::DataFrame,
        nodes::DataFrame
    ) where {T <: Real}
    node_dict = Dict{Union{String, Int}, Int}(id => idx for (id, idx) in zip(nodes.id, 1:nrow(nodes)))
    edge_dict = Dict{Union{String, Int}, Int}(id => idx for (id, idx) in zip(edges.id, 1:nrow(edges)))

    incidences = data[:incidences]

    for incidence in incidences
        node_idx = node_dict[incidence[:node]]
        edge_idx = edge_dict[incidence[:edge]]
        weight = get(incidence, :weight, one(T))
        direction = get(incidence, :direction, "none")

        if direction == "none"
            @warn "No direction given for hyperedge incidence. Ignoring; cannot include weight."
        else
            side = direction == "tail" ? 1 : 2
            hg[side, node_idx, edge_idx] = T(weight)
        end
    end
    return
end

"""
    _build_attr_dataframe(
	data::Dict{Symbol, Any},
	field::Symbol,
	V::Union{Type, Symbol},
	add_original_id_to_meta::Union{Symbol, Nothing},
	to_convert::Union{Nothing, Symbol, AbstractVector{Symbol}}
    )

THIS FUNCTION IS INTERNAL AND SHOULD NOT BE CALLED DIRECTLY.
Constructs vertex or directed hyperedge `DataFrame` objects based on attributes parsed from an
`HIF` (hypergraph interchange format) file.

This code is heavily inspired by the corresponding function in `SimpleHypergraphs.jl`; see
https://github.com/pszufe/SimpleHypergraphs.jl.

**Arguments**

* `data` : parsed `HIF` data in dictionary format
* `field` : category of data to be processed; either `:nodes` or `:edges`
* `V` : metadata type (or `:auto`), used to convert parsed data to appropriate type
* `add_original_id_to_meta` : ID field to be included in dataframe (to account for
    `SimpleDirectedHypergraphs` assuming sequential 1-based indexing for vertices/directed
    hyperedges)
* `to_convert` : symbols from the `HIF` attributes to be converted and included in the final
    dataframe

**Returns**

`items`, a `DataFrame` with relevant (meta)data

"""
function _build_attr_dataframe(
        data::Dict{Symbol, Any},
        field::Symbol,
        V::Union{Type, Symbol},
        add_original_id_to_meta::Union{Symbol, Nothing},
        to_convert::Union{Nothing, Symbol, AbstractVector{Symbol}}
    )
    @assert field ∈ (:nodes, :edges)
    fid = Symbol(string(field)[1:(end - 1)])  # :node or :edge

    # Make list of symbols to type-convert
    if isnothing(to_convert)
        to_convert = Symbol[]
    elseif isa(to_convert, Symbol)
        to_convert = [to_convert]
    end

    target_attr_type = Any
    dict_type = Dict{Symbol, Any}
    if V != :auto
        if isnothing(add_original_id_to_meta) && length(to_convert) <= 1
            target_attr_type = V
        elseif isnothing(add_original_id_to_meta) && length(to_convert) > 1
            # Vertex/dihyperedge attributes will (at least initially) be dict if the input type is
            # a dict in need of conversion to the type `V`
            target_attr_type = Dict{Symbol, V}
            dict_type = target_attr_type
        else
            # Including ID type
            target_attr_type = Dict{Symbol, Union{V, Int, String}}
            dict_type = target_attr_type
        end
    end

    items = DataFrame(;
        id = Union{String, Int}[],
        attrs = target_attr_type[]
    )
    if !haskey(data, field)
        return items
    end

    seen = Set{Union{Int, String}}()
    for item in data[field]
        id = item[fid]
        if id ∈ seen
            continue
        end
        val = get(item, :attrs, nothing)

        # Symbol-value pairs
        # This works if val is a single metadata value or (if using `hg_save`) a dictionary
        kvs = Tuple{Symbol, Any}[]
        if length(to_convert) == 0
            push!(kvs, (:only, val))
        else
            for tc in to_convert
                push!(kvs, (tc, get(val, tc, nothing)))
            end
        end

        # Gather up everything that's going to be a metadata "value" for this vertex/dihyperedge
        for_attrs = Tuple{Symbol, Any}[]
        for (k, v) in kvs
            if isnothing(v)
                continue
            end

            if V == String && !(isa(v, String))
                push!(for_attrs, (k, JSON.json(v)))
            elseif V != :auto
                push!(for_attrs, (k, convert(V, v)))
            else
                push!(for_attrs, (k, v))
            end
        end

        if !isnothing(add_original_id_to_meta)
            push!(for_attrs, (add_original_id_to_meta, id))
        end

        if length(for_attrs) == 1 && for_attrs[1][1] in [:only, :value]
            val = for_attrs[1][2]
        else
            val = dict_type()
            for (k, v) in for_attrs
                val[k] = v
            end
        end

        push!(items, [id, val])
        push!(seen, id)
    end

    return items
end

"""
    _separate_tail_head_meta(data::Vector{X}) where {X}

THIS FUNCTION IS INTERNAL AND SHOULD NOT BE CALLED DIRECTLY.
If possible, splits the edge `DataFrame` (generated by `_build_attr_dataframe`) into a collection
for tail metadata and a collection for head metadata

TODO: maybe this is too storage-inefficient? Maybe we put it all on one side?
Maybe we reassess if having split tail and head metadata is a good idea?

This code is heavily inspired by the corresponding function in `SimpleHypergraphs.jl`; see
https://github.com/pszufe/SimpleHypergraphs.jl.

**Arguments**

* `data` : a vector of metadata entries, each of which could be a dictionary, a 2-tuple, or a
    single value

**Returns**

A 2-tuple with the first entry containing dihyperedge tail metadata and the second entry containin
dihyperedge head metadata

"""
function _separate_tail_head_meta(data::Vector{X}) where {X}
    # Simplest case - "tail" and "head" as separate columns
    if X <: AbstractDict
        if !(all(x -> (:tail in keys(x) && :head in keys(x)), data))
            # No clearly marked tail and head data
            return (data, data)
        else

            t = X[]
            h = X[]

            for d in data
                thist = X()
                thish = X()

                for (k, v) in d
                    if k != :tail
                        thish[k] = v
                    end
                    if k != :head
                        thist[k] = v
                    end
                end
                push!(t, thist)
                push!(h, thish)
            end

            return (t, h)
        end
    elseif X <: NTuple{2, Any}
        # If each edge has two values, assume these are the tail and head metadata
        return ([e[1] for e in data], [e[2] for e in data])
    else
        # If there aren't multiple keys, then assume that data is overall tail and head metadata
        # Cannot separate
        return (data, data)
    end
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
    ) where {H <: AbstractDirectedHypergraph, U <: Real}

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

**Arguments**

* `io` : input-output stream to be parsed
* `format` : file format to be parsed (here, `HIF_Format`)
* `HType` : directed hypergraph type (default is `DirectedHypergraph`)
* `T` : real-valued data type for incidence weights (default is `Bool`)
* `D` : dictionary data type for internal storage (default is `Dict{Int, T}`)
* `V` : vertex metadata type, or `:auto` (default is `:auto`)
* `E` : hyperedge metadata type, or `:auto` (default is `:auto`)
* `sort_by_id` : if `true` (default `false`), then the vertex and hyperedge `DataFrames` generated
    during parsing will be sorted by original ID labels
* `add_original_id_to_meta` : ID field to be included in dataframe (to account for
    `SimpleDirectedHypergraphs` assuming sequential 1-based indexing for vertices/directed
    hyperedges)
* `show_warning` : If `true` (default), a warning will be printed if vertex or edge IDs are changed
    during parsing
* `to_convert` : symbols from the `HIF` attributes to be converted and included in the final
    dataframe

**Returns**

A dihypergraph of type `HType` with weight type `T`, dict type `D`, and vertex and dihyperedge type
dependent on `V`, `E`, and the `DataFrames` resulting from parsing and processing

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
        to_convert::Union{Nothing, Symbol, AbstractVector{Symbol}} = :auto
    ) where {H <: AbstractDirectedHypergraph, U <: Real}
    data = JSON.parse(io; dicttype = Dict{Symbol, Any})

    # Basic, minimal validation
    # "incidences" is the only required key in the HIF standard
    haskey(data, :incidences) || throw(ArgumentError("Missing required attribute 'incidences'"))

    if isempty(data[:incidences])
        if isempty(get(data, :edges, [])) && isempty(get(data, :nodes, []))
            return HType{
                T,
                V == :auto ? Nothing : V,
                E == :auto ? Nothing : E,
                D,
            }(0, 0)
        end
    end

    nodesdf = _build_attr_dataframe(
        data,
        :nodes,
        V,
        add_original_id_to_meta,
        (to_convert == :auto ? :value : to_convert)
    )
    edgesdf = _build_attr_dataframe(
        data,
        :edges,
        E,
        add_original_id_to_meta,
        (to_convert == :auto ? [:tail, :head] : to_convert)
    )

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

    SimpleHypergraphs._add_nodes_and_edges_from_incidences!(data, nodesdf, edgesdf, add_original_id_to_meta)

    # narrow types for attrs if V or E is :auto
    if V == :auto
        SimpleHypergraphs._sanitize_types_items!(nodesdf)
    end
    if E == :auto
        SimpleHypergraphs._sanitize_types_items!(edgesdf)
    end

    SimpleHypergraphs._add_id_sort_column!(nodesdf)
    SimpleHypergraphs._add_id_sort_column!(edgesdf)

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

    tail_meta, head_meta = _separate_tail_head_meta(edgesdf.attrs)

    hg = HType{T, eltype(nodesdf.attrs), eltype(tail_meta), D}(
        nrow(nodesdf), nrow(edgesdf), nodesdf.attrs, tail_meta, head_meta
    )

    _add_weights_from_incidences!(data, hg, edgesdf, nodesdf)

    return hg
end


"""
    dhg_load(
        fname::AbstractString;
        format::Abstract_HG_format = EHGF_Format(),
        HType::Type{H} = DirectedHypergraph,
        T::Type{U} = Bool,
        D::Type{<:AbstractDict{Int, U}} = Dict{Int, T},
        V = Nothing,
        E = Nothing
    ) where {U <: Real, H <: AbstractDirectedHypergraph}

Loads a hypergraph from a file `fname`.
The default saving format is `ehgf`.

**Arguments**

* `fname` : filename
* `format` : file format (default is `EHGF_Format`, the "extended hypergraph format")
* `HType` : type of hypergraph to store data in
* `T` : type of weight values stored in the hypergraph's adjacency matrix
* `D` : dictionary for storing values the default is `Dict{Int, T}`
* `V` : type of values stored in the vertices of the hypergraph
* `E` : type of values stored in the edges of the hypergraph

**Returns**

A dihypergraph (of type `HType` in general, or of type `DirectedHypergraph` for `EHGF`-formatted
inputs), with weight type `T`, dictionary type `D`, and vertex and metadata types `V` and `E`,
respectively.

Note that, for parsing from `HIF`-formatted data, the metadata types may change

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

    return if format == EHGF_Format()
        if HType == DirectedHypergraph
            open(io -> dhg_load(io, format; HType = HType, T = T, D = D), fname, "r")
        else
            error("EHGF loading only implemented for DirectedHypergraph")
        end
    else
        open(io -> dhg_load(io, format; HType = HType, T = T, D = D, V = V, E = E), fname, "r")
    end

end
