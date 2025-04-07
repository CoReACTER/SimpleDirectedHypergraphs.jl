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
    hg_save(io::IO, h::DirectedHypergraph, format::JSON_Format)

Saves a directed hypergraph `h` to an output stream `io` in `json` format.

If `h` has `Composite Types` either for vertex metadata or hyperedges metadata,
the user has to explicit tell the JSON3 package about it, for instance using:

`JSON3.StructType(::Type{MyType}) = JSON3.Struct()`.

See the (JSON3.jl documentation)[https://github.com/quinnj/JSON3.jl] for more details.

The `json` in output contains the following information (keys):

* `n` : number of vertices
* `k` : number of hyperedges
* `tail` : a matrix representation of the tails of `h`, where rows are vertices and columns are hyperedges
* `head` : a matrix representation of the heads of `h`, where rows are vertices and columns are hyperedges
* `v_meta` : vertices metadata
* `he_meta_tail` : metadata for hyperedge tails
* `he_meta_head` : metadata for hyperedge heads

"""
function SimpleHypergraphs.hg_save(io::IO, h::DirectedHypergraph, format::JSON_Format)
    json_hg = Dict{Symbol, Any}()

    json_hg[:n] = nhv(h)
    json_hg[:k] = nhe(h)

    json_hg[:tail] = JSON3.write(Matrix(h.hg_tail))
    json_hg[:head] = JSON3.write(Matrix(h.hg_head))
    
    json_hg[:v_meta] = JSON3.write(h.v_meta)
    json_hg[:he_meta_tail] = JSON3.write(h.he_meta_tail)
    json_hg[:he_meta_head] = JSON3.write(h.he_meta_head)

    JSON3.write(io, json_hg)
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
    json_hg = JSON3.read(readline(io))

    m_tail = reshape(JSON3.read(json_hg.tail, Array{Union{T, Nothing}}), json_hg.n, json_hg.k)
    m_head = reshape(JSON3.read(json_hg.head, Array{Union{T, Nothing}}), json_hg.n, json_hg.k)

    if V != Nothing && E != Nothing
        v_meta = JSON3.read(json_hg.v_meta, Array{Union{V, Nothing}})
        he_meta_tail = JSON3.read(json_hg.he_meta_tail, Array{Union{E, Nothing}})
        he_meta_head = JSON3.read(json_hg.he_meta_head, Array{Union{E, Nothing}})
        h = HType{T, V, E, D}(m_tail, m_head; v_meta=v_meta, he_meta_tail=he_meta_tail, he_meta_head=he_meta_head)
    else
        h = HType{T, D}(m_tail, m_head)
    end

    h
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
