using SimpleDirectedHypergraphs
using SimpleHypergraphs
using StatsBase
using Random
using DataStructures
using Graphs
using Test
# using JET

dh1 = DirectedHypergraph{Float64,Int,String}(7, 6)
dh1[1, 1, 1] = 1.0
dh1.hg_head[2:3, 1] .= 2.5  # Assignment on a directed hypergraph directly with slices is currently awkward
dh1[1, 3, 3] = 4.0
dh1[2, 4, 3] = 5.5
dh1[1, 5, 4] = 7.0
dh1[2, 6, 4] = 8.5
dh1[1, 6, 5] = 10.0
dh1[2, 7, 5] = -1.5
dh1[1, 7, 6] = 0.0
dh1[2, 5, 6] = 1.5

tail_2 = [
    true nothing nothing nothing nothing nothing nothing nothing nothing
    nothing true nothing nothing nothing true nothing nothing nothing
    nothing nothing true nothing nothing nothing true nothing nothing
    nothing nothing nothing true nothing true true nothing nothing
    nothing nothing nothing nothing true nothing nothing true true
]
head_2 = [
    nothing nothing nothing nothing nothing true true nothing nothing
    true nothing nothing nothing nothing nothing nothing true nothing
    true nothing nothing nothing nothing nothing nothing nothing true
    nothing true true nothing true nothing nothing nothing nothing
    nothing nothing nothing true nothing nothing nothing nothing nothing
]
dh2 = DirectedHypergraph(tail_2, head_2)

# @testset "SimpleDirectedHypergraphs Code linting (JET.jl)" begin
#     JET.test_package(SimpleDirectedHypergraphs; target_defined_modules = true)
# end


@testset "SimpleDirectedHypergraphs DirectedHypergraph     " begin
    h = dhg_load("data/test_dhg.ehgf"; format=EHGF_Format(), T=Int, HType=DirectedHypergraph)
    @test size(h) == (6, 3)
    @test nhv(h) == 6
    @test nhe(h) == 3
    m = Matrix(h)
    @test m == h
    @test h == [
        (1, nothing) (nothing, nothing) (nothing, nothing)
        (2, nothing) (3, nothing) (nothing, nothing)
        (nothing, nothing) (nothing, 0) (nothing, nothing)
        (nothing, 4) (nothing, nothing) (1, nothing)
        (nothing, 5) (12, nothing) (nothing, nothing)
        (nothing, nothing) (nothing, nothing) (nothing, 4)
    ]
    mktemp("data") do path, _
        println(path)
        SimpleHypergraphs.hg_save(path, h; format=EHGF_Format())

        loaded_hg = replace(read(path, String), r"\n*$" => "")

        @test loaded_hg ==
              reduce(replace,
            ["\r\n" => "\n",
                r"^\"\"\"(?s).*\"\"\"\n" => "", #remove initial comments
                r"\n*$" => ""], #remove final \n*
            init=read("data/test_dhg.ehgf", String)) #no comments

        @test loaded_hg ==
              reduce(replace,
            ["\r\n" => "\n",
                r"^\"\"\"(?s).*\"\"\"\n" => "", #remove initial comments
                r"\n*$" => ""], #remove final \n*
            init=read("data/singlelinecomment.ehgf", String)) #single line comment

        @test loaded_hg ==
              reduce(replace,
            ["\r\n" => "\n",
                r"^\"\"\"(?s).*\"\"\"\n" => "", #remove initial comments
                r"\n*$" => ""], #remove final \n*
            init=read("data/multilinecomment.ehgf", String)) #multiple lines comment

        for v = 1:nhv(dh1)
            set_vertex_meta!(dh1, v, v)
        end

        for he = 1:nhe(dh1)
            set_hyperedge_meta!(dh1, string(he), string(he), he)
        end

        SimpleHypergraphs.hg_save(path, dh1; format=JSON_Format())
        loaded_hg = dhg_load(path; format=JSON_Format(), HType=DirectedHypergraph, T=Float64, V=Int, E=String)

        @test dh1 == loaded_hg
        @test dh1.v_meta == loaded_hg.v_meta
        @test dh1.he_meta_tail == loaded_hg.he_meta_tail
        @test dh1.he_meta_head == loaded_hg.he_meta_head

        @test get_vertex_meta(dh1, 1) == get_vertex_meta(loaded_hg, 1)
        @test get_hyperedge_meta(dh1, 2) == get_hyperedge_meta(loaded_hg, 2)

    end

    @test_throws ArgumentError dhg_load("data/malformedcomment.ehgf"; format=EHGF_Format(), HType=DirectedHypergraph, T=Int)
    @test_throws ArgumentError dhg_load("data/argumenterror.ehgf"; format=EHGF_Format(), HType=DirectedHypergraph, T=Int)

    dh2 = DirectedHypergraph{Float64}(0, 0)
    @test dh2 == DirectedHypergraph{Float64,Nothing}(0, 0)
    @test dh2 == DirectedHypergraph{Float64,Nothing,Nothing}(0, 0)
    @test dh2 == DirectedHypergraph{Float64,Nothing,Nothing,Dict{Int,Float64}}(0, 0)

    dh3 = DirectedHypergraph(0, 0)
    @test dh3 == DirectedHypergraph{Bool,Nothing,Nothing,Dict{Int,Bool}}(0, 0)

    for i in 1:6
        SimpleHypergraphs.add_vertex!(dh2)
    end
    SimpleHypergraphs.add_hyperedge!(dh2; vertices_tail=Dict(1 => 1.0), vertices_head=Dict(2:3 .=> 2.5))
    SimpleHypergraphs.add_hyperedge!(dh2)
    SimpleHypergraphs.add_hyperedge!(dh2; vertices_tail=Dict(3 => 4.0), vertices_head=Dict(4 => 5.5))
    SimpleHypergraphs.add_hyperedge!(dh2; vertices_tail=Dict(5 => 7.0), vertices_head=Dict(6 => 8.5))
    SimpleHypergraphs.add_hyperedge!(dh2; vertices_tail=Dict(6 => 10.0))
    SimpleHypergraphs.add_hyperedge!(dh2; vertices_head=Dict(5 => 1.5))
    SimpleHypergraphs.add_vertex!(dh2; hyperedges_tail=Dict(6 => 0.0), hyperedges_head=Dict(5 => -1.5))
    @test dh1 == dh2
    mtail = Matrix(dh1.hg_tail)
    mhead = Matrix(dh1.hg_head)
    @test mtail == Matrix(dh2.hg_tail)
    @test mhead == Matrix(dh2.hg_head)
    @test dh1 == DirectedHypergraph(mtail, mhead)
    @test dh1 == DirectedHypergraph{Float64}(mtail, mhead)
    @test dh1 == DirectedHypergraph{Float64,Nothing}(mtail, mhead)
    @test dh1 == DirectedHypergraph{Float64,Nothing,Nothing}(mtail, mhead)
    @test dh1 == DirectedHypergraph{Float64,Nothing,Nothing,Dict{Int,Float64}}(mtail, mhead)
    @test all(Matrix(dh1.hg_tail) .== Matrix(
        DirectedHypergraph{Float64,Nothing,Nothing,SortedDict{Int,Float64}}(mtail, mhead).hg_tail)
    )
    @test all(Matrix(dh1.hg_head) .== Matrix(
        DirectedHypergraph{Float64,Nothing,Nothing,SortedDict{Int,Float64}}(mtail, mhead).hg_head)
    )
    @test getindex(dh1, 5, 4) == (7.0, nothing)

    dh4 = DirectedHypergraph{Float64,String,Nothing}(1, 1)
    @test SimpleHypergraphs.add_vertex!(dh4; v_meta="test") == 2
    @test SimpleHypergraphs.set_vertex_meta!(dh4, "t", 1) == ["t", "test"]
    @test SimpleHypergraphs.get_vertex_meta(dh4, 2) == "test"
    @test get_hyperedge_meta(dh4, 1) == (nothing, nothing)
    @test_throws BoundsError get_hyperedge_meta(dh4, 2)

    dh5 = DirectedHypergraph{Float64,Nothing,String}(1, 1)
    @test SimpleHypergraphs.add_hyperedge!(dh5; he_meta_tail="test") == 2
    @test SimpleHypergraphs.set_hyperedge_meta!(dh5, "t", "h", 1) == (["t", "test"], ["h", nothing])
    @test get_hyperedge_meta(dh5, 2) == ("test", nothing)
    @test get_vertex_meta(dh5, 1) === nothing
    @test_throws BoundsError get_vertex_meta(dh5, 2)

    dh6 = DirectedHypergraph{Float64,String,String,SortedDict{Int,Float64}}(1, 1)
    @test typeof(dh6.hg_tail.v2he[1]) <: AbstractDict{Int,Float64}
    @test typeof(dh6.hg_head.v2he[1]) <: AbstractDict{Int,Float64}
    @test typeof(dh6.hg_tail.he2v[1]) <: AbstractDict{Int,Float64}
    @test typeof(dh6.hg_head.he2v[1]) <: AbstractDict{Int,Float64}
    @test SimpleHypergraphs.add_vertex!(dh6; v_meta="test") == 2
    @test SimpleHypergraphs.set_vertex_meta!(dh6, "t", 1) == ["t", "test"]
    @test get_vertex_meta(dh6, 2) == "test"
    @test get_hyperedge_meta(dh6, 1) == (nothing, nothing)
    @test SimpleHypergraphs.add_hyperedge!(dh6; he_meta_tail="test") == 2
    @test SimpleHypergraphs.set_hyperedge_meta!(dh6, "t", "h", 1) == (["t", "test"], ["h", nothing])
    @test get_hyperedge_meta(dh6, 2) == ("test", nothing)
    @test_throws BoundsError get_vertex_meta(dh6, 3)
    @test_throws BoundsError get_hyperedge_meta(dh6, 3)
    dh6.hg_tail .= [1.0 2.0; 3.0 4.0]
    @test dh6[2, 2] == (4.0, nothing)

    dh1_0 = deepcopy(dh1)
    @test SimpleHypergraphs.add_vertex!(dh1_0) == 8
    dh1_0.hg_tail[8, :] = dh1_0.hg_tail[7, :]
    dh1_0.hg_head[8, :] = dh1_0.hg_head[7, :]
    @test SimpleHypergraphs.remove_vertex!(dh1_0, 8) == dh1
    setindex!(dh1_0, nothing, 1, 1)
    @test dh1_0[1, 1] == (nothing, nothing)
    @test_throws BoundsError setindex!(dh1_0, nothing, 10, 9)

    dh1_1 = DirectedHypergraph(
        [
            1 nothing nothing
            nothing nothing 2
            nothing nothing nothing
        ],
        [
            nothing nothing 2
            2 nothing 1
            nothing nothing nothing
        ]
    )
    @test SimpleHypergraphs.add_hyperedge!(dh1_1) == 4
    @test size(SimpleHypergraphs.remove_hyperedge!(dh1_1, 4))[2] == 3
    @test SimpleHypergraphs.add_vertex!(dh1_1) == 4
    @test SimpleHypergraphs.add_hyperedge!(dh1_1) == 4
    hp = SimpleHypergraphs.prune_hypergraph(dh1_1)
    @test size(hp)[1] == 2 && size(dh1_1)[1] == 4
    @test size(hp)[2] == 2 && size(dh1_1)[2] == 4
    SimpleHypergraphs.prune_hypergraph!(dh1_1)
    @test size(dh1_1)[1] == 2
    @test size(dh1_1)[2] == 2
end;

# TODO: you are here
@testset "SimpleDirectedHypergraphs BipartiteView          " begin
    dh2 = deepcopy(dh1)

    @test Graphs.nv(Graphs.zero(BipartiteView{DirectedHypergraph{Int}})) == 0

    b = BipartiteView(dh2)
    @test Graphs.edgetype(b) == Graphs.SimpleGraphs.SimpleEdge{Int}
    @test Graphs.has_vertex(b, 0) == false
    @test Graphs.has_vertex(b, 1) == true
    @test Graphs.has_edge(b, 1, 1) == false
    @test Graphs.nv(Graphs.zero(b)) == 0

    @test Graphs.is_directed(b) == true
    @test Graphs.is_directed(typeof(b)) == true
    @test Graphs.eltype(b) == Int


    @test sum(Graphs.adjacency_matrix(Graphs.SimpleDiGraph(b))) == 11

    @test sort(collect(Graphs.outneighbors(b, 5))) == [11]
    @test sort(collect(Graphs.outneighbors(b, 1))) == [8]
    @test sort(collect(Graphs.inneighbors(b, 2))) == [8]

    @test Set(Graphs.vertices(b)) == Set(1:Graphs.nv(b))

    @test SimpleHypergraphs.shortest_path(b, 1, 4) == [1, 3, 4]
    @test SimpleHypergraphs.shortest_path(b, 1, 5) == Int64[]
    @test Graphs.is_weakly_connected(b) == false

    @test SimpleHypergraphs.add_vertex!(dh2) == 8
    @test SimpleHypergraphs.add_hyperedge!(dh2) == 7
    dh2[1, 4, 7] = 1
    dh2[2, 8, 7] = 1

    @test SimpleHypergraphs.shortest_path(b, 1, 8) == [1, 3, 4, 8]

    bipartite_graph = Graphs.SimpleDiGraph(b)

    @test Graphs.SimpleGraphs.fadj(bipartite_graph) == Graphs.SimpleGraphs.fadj(b)
    @test Graphs.nv(b) == 15
    @test Graphs.ne(b) == 13

    @test sort!(Graphs.SimpleGraphs.fadj(b, 1)) == [9]
    @test sort!(Graphs.SimpleGraphs.fadj(b, 2)) == Int64[]
    @test sort!(Graphs.SimpleGraphs.badj(b, 2)) == [9]
end;


@testset "SimpleDirectedHypergraphs TwoSectionView         " begin
    ht = DirectedHypergraph{Float64}(3, 3)
    ht[1, 1, 1] = 1
    ht.hg_head[2:3, 1] .= 2
    ht[1, 2, 2] = 2
    ht[2, 1, 2] = 2
    ht[1, 3, 3] = 3
    ht[2, 1, 3] = 3

    @test Graphs.nv(Graphs.zero(TwoSectionView{DirectedHypergraph{Int64}})) == 0

    t = TwoSectionView(dh1)
    @test Graphs.edgetype(t) == Graphs.SimpleGraphs.SimpleEdge{Int}
    @test Graphs.has_vertex(t, 0) == false
    @test Graphs.has_vertex(t, 1) == true
    @test Graphs.nv(Graphs.zero(t)) == 0

    @test Graphs.is_directed(t) == true
    @test Graphs.is_directed(typeof(t)) == true
    @test Graphs.eltype(t) == Int

    @test Graphs.nv(t) == 7
    @test Graphs.ne(t) == 6

    @test sort(Graphs.all_neighbors(t, 1)) == [2, 3]
    @test sort(Graphs.outneighbors(t, 5)) == [6]
    @test sort(Graphs.inneighbors(t, 4)) == [3]
    @inferred Graphs.all_neighbors(t, 1)

    @test Graphs.has_edge(t, 1, 2) == true
    @test Graphs.has_edge(t, 1, 5) == false

    @test sum(Graphs.adjacency_matrix(Graphs.SimpleDiGraph(t))) == 6
    @test SimpleHypergraphs.shortest_path(t, 1, 4) == [1, 3, 4]
    @test Graphs.is_weakly_connected(t) == false
    @test Graphs.is_strongly_connected(t) == false

    SimpleHypergraphs.add_vertex!(dh1)
    SimpleHypergraphs.add_hyperedge!(dh1)
    dh1[1, 4, 7] = 1
    dh1[2, 8, 7] = 1

    @test Graphs.ne(t) == 7
    @test Graphs.nv(t) == 8
    @test sort(Graphs.outneighbors(t, 4)) == [8]

    @test SimpleHypergraphs.shortest_path(t, 1, 8) == [1, 3, 4, 8]

    @test sum(Graphs.adjacency_matrix(Graphs.SimpleDiGraph(t))) == 7

    Random.seed!(0)
    g = Graphs.erdos_renyi(8, 0.3; is_directed=true)
    h_from_g = DirectedHypergraph(g)
    @test Graphs.adjacency_matrix(g) == Graphs.adjacency_matrix(TwoSectionView(h_from_g))
    @test minimum([sum((h_from_g.hg_tail.==true)[:, n]) for n in 1:6] .== 1)
    @test minimum([sum((h_from_g.hg_head.==true)[:, n]) for n in 1:6] .== 1)
    @test Graphs.SimpleGraphs.fadj(g) == Graphs.SimpleGraphs.fadj(TwoSectionView(h_from_g))
    @test Graphs.SimpleGraphs.badj(g) == Graphs.SimpleGraphs.badj(TwoSectionView(h_from_g))
end;


@testset "SimpleDirectedHypergraphs random-models          " begin
    DHᵣ = random_model(5, 5, DirectedHypergraph)
    @test nhv(DHᵣ) == 5
    @test nhe(DHᵣ) == 5
    @test all(length.(DHᵣ.hg_tail.v2he) .> 0)
    @test all(length.(DHᵣ.hg_head.v2he) .> 0)
    @test all(length.(DHᵣ.hg_tail.v2he) .<= 5)
    @test all(length.(DHᵣ.hg_head.v2he) .<= 5)

    @test_throws ErrorException random_model(1, 2, DirectedHypergraph; no_self_loops=true)

    DHr_nsl = random_model(5, 5, DirectedHypergraph; no_self_loops=true)
    for i in 1:5
        @test length(intersect(keys(DHr_nsl.hg_tail.v2he[i]), keys(DHr_nsl.hg_head.v2he[i]))) == 0
    end

    DHᵣ2 = random_model(5, 0, DirectedHypergraph)
    add_hyperedge!(DHᵣ2; vertices_tail=Dict(2 => true, 4 => true), vertices_head=Dict(1 => true, 5 => true))
    @test nhv(DHᵣ2) == 5
    @test nhe(DHᵣ2) == 1

    DHκ = random_kuniform_model(5, 5, 3, DirectedHypergraph)
    @test nhv(DHκ) == 5
    @test nhe(DHκ) == 5
    @test all(length.(DHκ.hg_tail.he2v) .+ length.(DHκ.hg_head.he2v) .== 3)

    @test_throws ErrorException random_kuniform_model(1, 3, 1, DirectedHypergraph; no_self_loops=true)

    DHκ_nsl = random_kuniform_model(5, 5, 3, DirectedHypergraph; no_self_loops=true)
    for i in 1:5
        @test length(intersect(keys(DHκ_nsl.hg_tail.v2he[i]), keys(DHκ_nsl.hg_head.v2he[i]))) == 0
    end

    DHδ = random_dregular_model(5, 5, 3, DirectedHypergraph)
    @test nhv(DHδ) == 5
    @test nhe(DHδ) == 5
    @test all(length.(DHδ.hg_tail.v2he) .+ length.(DHδ.hg_head.v2he) .== 3)

    @test_throws ErrorException random_dregular_model(1, 3, 1, DirectedHypergraph; no_self_loops=true)

    DHδ_nsl = random_kuniform_model(5, 5, 3, DirectedHypergraph; no_self_loops=true)
    for i in 1:5
        @test length(intersect(keys(DHκ_nsl.hg_tail.v2he[i]), keys(DHκ_nsl.hg_head.v2he[i]))) == 0
    end
end;

@testset "SimpleDirectedHypergraphs randomwalk             " begin
    dh = DirectedHypergraph{Float64}(8, 9)
    dh[1, 1, 1] = 1.0
    dh.hg_head[2:3, 1] .= 2.5
    dh[1, 3, 2] = 0.5
    dh[2, 4, 2] = 1.5
    dh.hg_tail[3:4, 3] .= 2.0
    dh[2, 1, 3] = 0.5
    dh[1, 4, 4] = 2.5
    dh[2, 5, 4] = 1.0
    dh[1, 3, 5] = 2.0
    dh[2, 5, 5] = 1.0
    dh[1, 5, 6] = 1.5
    dh[2, 3, 6] = 0.5
    dh[1, 6, 7] = 1.0
    dh[2, 7, 7] = 2.5
    dh[1, 7, 8] = 2.0
    dh[2, 8, 8] = 1.5
    dh[1, 8, 9] = 1.0
    dh[2, 6, 9] = 1.5

    dw1 = countmap([random_walk(dh, 1) for _ in 1:10^6])
    @test keys(dw1) == Set([2, 3])
    @test -(extrema(values(dw1))...) > -10000

    dw1_rev = countmap([random_walk(dh, 1; reverse=true) for _ in 1:10^6])
    @test keys(dw1_rev) == Set([3, 4])
    @test -(extrema(values(dw1_rev))...) > -10000

    dw2 = countmap([random_walk(dh, 6) for _ in 1:10^6])
    @test keys(dw2) == Set([7])
    @test -(extrema(values(dw2))...) > -10000

    @test_throws ArgumentError random_walk(dh, 0)

end;

@testset "SimpleDirectedHypergraphs connected components" begin
    weak_conn = sort!(get_weakly_connected_components(dh1))
    @test length(weak_conn) == 2
    @test weak_conn[1] == [1, 2, 3, 4, 8]
    @test weak_conn[2] == [5, 6, 7]

    strong_conn = sort!(get_strongly_connected_components(dh1))
    @test length(strong_conn) == 6
    @test strong_conn[1] == [1]
    @test strong_conn[2] == [2]
    @test strong_conn[3] == [3]
    @test strong_conn[4] == [4]
    @test strong_conn[5] == [5, 6, 7]
    @test strong_conn[6] == [8]
end;


@testset "SimpleDirectedHypergraphs dual                   " begin
    m_tail = [
        1 nothing nothing nothing
        nothing 2 nothing nothing
        nothing nothing 3 3
        4 nothing nothing 4
    ]

    m_head = [
        nothing 5 5 nothing
        nothing nothing nothing 6
        7 nothing nothing nothing
        nothing nothing 8 nothing
    ]

    v_meta = Array{Union{Nothing,Char},1}(collect('a':'d'))
    he_meta_tail = Array{Union{Nothing,Symbol},1}(Symbol.(collect('A':'D')))
    he_meta_head = Array{Union{Nothing,Symbol},1}(Symbol.(collect('E':'H')))

    dh = DirectedHypergraph{Int,Char,Symbol}(m_tail, m_head; v_meta=v_meta, he_meta_tail=he_meta_tail, he_meta_head=he_meta_head)
    dh_dual = dual(dh)

    @test nhv(dh_dual) == nhe(dh)
    @test nhe(dh_dual) == nhv(dh)

    m_dual_tail = Matrix(dh_dual.hg_tail)
    m_dual_head = Matrix(dh_dual.hg_head)

    m_dual_tail[m_dual_tail.==nothing] .= 0
    m_dual_head[m_dual_head.==nothing] .= 0
    m_tail[m_tail.==nothing] .= 0
    m_head[m_head.==nothing] .= 0

    @test m_tail == transpose(m_dual_tail)
    @test m_head == transpose(m_dual_head)

    @test dh_dual.he_meta_tail == dh.v_meta
    @test dh_dual.he_meta_head == dh.v_meta
    @test all(dh_dual.v_meta .== [(dh.he_meta_tail[i], dh.he_meta_head[i]) for i in 1:nhe(dh)])

    @test_throws AssertionError dual(DirectedHypergraph(0, 0))
end;


@testset "SimpleDirectedHypergraphs paths/distance         " begin
    # Custom, strongly-connected example with loops
    w2 = [3, 7, 19, 10, 13, 11, 1, 12, 8]

    # Traversing forward from vertices to vertices/hyperedges
    fr = forward_reachable(dh2, 1)
    @test fr[1] == Set{Int}(1:5)
    @test fr[2] == Set{Int}(1:9)
    fr = forward_reachable(dh2, 5)
    @test fr[1] == Set{Int}(1:5)
    @test fr[2] == Set{Int}(1:9)

    # Back-tracing from vertices to vertices/hyperedges
    bt = backward_traceable(dh2, 1)
    @test bt[1] == Set{Int}(1:5)
    @test bt[2] == Set{Int}(1:9)
    bt = backward_traceable(dh2, 5)
    @test bt[1] == Set{Int}(1:5)
    @test bt[2] == Set{Int}(1:9)

    for i in 1:nhv(dh2)
        # In this example, all nodes (and hyperedges) can be reached from anywhere
        for j in 1:nhv(dh2)
            @test is_reachable(dh2, i, j, :vertex)
        end
        for e in 1:nhe(dh2)
            @test is_reachable(dh2, i, e, :hyperedge)
        end
    end
    @test_throws AssertionError is_reachable(dh2, 1, 2, :test)

    # Test `short_hyperpath_vhe` for greedy hyperpath generation
    he_inedges = [
        Set{Int}([6,7]),
        Set{Int}([1,8]),
        Set{Int}([1,9]),
        Set{Int}([2,3,5]),
        Set{Int}(4),
        Set{Int}([1,2,3,5,8]),
        Set{Int}([1,2,3,5,9]),
        Set{Int}(4),
        Set{Int}(4)
    ]

    # Test method to find a reasonably short path from a vertex to a hyperedge
    short_path = SimpleDirectedHypergraphs.short_hyperpath_vhe(dh2, 1, 4, he_inedges, w2)
    @test short_path == Set{Int}([1,2,4]) || short_path == Set{Int}([1,3,4])

    # Test heuristic (but usually accurate) shortest-path algorithm
    @test shortest_hyperpath_kk_heuristic(dh2, 1, 5, w2) == Set{Int}([1,2,4])
    @test shortest_hyperpath_kk_heuristic(dh2, 1, Set{Int}([3,4]), w2) == Set{Int}([1,2])
    @test shortest_hyperpath_kk_heuristic(dh2, Set{Int}([4,5]), 1, w2) == Set{Int}([7,9])
    @test shortest_hyperpath_kk_heuristic(dh2, Set{Int}([2, 3]), Set{Int}([4, 1]), w2) == Set{Int}([2,7])

    # Test generation of all possible paths
    @test length(all_hyperpaths(dh2, 1, 5)) == 2

    # Example adapted from Blau et al., DOI: 10.1039/D0SC05647B
    tail_3 = [
        true true true true nothing nothing nothing nothing nothing nothing nothing
        nothing nothing nothing nothing true nothing nothing nothing nothing nothing nothing
        nothing nothing nothing nothing nothing true nothing nothing nothing nothing nothing
        nothing nothing nothing nothing nothing nothing true nothing nothing nothing nothing
        nothing nothing nothing nothing nothing nothing nothing true nothing nothing nothing
        nothing nothing nothing nothing nothing nothing nothing nothing true nothing nothing
        nothing nothing nothing nothing nothing nothing nothing nothing nothing true nothing
        nothing nothing nothing nothing nothing nothing nothing nothing nothing nothing true
        nothing nothing nothing nothing nothing nothing nothing nothing nothing nothing nothing
    ]
    head_3 = [
        nothing nothing nothing nothing nothing nothing nothing nothing nothing nothing nothing
        nothing true nothing nothing nothing nothing nothing nothing nothing nothing nothing
        true nothing nothing nothing nothing nothing nothing nothing nothing nothing nothing
        true nothing nothing nothing nothing nothing nothing nothing nothing nothing nothing
        nothing nothing true nothing nothing nothing nothing nothing nothing nothing nothing
        nothing nothing true nothing nothing nothing nothing nothing nothing nothing nothing
        nothing nothing nothing true nothing nothing true true nothing nothing nothing
        nothing nothing nothing nothing nothing true nothing nothing true nothing nothing
        nothing nothing nothing nothing true nothing nothing nothing nothing true true
    ]
    dh3 = DirectedHypergraph(tail_3, head_3)

    w3 = [1, 6, 2, 4, 3, 1, 2, 2, 4, 1, 1]

    # Traversing forward from vertices to vertices/hyperedges
    fr = forward_reachable(dh3, 1)
    @test fr[1] == Set{Int}(1:9)
    @test fr[2] == Set{Int}(1:11)
    fr = forward_reachable(dh3, 6)
    @test fr[1] == Set{Int}([6, 8, 9])
    @test fr[2] == Set{Int}([9, 11])
    fr = forward_reachable(dh3, 9)
    @test fr[1] == Set{Int}(9)
    @test fr[2] == Set{Int}()

    # Back-tracing from vertices to vertices/hyperedges
    bt = backward_traceable(dh3, 1)
    @test bt[1] == Set{Int}(1)
    @test bt[2] == Set{Int}()
    bt = backward_traceable(dh3, 6)
    @test bt[1] == Set{Int}([6, 1])
    @test bt[2] == Set{Int}(3)
    bt = backward_traceable(dh3, 9)
    @test bt[1] == Set{Int}(1:9)
    @test bt[2] == Set{Int}(1:11)

    # Test heuristic (but usually accurate) shortest-path algorithm
    @test shortest_hyperpath_kk_heuristic(dh3, 1, 9, w3) == Set{Int}([1, 6, 11])

    # Test generation of all possible paths
    @test length(all_hyperpaths(dh3, 1, 9)) == 6
    @test length(all_hyperpaths(dh3, 1, Set{Int}([8, 9]))) == 9
    @test length(all_hyperpaths(dh3, Set{Int}([1,8]), 9)) == 1
    @test length(all_hyperpaths(dh3, Set{Int}([2,6]), Set{Int}([8,9]))) == 2
end;

@testset "SimpleDirectedHypergraphs diameter               " begin

end;