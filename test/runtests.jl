using SimpleDirectedHypergraphs
using Test
using JET

@testset "SimpleDirectedHypergraphs.jl" begin
    @testset "Code linting (JET.jl)" begin
        JET.test_package(SimpleDirectedHypergraphs; target_defined_modules = true)
    end
    # Write your tests here.
end
