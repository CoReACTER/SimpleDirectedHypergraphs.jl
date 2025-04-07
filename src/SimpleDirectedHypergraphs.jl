module SimpleDirectedHypergraphs

using SimpleHypergraphs
using Graphs
using StatsBase
using DataStructures
# using PyPlot
using JSON3
using Random
using LinearAlgebra
using SimpleTraits

export AbstractDirectedHypergraph
export DirectedHypergraph, HyperedgeDirection

export EHGF_Format
export dhg_load

export to_undirected

export get_weakly_connected_components, get_strongly_connected_components


include("abstracttypes.jl")
include("dihypergraph.jl")
include("io.jl")

include("models/bipartite.jl")
include("models/twosection.jl")
include("models/random-models.jl")
include("models/dual.jl")

end # module
