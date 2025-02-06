using SimpleDirectedHypergraphs
using Documenter

DocMeta.setdocmeta!(SimpleDirectedHypergraphs, :DocTestSetup, :(using SimpleDirectedHypergraphs); recursive=true)

makedocs(;
    modules=[SimpleDirectedHypergraphs],
    authors="Evan Walter Clark Spotte-Smith",
    sitename="SimpleDirectedHypergraphs.jl",
    format=Documenter.HTML(;
        canonical="https://espottesmith.github.io/SimpleDirectedHypergraphs.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/espottesmith/SimpleDirectedHypergraphs.jl",
    devbranch="main",
)
