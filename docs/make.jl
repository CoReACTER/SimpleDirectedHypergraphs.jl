using SimpleDirectedHypergraphs
using Documenter

DocMeta.setdocmeta!(SimpleDirectedHypergraphs, :DocTestSetup, :(using SimpleDirectedHypergraphs); recursive=true)

makedocs(;
    modules=[SimpleDirectedHypergraphs],
    authors="Evan Walter Clark Spotte-Smith",
    sitename="SimpleDirectedHypergraphs.jl",
    format=Documenter.HTML(;
        canonical="https://coreacter.codeberg.page/SimpleDirectedHypergraphs.jl/@main",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
    checkdocs=:exports
)

deploydocs(;
    repo="codeberg.org/CoReACTER/SimpleDirectedHypergraphs.jl",
    devbranch="main",
)
