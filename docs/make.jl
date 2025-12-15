using SimpleDirectedHypergraphs
using Documenter
using OpenSSL_jll


DocMeta.setdocmeta!(SimpleDirectedHypergraphs, :DocTestSetup, :(using SimpleDirectedHypergraphs); recursive=true)

doctest(SimpleDirectedHypergraphs)

makedocs(;
    modules=[SimpleDirectedHypergraphs],
    authors="Evan Walter Clark Spotte-Smith, CoReACTER",
    sitename="SimpleDirectedHypergraphs.jl",
    format=Documenter.HTML(;
        canonical="coreacter.codeberg.page/SimpleDirectedHypergraphs.jl/",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
    checkdocs=:exports,
    repo="codeberg.org/CoReACTER/SimpleDirectedHypergraphs.jl"
)

deploydocs(;
    repo="codeberg.org/CoReACTER/SimpleDirectedHypergraphs.jl",
    branch="pages",
    devbranch="main",
)
