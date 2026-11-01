using ScalingCollapse
using Documenter
using DocumenterVitepress

DocMeta.setdocmeta!(ScalingCollapse, :DocTestSetup, :(using ScalingCollapse); recursive=true)

makedocs(;
    modules=[ScalingCollapse],
    authors="Malte Pütz",
    sitename="ScalingCollapse.jl",
    checkdocs=:exports,
    format=DocumenterVitepress.MarkdownVitepress(;
        repo="github.com/maltepuetz/ScalingCollapse.jl",
        devbranch="main",
        devurl="dev",
        deploy_url="https://maltepuetz.github.io/ScalingCollapse.jl",
    ),
    pages=[
        "Home" => "index.md",
        "Manual" => [
            "Getting Started" => "getting-started.md",
        ],
        "API Reference" => "api.md",
    ],
)

DocumenterVitepress.deploydocs(;
    repo="github.com/maltepuetz/ScalingCollapse.jl",
    target=joinpath(@__DIR__, "build"),
    branch="gh-pages",
    devbranch="main",
    push_preview=true,
)
