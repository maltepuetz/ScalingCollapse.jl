using ScalingCollapse
using Documenter
using DocumenterVitepress

include("figure_theme.jl")

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
        "Getting Started" => "getting-started.md",
        "Manual" => [
            "Finite-size scaling" => "manual/background.md",
            "Input data" => "manual/data.md",
            "Scaling functions" => "manual/scaling-functions.md",
            "Fitting" => "manual/fitting.md",
            "Error analysis" => "manual/errors.md",
            "Quality functions" => "manual/quality-functions.md",
        ],
        "Examples" => [
            "2D Ising model" => "examples/ising.md",
            "Site percolation with error bars" => "examples/percolation.md",
            "Specific heat: a custom scaling function" => "examples/specific-heat.md",
            "Mean-field model: exact data" => "examples/curie-weiss.md",
            "Uniform grids and quality landscapes" => "examples/landscapes.md",
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
