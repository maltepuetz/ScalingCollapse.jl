using ScalingCollapse
using Documenter

DocMeta.setdocmeta!(ScalingCollapse, :DocTestSetup, :(using ScalingCollapse); recursive=true)

makedocs(;
    modules=[ScalingCollapse],
    authors="Malte Pütz",
    repo="https://github.com/maltepuetz/ScalingCollapse.jl/blob/{commit}{path}#{line}",
    sitename="ScalingCollapse.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://maltepuetz.github.io/ScalingCollapse.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/maltepuetz/ScalingCollapse.jl",
    devbranch="main",
)
