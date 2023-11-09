using Scaling
using Documenter

DocMeta.setdocmeta!(Scaling, :DocTestSetup, :(using Scaling); recursive=true)

makedocs(;
    modules=[Scaling],
    authors="Malte Pütz",
    repo="https://github.com/maltepuetz/Scaling.jl/blob/{commit}{path}#{line}",
    sitename="Scaling.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://maltepuetz.github.io/Scaling.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/maltepuetz/Scaling.jl",
    devbranch="main",
)
