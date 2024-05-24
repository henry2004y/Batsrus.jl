using Batsrus, Documenter
using DemoCards

branch = "master"
# generate demo files
demos, postprocess_cb, demo_assets = makedemos("examples"; branch)
# if there are generated css assets, pass it to Documenter.HTML
assets = String[]
isnothing(demo_assets) || (push!(assets, demo_assets))

makedocs(;
    modules=[Batsrus, Batsrus.UnitfulBatsrus, Batsrus.HDF],
    authors="Hongyang Zhou <hyzhou@umich.edu>",
    sitename="Batsrus.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://henry2004y.github.io/Batsrus.jl",
        assets,
        sidebar_sitename=false
    ),
    pages=[
        "Home" => "index.md",
        "Manual" => "man/manual.md",
        "Examples" => demos,
        "Internal" => "man/internal.md",
        "Log" => "man/log.md",
    ],
)

deploydocs(;
    repo="github.com/henry2004y/Batsrus.jl",
)
