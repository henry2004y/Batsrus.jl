using Batsrus, Documenter, PyPlot
using Batsrus.UnitfulBatsrus

makedocs(;
    modules=[Batsrus],
    authors="Hongyang Zhou <hyzhou@umich.edu>",
    repo="https://github.com/henry2004y/Batsrus.jl/blob/{commit}{path}#L{line}",
    sitename="Batsrus.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://henry2004y.github.io/Batsrus.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Example" => "man/examples.md",
        "Internal" => "man/internal.md",
        "Log" => "man/log.md",
    ],
)

deploydocs(;
    repo="github.com/henry2004y/Batsrus.jl",
)
