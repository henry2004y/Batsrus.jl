using Documenter, VisAna

push!(LOAD_PATH,"../src/")

makedocs(
    sitename="VisAna Documentation",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    )
)

deploydocs(
    repo = "github.com/henry2004y/VisAnaJulia",
    target = "build",
    branch = "gh-pages"
)
