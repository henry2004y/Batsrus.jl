using Documenter, VisAna

#push!(LOAD_PATH,"../src/")

makedocs(
    sitename="VisAna",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    )
)

deploydocs(
    repo = "github.com/henry2004y/VisAnaJulia.git",
    target = "build",
    branch = "gh-pages"
)
