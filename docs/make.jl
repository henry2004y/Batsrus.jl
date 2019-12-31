using Documenter, VisAna

push!(LOAD_PATH,"../src/")

makedocs(sitename="VisAna Documentation")

deploydocs(
    repo = "github.com/henry2004y/VisAnaJulia",
    target = "build",
    branch = "gh-pages"
)
