using Documenter, Example

push!(LOAD_PATH,"../src/")

makedocs(sitename="My Documentation")

deploydocs(
    repo = "github.com/henry2004y/VisAnaJulia.git",
)
