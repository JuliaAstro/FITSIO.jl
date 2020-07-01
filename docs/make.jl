using Documenter, FITSIO

makedocs(
    modules = [FITSIO],
    sitename = "FITSIO.jl",
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    pages    = [
        "Introduction" => "index.md",
        "API Reference" => "api.md"
    ]
)

deploydocs(repo = "github.com/JuliaAstro/FITSIO.jl.git", push_preview=true)
