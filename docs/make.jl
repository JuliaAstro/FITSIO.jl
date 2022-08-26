using Documenter, FITSIO

include("pages.jl")
makedocs(
    modules = [FITSIO],
    sitename = "FITSIO.jl",
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    pages = pages
)

deploydocs(repo = "github.com/JuliaAstro/FITSIO.jl.git", push_preview=true)
