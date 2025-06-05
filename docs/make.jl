using Documenter, FITSIO

include("pages.jl")
makedocs(
    modules = [FITSIO],
    sitename = "FITSIO.jl",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://juliaastro.org/FITSIO/stable",
    ),
    pages = pages,
    checkdocs = :exports,
)

deploydocs(
    repo = "github.com/JuliaAstro/FITSIO.jl.git",
    push_preview = true,
    versions = ["stable" => "v^", "v#.#"], # Restrict to minor releases
)
