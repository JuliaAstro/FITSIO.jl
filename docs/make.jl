using Documenter, DocumenterInterLinks, FITSIO
using DataFrames, WCS # Precompile package extensions

links = InterLinks(
    "WCS" => "https://juliaastro.org/WCS/stable/",
)

include("pages.jl")
makedocs(;
    modules = [FITSIO, Base.get_extension(FITSIO, :WCSExt)],
    sitename = "FITSIO.jl",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://juliaastro.org/FITSIO/stable/",
    ),
    pages = pages,
    checkdocs = :exports,
    plugins = [links],
)

deploydocs(;
    repo = "github.com/JuliaAstro/FITSIO.jl.git",
    push_preview = true,
    versions = ["stable" => "v^", "v#.#"], # Restrict to minor releases
)
