using Documenter, FITSIO

makedocs(
    modules = [FITSIO],
    format = :html,
    sitename = "FITSIO",
    pages    = Any[
        "Introduction"   => "index.md",
    ]
)

deploydocs(
    repo = "github.com/JuliaAstro/FITSIO.jl.git",
    target = "build",
)
