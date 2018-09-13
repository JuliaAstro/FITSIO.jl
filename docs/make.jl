using Documenter, FITSIO

makedocs(
    modules = [FITSIO],
    format = :html,
    sitename = "FITSIO.jl",
    pages    = Any[
        "Introduction" => "index.md",
        "API Reference" => "api.md",
        "Libcfitsio Submodule" => "libcfitsio.md"
    ]
)

deploydocs(
    repo = "github.com/JuliaAstro/FITSIO.jl.git",
    target = "build",
    deps = nothing,
    make = nothing,
    julia  = "1.0",
)
