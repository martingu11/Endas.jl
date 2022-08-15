
#push!(LOAD_PATH,"../src/")

using Documenter
using Endas


makedocs(
    sitename="Endas.jl",
    format = Documenter.HTML(prettyurls = false),         
    modules = [Endas],
    pages = [
        "index.md",
        "Manual" => [
            "usage.md",
            "observations.md",
            "covariance.md",
            "models.md",
        ],
        "includedmodels.md",
        "examples.md",
        "reference.md"
    ]
)
