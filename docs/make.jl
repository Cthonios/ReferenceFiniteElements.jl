using Adapt
using Exodus
using LaTeXStrings
using RecipesBase
using ReferenceFiniteElements
using Symbolics
using Documenter

DocMeta.setdocmeta!(ReferenceFiniteElements, :DocTestSetup, :(using ReferenceFiniteElements); recursive=true)

makedocs(;
    modules=[ReferenceFiniteElements],
    authors="Craig M. Hamel <cmhamel32@gmail.com> and contributors",
    repo="https://github.com/Cthonios/ReferenceFiniteElements.jl/blob/{commit}{path}#{line}",
    sitename="ReferenceFiniteElements.jl",
    format=Documenter.HTML(;
        repolink="https://github.com/Cthonios/ReferenceFiniteElements.jl",
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://cthonios.github.io/ReferenceFiniteElements.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home"       => [
            "Installation"  => "installation.md",
            # "Quick start"   => "quick_start.md",
            # "Storage types" => "storage_types.md"
        ],
        # "Extensions" => "extensions.md",
        "Abstract Types" => "abstract_types.md",
        "Element Types" => "element_types.md",
        "Developer"  => "developer.md",
        "Index"      => "index.md"
    ],
)

deploydocs(;
    repo="github.com/Cthonios/ReferenceFiniteElements.jl",
    devbranch="main",
)
