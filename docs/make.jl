using ReferenceFiniteElements
using Documenter

DocMeta.setdocmeta!(ReferenceFiniteElements, :DocTestSetup, :(using ReferenceFiniteElements); recursive=true)

makedocs(;
    modules=[ReferenceFiniteElements],
    authors="Craig M. Hamel <cmhamel32@gmail.com> and contributors",
    repo="https://github.com/Cthonios/ReferenceFiniteElements.jl/blob/{commit}{path}#{line}",
    sitename="ReferenceFiniteElements.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://Cthonios.github.io/ReferenceFiniteElements.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/Cthonios/ReferenceFiniteElements.jl",
    devbranch="main",
)
