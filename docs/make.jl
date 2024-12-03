using AiidaDFTK
using Documenter

DocMeta.setdocmeta!(AiidaDFTK, :DocTestSetup, :(using AiidaDFTK); recursive=true)

makedocs(;
    modules=[AiidaDFTK],
    authors="Michael F. Herbst <info@michael-herbst.com> and contributors",
    repo="https://github.com/epfl-matmat/AiidaDFTK.jl/blob/{commit}{path}#{line}",
    sitename="AiidaDFTK.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://epfl-matmat.github.io/AiidaDFTK.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "input_output.md",
    ],
)

deploydocs(;
    repo="github.com/epfl-matmat/AiidaDFTK.jl",
    devbranch="master",
)
