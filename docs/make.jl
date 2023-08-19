using AiidaDFTK
using Documenter

DocMeta.setdocmeta!(AiidaDFTK, :DocTestSetup, :(using AiidaDFTK); recursive=true)

makedocs(;
    modules=[AiidaDFTK],
    authors="Michael F. Herbst <info@michael-herbst.com> and contributors",
    repo="https://github.com/mfherbst/AiidaDFTK.jl/blob/{commit}{path}#{line}",
    sitename="AiidaDFTK.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://mfherbst.github.io/AiidaDFTK.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/mfherbst/AiidaDFTK.jl",
    devbranch="master",
)
