using FFTGrp
using Documenter

DocMeta.setdocmeta!(FFTGrp, :DocTestSetup, :(using FFTGrp); recursive=true)

makedocs(;
    modules=[FFTGrp],
    authors="Naoki Maeda2 <fnkyksj@gmail.com> and contributors",
    repo="https://github.com/nmaedajp/FFTGrp.jl/blob/{commit}{path}#{line}",
    sitename="FFTGrp.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://nmaedajp.github.io/FFTGrp.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/nmaedajp/FFTGrp.jl",
    devbranch="main",
)
