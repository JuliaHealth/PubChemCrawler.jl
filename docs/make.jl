using PubChemCrawler
using Documenter

makedocs(;
    modules=[PubChemCrawler],
    authors="Tim Holy <tim.holy@gmail.com> and contributors",
    sitename="PubChemCrawler.jl",
    format=Documenter.HTML(;
        canonical="https://JuliaHealth.github.io/PubChemCrawler.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/JuliaHealth/PubChemCrawler.jl",
    push_preview=true,
)
