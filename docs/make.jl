using Documenter, Hazma

makedocs(;
    modules=[Hazma],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/LoganAMorrison/Hazma.jl/blob/{commit}{path}#L{line}",
    sitename="Hazma.jl",
    authors="Logan A Morrison and Adam Coogan",
    assets=String[],
)

deploydocs(;
    repo="github.com/LoganAMorrison/Hazma.jl",
)
