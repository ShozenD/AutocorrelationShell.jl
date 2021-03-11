using
    Documenter,
    AutocorrelationShell

DocMeta.setdocmeta!(AutocorrelationShell, :DocTestSetup, :(using AutocorrelationShell); recursive=true)

makedocs(
    sitename="AutocorrelationShell.jl",
    modules = [AutocorrelationShell],
    authors = "Shozen Dan, Rishi Subramanian, Christina Chang, and Naoki Saito",
    repo="https://github.com/ShozenD/AutocorrelationShell.jl/blob/{commit}{path}#{line}",

    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://ShozenD.github.io/AutocorrelationShell.jl",
        assets=String[],
    ),

    pages = [
        "Home" => "index.md",
        "Getting Started" => "getting_started.md",
        "AC Wavelet 1D" => "acw1d.md",
        "AC Wavelet 2D" => "acw2d.md",
        "AC Wavelet Packets" => "acwpt.md",
        "AC Wavelet Utils" => "acwutil.md"
    ]
)

deploydocs(
    repo = "https://github.com/ShozenD/AutocorrelationShell.jl"
)

println("\n", base64encode(read(filename, String)), "\n")
