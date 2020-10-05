push!(LOAD_PATH, "../src/")

using
    Documenter,
    AutocorrelationShell

makedocs(
    sitename="AutocorrelationShell.jl",
    modules = [AutocorrelationShell],
    authors = "Naoki Saito, Rishi Subramanian, Christina Chang, and Shozen Dan",

    doctest = true,

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
    repo = "https://github.com/ShozenD/AutocorrelationShell.git"
)
