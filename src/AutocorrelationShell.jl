__precompile__()

module AutocorrelationShell
@warn("The AutocorrelationShell.jl is deprecated and will be part of WaveletExts.jl from April 23, 2021:
        pkg> update
        pkg> add WaveletsExt")

using SpecialFunctions, Reexport

include("mod/ACWT.jl")
include("mod/ACTransforms.jl")
include("mod/ACPlots.jl")

@reexport using .ACWT, .ACTransforms, .ACPlots

end # module
