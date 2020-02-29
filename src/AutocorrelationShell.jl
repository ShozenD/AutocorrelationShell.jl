__precompile__()

module AutocorrelationShell
using SpecialFunctions

include("mod/AC1D.jl")
include("mod/AC2D.jl")
include("mod/ACThreshold.jl")

using Reexport
@reexport using .AC1D, .AC2D, .ACThreshold

end # module
