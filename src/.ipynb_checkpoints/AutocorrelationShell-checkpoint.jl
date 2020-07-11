__precompile__()

module AutocorrelationShell
using SpecialFunctions

include("mod/AC1D.jl")
include("mod/AC2D.jl")
include("mod/ACThreshold.jl")
include("mod/ACUtil.jl")

using Reexport
@reexport using .AC1D, .AC2D, .ACThreshold, .ACUtil

end # module
