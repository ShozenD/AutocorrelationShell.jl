__precompile__()

module AutocorrelationShell
using SpecialFunctions

include("mod/AC1D.jl")
include("mod/AC2D.jl")
include("mod/ACThreshold.jl")
include("mod/ACUtil.jl")
include("mod/ACWPT1D.jl")

using Reexport
@reexport using .AC1D, .AC2D, .ACThreshold, .ACUtil, .ACWPT1D

end # module
