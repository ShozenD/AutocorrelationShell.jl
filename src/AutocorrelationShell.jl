__precompile__()

module AutocorrelationShell
using SpecialFunctions

include("mod/ACWT.jl")
include("mod/ACUtil.jl")
include("mod/ACTransforms.jl")
include("mod/ACPlots.jl")

using Reexport
@reexport using .ACWT, .ACTransforms, .ACUtil, .ACPlots

end # module
