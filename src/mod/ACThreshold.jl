module ACThreshold
export
    # threshold
    acthreshold

using Wavelets
using LinearAlgebra
using Statistics

# hard
function HardThreshold!(x::AbstractArray{<:Number}, t::Real)
    @assert t>=0
    @inbounds begin
        for i in eachindex(x)
            if abs(x[i]) <= t
                x[i] = 0
            end
        end
    end
    return x
end

# soft
function SoftThreshold!(x::AbstractArray{<:Number}, t::Real)
    @assert t>=0
    @inbounds begin
        for i in eachindex(x)
            sh = abs(x[i]) - t
            if sh < 0
                x[i] = 0
            else
                x[i] = sign(x[i])*sh
            end
        end
    end
    return x
end

# Overall function
"""
    acthreshold(x, type, t)

    Thresholds the ac2d output with either hard of soft thresholding.

    ### Arguments
    `x`: ac2d output.
    `type`: Threshold type. "hard" or "soft".
    `t`: threshold value.
"""
function acthreshold(x, type, t::Real)
    y = deepcopy(x) # to prevent inplace behavior but slows speed
    if type=="hard"
        HardThreshold!(y, t)
    elseif type=="soft"
        SoftThreshold!(y, t)
    end
    return y
end

end # module
