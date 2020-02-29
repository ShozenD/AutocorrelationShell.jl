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
    n = size(x)[1]
    y = deepcopy(x) # to prevent inplace but slows speed. 
    for i in 1:n
        for j in 1:n
            if type=="hard"
                HardThreshold!(y[i][j], t)
            elseif type=="soft"
                SoftThreshold!(y[i][j], t)
            end
        end
    end
    return y
end

end # module
