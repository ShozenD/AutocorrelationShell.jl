module ACThreshold
export
    # threshold
    acthreshold

using Wavelets
using LinearAlgebra
using Statistics

# hard
function HardThreshold(x::AbstractArray{<:Number}, t::Real)
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
function SoftThreshold(x::AbstractArray{<:Number}, t::Real)
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
function acthreshold(ac2d, type, t::Real)
    n = size(ac2d)[1]
    for i in 1:n
        for j in 1:n
            if type=="hard"
                ac2d[i][j] = HardThreshold(ac2d[i][j], t)
            elseif type=="soft"
                ac2d[i][j] = SoftThreshold(ac2d[i][j], t)
            end
        end
    end
    return ac2d
end

end # module
