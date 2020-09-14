"""
    make_noisy(x::AbstractArray{<:Number}, rng::MersenneTwister, a::Float64)

Adds gaussian white noise to an image.

# Arguments
`x::AbstractArray{<:Number}`: image
`rng::MersenneTwister`: random number generator ex) rng = MersenneTwister(123)
`a::Float64`: A scaling parameter to adjust level of noise.

# Example
```julia
using AutocorrelationShell, Images, FileIO, Random

# Load test image
img = load("./test/pictures/boat.jpg")
img = Float64.(Gray.(img))

rng = MersenneTwister(123)
noisy_image = make_noisy(img, rng, 0.7)
```
"""
function make_noisy(x::AbstractArray{<:Number}, rng::MersenneTwister, a::Float64)
    noisy = x + a * randn(rng, Float64, size(x))
    return noisy
end

"""
    snr(f::AbstractArray{T}, g::AbstractArray{T}) where T<:Number

Computes the signal to noise ratio between two signals. Returns signal to noise ratio(dB).

# Arguments
`f::AbstractArray{<:Number}`: Original signal
`g::AbstractArray{<:Number}`: Noisy signal
"""
function snr(f::AbstractArray{T}, g::AbstractArray{T}) where T<:Number
    L2_f = norm(f)
    L2_fg = norm(f-g)
    return 20*log10(L2_f/L2_fg)
end

"""
    get_snr(y::AbstractArray{T, Integer}, x::AbstractArray{T, Integer}; type::AbstractString="hard", step::Float64=0.5) where T<:Number

Thresholds the wavelet coefficients matrix of a noisy signal, reconstructs it, then
computes the signal to noise ratio between the reconstructed signal and the
original non-noisy signal. Returns a list of the ratio of non-zero coefficient at
each threshold and a list containing the signal to noise ratio.

# Arguments
- `y`: original 2D signal
- `x`: wavelet coefficient matrix of a noisy 2D signal
- `type`: type of thresholding (soft or hard)
- `step`: step at which to increase the threshold

# Example
```julia
img = load("./test/pictures/lenna.jpg")
img = Float64.(Gray.(img))

# Add noise to image
rng = MersenneTwister(123)
noisy = make_noisy(img, rng, 0.7)

# Apply Wavelet transform
H = wavelet(WT.db2)
L = 2
Q = qfilter(H)
P = pfilter(H)

ac_noisy = ac2d(noisy, L, P, Q)

coef_ratio, snr_list = get_snr(img, ac_noisy, type="soft", step=0.5)
```
"""
function get_snr(y::AbstractArray{T, Integer}, x::AbstractArray{T, Integer}; type::AbstractString="hard", step::Float64=0.5) where T<:Number
    max_coef = maximum(abs.(x)) # find largest coefficient
    num_coef = length(x)
    snr_list = zeros(0)
    coef_ratio_list = zeros(0)
    @inbounds begin
        for i in range(0, stop=round(max_coef), step=step)
            thresh = acthreshold(x, type, i)
            reconst = iac2d(thresh)
            num_nonzero_coef = sum(abs.(thresh) .> 0)
            append!(snr_list, snr(y, reconst))
            append!(coef_ratio_list, num_nonzero_coef/num_coef)
            println(i/round(max_coef) * 100)
        end
    end
    return coef_ratio_list, snr_list
end

"""
    get_ssim(y::AbstractArray{T}, x::AbstractArray{T}; type::AbstractString="hard", step::Float64=0.5) where T<:Number

Thresholds the wavelet coefficients matrix of a noisy signal, reconstructs it, then
computes the structural similarity index between the reconstructed signal and the
original non-noisy signal. Returns a list of the ratio of non-zero coefficient at
each threshold and a list containing the SSIM values.

# Arguments
- `y`: original 2D signal
- `x`: wavelet coefficient matrix of a noisy 2D signal
- `type`: type of thresholding (soft or hard)
- `step`: step at which to increase the threshold

# Example
```julia
using AutocorrelationShell, Wavelets, FileIO, Images, ImageQualityIndexes, Random

img = load("./test/pictures/lenna.jpg")
img = Float64.(Gray.(img))

# Add noise to image
rng = MersenneTwister(123)
noisy = make_noisy(img, rng, 0.7)

# Apply Wavelet transform
H = wavelet(WT.db2)
L = 2
Q = qfilter(H)
P = pfilter(H)

ac_noisy = acwt2d(noisy; L=2, P=P, Q=Q)

coef_ratios, ssim_values = get_ssim(img, ac_noisy; type="hard", step=0.5)
```
"""
function get_ssim(y::AbstractArray{T}, x::AbstractArray{T}; type::AbstractString="hard", step::Float64=0.5) where T<:Number
    max_coef = maximum(abs.(x)) # find largest coefficient
    num_coef = length(x)
    ssim_list = zeros(0)
    coef_ratio_list = zeros(0)
    @inbounds begin
        for i in range(0, stop=round(max_coef), step=step)
            thresh = acthreshold(x, type, i)
            reconst = iac2d(thresh)
            num_nonzero_coef = sum(abs.(thresh) .> 0)
            append!(ssim_list, assess(SSIM(), y, reconst))
            append!(coef_ratio_list, num_nonzero_coef/num_coef)
            println(i/round(max_coef) * 100)
        end
    end
    return coef_ratio_list, ssim_list
end

"""
    acwt_heatmap(x::AbstractArray{<:Number})

Takes the natural log of the absolute value of each cell, and plots a heatmap.
Can be used for visualizing an input image or the coefficient matrix of a 2D ACW decomposition.

# Arguments
- `x::AbstractArray{<:Number}`: A image or a 2D matrix
"""
function acwt_heatmap(x::AbstractArray{<:Number})
    heatmap(log.(abs.(x)),
            yflip=true,
            axis=nothing,
            colorbar_entry=false,
            aspect_ratio=:equal,
            showaxis=false)
end

## Threshold functions
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
    acthreshold(x::AbstractArray{<:Number}, type::AbstractString, t::Float64)

Thresholds a given array of coeffecients using either hard thresholding or soft thresholding.

# Arguments
- `x::AbstractArray{<:Number}`: Autocorrelation wavelet coefficient array
- `type::AbstractString`: Threshold type. "hard" or "soft"
- `t::Float64`: threshold value
"""
function acthreshold(x::AbstractArray{<:Number}, type::AbstractString, t::Float64)
    y = deepcopy(x) # to prevent inplace behavior but slows speed
    if type=="hard"
        HardThreshold!(y, t)
    elseif type=="soft"
        SoftThreshold!(y, t)
    end
    return y
end

## Entropy functions
# ENTROPY TYPES

# Extend Wavetes.Entropy to include norm and threshold entropy methods
struct NormEntropy <: Wavelets.Entropy end
struct ThresholdEntropy <: Wavelets.Entropy end # Should we implement?

# Entropy measures: Additive with wentropy(0) = 0
# all coefs assumed to be on [-1,1] after normalization with nrm
# given x and y, where x has "more concentrated energy" than y
# then coefentropy(x, et, norm) <= coefentropy(y, et, norm) should be satisfied.

# (non-normalized) Shannon Entropy
function wentropy(x::T, et::Wavelets.ShannonEntropy, nrm::T)  where T<:Number
    s = (x/nrm)^2
    if s == 0.0
        return -zero(T)
    else
        return -s*log(s)
    end
end

# log energy entropy
function wentropy(x::T, et::Wavelets.LogEnergyEntropy, nrm::T) where T<:Number
    s = (x/nrm)^2
    if s == 0.0
        return -zero(T)
    else
        return -log(s)
    end
end

function wentropy(x::T, et::NormEntropy, nrm::T; p::T=1) where T<:Number
    s = abs((x/nrm))
    return s^p
end

"""
    wentropy(x::AbstractArray{T}, et::Wavelets.Entropy, nrm::T=norm(x)) where T<:Number

Calculates the entropy of a given vector after normalizing it with the l2 norm.

# Arguments
- `x::AbstractArray{T}`: vector
- `et::Wavelets.Entropy`: Type of entropy. Available methods are ShannonEntropy(), LogEnergyEntropy(), NormEntropy()
- `nrm::T=norm(x)`: The norm of the given vector. *default*: L2 norm

# Examples
```{julia}
using AutocorrelationShell, Wavelets

wentropy(x, ShannonEntropy()) # Shannon Entropy

wentropy(x, LogEnergyEntropy()) # Log-energy Entropy

wentropy(x, NormEntropy()) # Norm Entropy
```
"""
function wentropy(x::AbstractArray{T}, et::Wavelets.Entropy, nrm::T=norm(x)) where T<:Number
    @assert nrm >= 0
    sum = zero(T)
    nrm == sum && return sum
    for i in eachindex(x)
        @inbounds sum += wentropy(x[i], et, nrm)
    end
    return sum
end

function wentropy(x::AbstractArray{T}, et::NormEntropy, nrm::T=norm(x); p=1.0) where T<:Number
    @assert nrm >= 0
    sum = zero(T)
    nrm == sum && return sum
    for i in eachindex(x)
        @inbounds sum += wentropy(x[i], et, nrm, p=p)
    end
    return sum
end

## Visualizations
"""
    wiggle(wav; taxis=1:size(wav,1), zaxis=1:size(wav,2), sc=1, EdgeColor=:black, FaceColor=:black, Orient=:across, Overlap=true, ZDir=:normal)

Plot a set of shaded wiggles

# Arguments
- `wav::Array`: matrix of waveform columns.
- `taxis::Array=1:size(wav,1)`: time axis vector
- `zaxis::Array=1:size(wav,2)`: space axis vector
- `sc::Float=1`: scale factor/magnification.
- `EdgeColor::Symbol=:black`: Sets edge of wiggles color.
- `FaceColor::Symbol=:black`: Sets shading color of wiggles.
- `Overlap::bool=true`: How signals are scaled.
        true  - Signals overlap (default);
        false - Signals are scaled so they do not overlap.
- `Orient::Symbol=:across`: Controls orientation of wiggles.
        :across - from left to right
        :down   - from top to down
- `ZDir::Symbol=:normal`: Direction of space axis.
        :normal  - First signal at bottom (default)
        :reverse - First signal at top.

Translated by Nicholas Hausch -- MATLAB file provided by Naoki Saito
The previous MATLAB version contributors are:
 Anthony K. Booer (SLB) and Bradley Marchand (NSWC-PC)
Revised by Naoki Saito, Feb. 05, 2018
"""
function wiggle(wav; taxis=1:size(wav,1), zaxis=1:size(wav,2), sc=1, EdgeColor=:black, FaceColor=:black, Overlap=true, Orient=:across, ZDir=:normal)

    # Set axes
    (n,m) = size(wav)

    # Sanity check
    if length(taxis) != n
        error("Inconsistent taxis dimension!")
    end
    if length(zaxis) != m
        error("Inconsistent zaxis dimension!")
    end

    # For calculation purposes
    maxrow = zeros(m); minrow = zeros(m)
    for k = 1:m
      maxrow[k] = maximum(wav[:,k]); minrow[k] = minimum(wav[:,k])
    end

    # Scale the data for plotting
    wamp = deepcopy(wav)
    dt = mean(diff(taxis))
    dz = mean(diff(zaxis))
    if Overlap
      wamp *= 2 * dz * (sc/maximum(maxrow-minrow))
    else
      wmax = maximum(maxrow); wmin = minimum(minrow);
      if wmax<=0
        wmax = 0
      end
      if wmin>=0
        wmin = 0
      end
        wamp = sc*wav/(wmax-wmin)
    end

    # Set initial plot
    t0 = minimum(taxis)
    t1 = maximum(taxis)
    z0 = minimum(zaxis)
    z1 = maximum(zaxis)
    if Orient == :down
     plot(xlims=(z0-dz,z1+dz), ylims=(t0,t1), yflip=true, legend=:none)
    else
     plot(xlims=(t0,t1), ylims=(z0-dz,z1+dz), legend=:none)
    end
    if ZDir == :reverse
        wamp = flipdim(wamp,2)
    end

    # Plot each wavelet
    for k = 1:m
      sig = wamp[:,k]
      t = deepcopy(taxis)
      w_sign = sign.(sig)
      for j=1:n-1
        if (w_sign[j]!=w_sign[j+1] && w_sign[j]!=0 && w_sign[j+1]!=0)
          sig = [sig; 0]
          t = [t; t[j]-sig[j]*(t[j+1]-t[j])/(sig[j+1]-sig[j])]
        end
      end
      IX = sortperm(t)
      t = t[IX]
      sig = sig[IX]
      len = length(t)
      len1 = collect(len:-1:1)
      indperm = [1:len;len1]
      inputx = t[indperm]
      inputy = zaxis[k] .+ [sig;min.(sig[len1],0)]
        # In the plot! functions below, theoretically speaking, either
        # fillrange = zaxis[k] or fillrange=[zaxis[k], zaxis[k]+dz] should be used.
        # However, those do not generate the desired plots as of O3/19/2018.
        # Somehow, the relative value of 0, i.e., fillrange=0, works well,
        # which is used temporarily.
      if Orient == :down
        plot!(inputy, inputx, fillrange=0, fillalpha=0.75, fillcolor=FaceColor, linecolor=EdgeColor, orientation=:v)
      else
        plot!(inputx, inputy, fillrange=0, fillalpha=0.75, fillcolor=FaceColor, linecolor=EdgeColor)
      end
    end
    plot!() # flushing the display.
end

"""
    wiggle!(wav; taxis=1:size(wav,1), zaxis=1:size(wav,2), sc=1, EdgeColor=:black, FaceColor=:black, Orient=:across, Overlap=true, ZDir=:normal)

Plot a set of shaded wiggles on the current displayed graphics

# Arguments
- `wav::Array`: matrix of waveform columns.
- `taxis::Array=1:size(wav,1)`: time axis vector
- `zaxis::Array=1:size(wav,2)`: space axis vector
- `sc::Float=1`: scale factor/magnification.
- `EdgeColor::Symbol=:black`: Sets edge of wiggles color.
- `FaceColor::Symbol=:black`: Sets shading color of wiggles.
- `Overlap::bool=true`: How signals are scaled.
        true  - Signals overlap (default);
        false - Signals are scaled so they do not overlap.
- `Orient::Symbol=:across`: Controls orientation of wiggles.
        :across - from left to right
        :down   - from top to down
- `ZDir::Symbol=:normal`: Direction of space axis.
        :normal  - First signal at bottom (default)
        :reverse - First signal at top.

Translated by Nicholas Hausch -- MATLAB file provided by Naoki Saito
The previous MATLAB version contributors are:
 Anthony K. Booer (SLB) and Bradley Marchand (NSWC-PC)
Revised by Naoki Saito, Feb. 05, 2018
"""
function wiggle!(wav; taxis=1:size(wav,1), zaxis=1:size(wav,2), sc=1, EdgeColor=:black, FaceColor=:black, Overlap=true, Orient=:across, ZDir=:normal)

    # Set axes
    (n,m) = size(wav)

    # Sanity check
    if length(taxis) != n
        error("Inconsistent taxis dimension!")
    end
    if length(zaxis) != m
        error("Inconsistent zaxis dimension!")
    end

    # For calculation purposes
    maxrow = zeros(m); minrow = zeros(m)
    for k = 1:m
      maxrow[k] = maximum(wav[:,k]); minrow[k] = minimum(wav[:,k])
    end

    # Scale the data for plotting
    wamp = deepcopy(wav)
    dt = mean(diff(taxis))
    dz = mean(diff(zaxis))
        if Overlap
      wamp *= 2 * dz * (sc/maximum(maxrow-minrow))
    else
      wmax = maximum(maxrow); wmin = minimum(minrow);
      if wmax<=0
        wmax = 0
      end
      if wmin>=0
        wmin = 0
      end
        wamp = sc*wav/(wmax-wmin)
    end

    # Set initial plot
    t0 = minimum(taxis)
    t1 = maximum(taxis)
    z0 = minimum(zaxis)
    z1 = maximum(zaxis)
    if Orient == :down
     plot!(xlims=(z0-dz,z1+dz), ylims=(t0,t1), yflip=true, legend=:none)
    else
     plot!(xlims=(t0,t1), ylims=(z0-dz,z1+dz), legend=:none)
    end
    if ZDir == :reverse
        wamp = flipdim(wamp,2)
    end

    # Plot each wavelet
    for k = 1:m
      sig = wamp[:,k]
      t = deepcopy(taxis)
      w_sign = sign.(sig)
      for j=1:n-1
        if (w_sign[j]!=w_sign[j+1] && w_sign[j]!=0 && w_sign[j+1]!=0)
          sig = [sig; 0]
          t = [t; t[j]-sig[j]*(t[j+1]-t[j])/(sig[j+1]-sig[j])]
        end
      end
      IX = sortperm(t)
      t = t[IX]
      sig = sig[IX]
      len = length(t)
      len1 = collect(len:-1:1)
      indperm = [1:len;len1]
      inputx = t[indperm]
      inputy = zaxis[k] .+ [sig;min.(sig[len1],0)]
        # In the plot! functions below, theoretically speaking, either
        # fillrange = zaxis[k] or fillrange=[zaxis[k], zaxis[k]+dz] should be used.
        # However, those do not generate the desired plots as of O3/19/2018.
        # Somehow, the relative value of 0, i.e., fillrange=0, works well,
        # which is used temporarily.
      if Orient == :down
        plot!(inputy, inputx, fillrange=0, fillalpha=0.75, fillcolor=FaceColor, linecolor=EdgeColor, orientation=:v)
      else
        plot!(inputx, inputy, fillrange=0, fillalpha=0.75, fillcolor=FaceColor, linecolor=EdgeColor)
      end
    end
    plot!() # flushing the display.
end
