```@meta
CurrentModule = AutocorrelationShell
DocTestSetup = quote
  using AutocorrelationShell
end
```

## Auto-correlation Wavelets Utilities
AutocorrelationShell.jl contains functions that can be used for adding
gaussian noise, computing signal to noise ratio, computing structural similarity index(SSIM), performing thresholding, visualizing matrix heat map, visualizing the decomposition of the 1D transform, or calculate the entropy of a vector.

```@docs
make_noisy(x::AbstractArray{<:Number}, rng::MersenneTwister, a::Float64)

snr(f::AbstractArray{T}, g::AbstractArray{T}) where T<:Number

acthreshold(x::AbstractArray{<:Number}, type::AbstractString, t::Float64)

get_snr(y::AbstractArray{T, Integer}, x::AbstractArray{T, Integer}; type::AbstractString="hard", step::Float64=0.5) where T<:Number

get_ssim(y::AbstractArray{T}, x::AbstractArray{T}; type::AbstractString="hard", step::Float64=0.5) where T<:Number

acwt_heatmap(x::AbstractArray{<:Number})

wentropy(x::AbstractArray{T}, et::Wavelets.Entropy, nrm::T=norm(x)) where T<:Number

wiggle(wav; taxis=1:size(wav,1), zaxis=1:size(wav,2), sc=1, EdgeColor=:black, FaceColor=:black, Overlap=true, Orient=:across, ZDir=:normal)

wiggle!(wav; taxis=1:size(wav,1), zaxis=1:size(wav,2), sc=1, EdgeColor=:black, FaceColor=:black, Overlap=true, Orient=:across, ZDir=:normal)
```
