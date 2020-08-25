```@meta
CurrentModule = AutocorrelationShell
DocTestSetup = quote
  using AutocorrelationShell
end
```

```@docs
make_noisy(x::AbstractArray{<:Number}, rng::MersenneTwister, a::Float64)

snr(f::AbstractArray{T}, g::AbstractArray{T}) where T<:Number

acthreshold(x::AbstractArray{<:Number}, type::AbstractString, t::Float64)

get_snr(y::AbstractArray{T}, x::AbstractArray{T}; type::AbstractString="hard", step::Float64=0.5) where T<:Number

get_ssim(y::AbstractArray{T}, x::AbstractArray{T}; type::AbstractString="hard", step::Float64=0.5) where T<:Number

acwt_heatmap(x::AbstractArray{<:Number})

wiggle(wav; taxis=1:size(wav,1), zaxis=1:size(wav,2), sc=1, EdgeColor=:black, FaceColor=:black, Overlap=true, Orient=:across, ZDir=:normal)

wiggle!(wav; taxis=1:size(wav,1), zaxis=1:size(wav,2), sc=1, EdgeColor=:black, FaceColor=:black, Overlap=true, Orient=:across, ZDir=:normal)
```
