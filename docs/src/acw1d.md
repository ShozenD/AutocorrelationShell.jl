```@meta
CurrentModule = AutocorrelationShell
DocTestSetup = quote
  using AutocorrelationShell
end
```

# 1D Autocorrelation Wavelets Transform
The 1D AC Wavelet transform will take a 1D signal such as a single time series
and decompose it into wavelet coefficients.

## Orthogonal Filters
To begin, one must first specify the high pass and low pass autocorrelation filters
to use in the decomposition. The types of filters available for use are `Haar`, `Coiflet`,
`Daubechies`, `Symlet`, `Battle`, `Beylkin`, `Vaidyanathan`, and `CDF`. These filters are
implemented within the [Wavelets.jl package](https://github.com/JuliaDSP/Wavelets.jl).

```@docs
pfilter(filter::OrthoFilter)

qfilter(filter::OrthoFilter)
```

## Forward AC Wavelet Transform
To perform the transform on a signal, one can use either the `fwt_ac` of `acwt` function.

The `fwt_ac` function is the original transform function, translated from matlab code by Rishi Subramanian and
Christina Chang in 2019. The `acwt` function is a wrapper around `fwt_ac` that improves the syntax. Generally, we
recommend using the `acwt` function.

```@docs
fwt_ac(x::Vector{T}, L::Integer, P::Vector{T}, Q::Vector{T}) where T<:Number

acwt(x::Vector{T}; L::Integer=1, P::Vector{T}, Q::Vector{T}) where T<:Real
```

## Inverse AC Wavelet Transform
To reconstruct the signal from a array of AC Wavelet coefficients, use the `iacwt` function.

```@docs
iacwt(acwt::AbstractArray{<:Number})
```
