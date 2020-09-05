```@meta
CurrentModule = AutocorrelationShell
DocTestSetup = quote
  using AutocorrelationShell
end
```

# 2D Autocorrelation Wavelets Transform
The 2D AC Wavelet transform will take a 2D signal such as an image
and decompose it into wavelet coefficients. First, the row-wise coefficients will
be computed using the 1D transform and then the column wise coefficients for each set
of row-wise coefficients will be computed to create a 4 dimensional tensor.

## Forward AC Wavelet Transform
There are 2 functions available for the 2D decomposition. The `ac2d` function is the original implementation and is kept for replication purposes. We generally recommend
using the `acwt2D` function, which is a wrapper around `ac2d` because it has a better
syntax.

```@docs
ac2d(x::AbstractArray{T,2}, L_row::Integer, L_col::Integer, P::Vector{T}, Q::Vector{T}) where T<:Number

acwt2D(x::AbstractArray; L_row::Integer=1, L_col::Integer=1, P::Vector{T}, Q::Vector{T}) where T<:Number
```

## Inverse AC Wavelet Transform
To reconstruct a 2D signal from a 4 dimensional tensor of wavelet coefficients, use the
`acwt2D` function.

```@docs
iacwt2D(x::AbstractArray{<:Number})
```
