```@meta
CurrentModule = AutocorrelationShell
DocTestSetup = quote
  using AutocorrelationShell
end
```

```@docs
ac2d(x::AbstractArray{T,2}, L_row::Integer, L_col::Integer, P::Vector{T}, Q::Vector{T}) where T<:Number

acwt2D(x::AbstractArray; L_row::Integer=1, L_col::Integer=1, P::Vector{T}, Q::Vector{T}) where T<:Number

iacwt2D(x::AbstractArray{<:Number})
```
