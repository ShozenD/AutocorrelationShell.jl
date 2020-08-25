```@meta
CurrentModule = AutocorrelationShell
DocTestSetup = quote
  using AutocorrelationShell
end
```

```@docs
pfilter(filter::OrthoFilter)

qfilter(filter::OrthoFilter)

fwt_ac(x::Vector{T}, L::Integer, P::Vector{T}, Q::Vector{T}) where T<:Number

acwt(x::Vector{T}; L::Integer=1, P::Vector{T}, Q::Vector{T}) where T<:Real

iacwt(acwt::AbstractArray{<:Number})
```
