```@meta
CurrentModule = AutocorrelationShell
DocTestSetup = quote
  using AutocorrelationShell
end
```

```@docs
acwpt(x::Vector{T}, P::Vector{T}, Q::Vector{T}) where T<:Real

aciwpt(tree::BinaryNode)

acwptBestBasisTree(node::BinaryNode; direction::AbstractString="right", et::Wavelets.Entropy=NormEntropy())
```
