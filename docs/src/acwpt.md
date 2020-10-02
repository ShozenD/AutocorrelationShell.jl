```@meta
CurrentModule = AutocorrelationShell
DocTestSetup = quote
  using AutocorrelationShell
end
```

## Autocorrelation Wavelet Packet transform
The wavelet packet transform is a wavelet transform where both the approximation and
detail coefficients are passed through low and high pass quadrature mirror filters to
create a full binary tree. The main distinction between normal wavelet packets
and the autocorrelation wavelet packets is that the autocorrelation transform is
redundant. Therefore the number of wavelet coefficients will double at each subsequent
level of the decomposition tree. To perform the autocorrelation wavelet transform, use the
`acwpt` function.

```@docs
acwpt(x::Vector{T}, P::Vector{T}, Q::Vector{T}) where T<:Real
```

By design, the original signal is contained in the root node of the binary decomposition tree.
However, One can reconstruct the signal from the wavelet coefficents, perhaps after pruning the tree, using the `aciwpt` function.
```@docs
iacwpt(tree::BinaryNode)
```

For a wavelet packet decomposition, it is interesting to find an optimal decomposition with
respect to a convenient criterion. For a more detailed explanation on how an optimal decomposition is chosen for a wavelet packet decomposition, refer to the **Choosing the Optimal Decomposition** section of the [Wavelet Packets](https://www.mathworks.com/help/wavelet/ug/wavelet-packets.html) documentation on the MathWorks website.

The `acwptPostOrderBestBasis` traverses the binary tree in an bottom-up order and is therefore the most democratic way of choosing the optimal wavelet coefficients. On the otherhand, `acwptPreOrderBestBasis` will traverse the binary tree in an top-bottom order and
is therefore a "greedy" way of choosing the optimal wavelet coefficients. The two function
will often return a different binary tree for the same decomposition.

We use different entropy methods as our criteria because the describe information-related properties for an accurate representation of a given signal. The entropy methods currently available within our package are, `ShannonEntropy`, `LogEnergyEntropy`, and `NormEntropy`. For a more detailed description of each method, refer to the [AC Wavelet Utils documentation](https://boundaryvalueproblems.gitlab.io/autocorrelation-shell/acwutil/).

```@docs
acwptPostOrderBestBasis(tree::BinaryNode; et::Wavelets.Entropy=NormEntropy())

acwptPreOrderBestBasis(tree::BinaryNode; et::Wavelets.Entropy=NormEntropy())
```
