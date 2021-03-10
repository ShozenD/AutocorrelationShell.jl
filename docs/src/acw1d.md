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

## Forward AC Wavelet Transform
To perform the transform on a signal, use the `acwt` function.

```@docs
acwt
```

## Inverse AC Wavelet Transform
To reconstruct the signal from a array of AC Wavelet coefficients, use the `iacwt` function.

```@docs
iacwt
```
