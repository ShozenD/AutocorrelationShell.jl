[![AutocorrelationShell.jl](figures/autocorrelation_shell_logo.png)](https://ShozenD.github.io/AutocorrelationShell.jl/stable)

[![CI](https://github.com/ShozenD/AutocorrelationShell.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/ShozenD/AutocorrelationShell.jl/actions)
[![](https://gitlab.com/BoundaryValueProblems/autocorrelation-shell/badges/master/pipeline.svg)](https://gitlab.com/BoundaryValueProblems/autocorrelation-shell/-/commits/master)
[![codecov](https://codecov.io/gh/ShozenD/AutocorrelationShell.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/ShozenD/AutocorrelationShell.jl)


This package is a [Julia](https://github.com/JuliaLang/julia) implementation of autocorrelation wavelets. The package includes the 1D autocorrelation wavelet transform, 2D autocorrelation wavelet transform, and autocorrelation wavelet packet transform.

Signal representations using autocorrelation wavelets are redundant and non-orthogonal. Some desirable properties of autocorrelation wavelet transforms are symmetry without losing vanishing moments, edge detection and characterization capabilities, and shift invariance. Autocorrelation wavelets can be used as a tool for data analysis such as time series analysis and image analysis.

## Authors
This package was first translated from Matlab code by Rishi Subramanian, and was extended by Christina Chang, and currently maintained by Shozen Dan under the supervision of Professor Naoki Saito at University of California, Davis.

## Installation
The package is part of the official Julia Registry. It can be install via the Julia REPL
```
(@1.x) pkg> add AutocorrelationShell
```
or
```
julia> using Pkg; Pkg.add("AutocorrelationShell")
```

## Usage
Load `AutocorrelationShell.jl` with the `Wavelets.jl` package
```{julia}
using Wavelets, AutocorrelationShell
```

## 1D Autocorrelation Wavelet Transform
```{julia}
# Forward 1D Autocorrelation Wavelet Transform
y = acwt(x, wavelet(WT.db4))

# Inverse 1D Autocorrelation Wavelet Transform
iacwt(y)
```

### Example
Perform forward autocorrelation wavelet transform on the vector x
```{julia}
x = zeros(256); x[128] = 1;

# Decompose signal
y = acwt(x, wavelet(WT.db4))

# Display decomposition
wiggle(y)
```

Result:

![Result](figures/auto_decomposition.png)

## 2D Autocorrelation Wavelet Transform
```{julia}
# Forward 2D Autocorrelation Wavelet Transform
y = acwt(img, wavelet(WT.db4))
```
The `acwt` function performs a forward wavelet transformation on 2D signals such as images. It returns a 4 dimensional tensor with the dimensions (num_row, num_col, levels_of_decomp_row, levels_of_decomp_col).

<img src="figures/ac2d_decomp_heatmap.png" alt="AC2D transform example" width="600" />

```{julia}
# Inverse 2D Autocorrelation Wavelet Transform
iacwt(y)
```
The `iacwt` function is the inverse function of `acwt`. It takes an array of autocorrelation wavelet coefficients and reconstructs the original signal.

### Example
```{julia}
X = load(../test/pictures/boat.jpg)
X = Float64.(Gray.(X))

Y = acwt(X, wavelet(WT.db4))

# Revert to original signal
Z = iacwt(Y)
```

## Autocorrelation Wavelet Packet Transform
```{julia}
# Autocorrelation Wavelet Packet Transform
acwpt(x, wavelet(WT.db4))
```
The `acwpt` function computes the autocorrelation wavelet packet transform for 1 dimensional signal. It returns a binary tree where the root node contains the original signal, and each child node contains a vector of 1 dimensional autocorrelation wavelet transform coefficients.

<img src="figures/acwpt_diagram.png" alt="AC Wavelet Packet Transform Diagram" width="600" />

### Example
```{julia}
using Wavelets

X₁ = randn(4); # length 4 random signal
y = acwpt(X₁, wavelet(WT.db4))
```
