# Getting Started

In this section we will provide a condensed overview of the
package. In order to keep this overview concise, we will not
discuss any background information or theory here in detail.

## Installation

To install
[AutocorrelationShell.jl](https://gitlab.com/BoundaryValueProblems/autocorrelation-shell), start up
Julia and type the following code-snipped into the REPL. It makes
use of the native Julia package manger.

```julia
Pkg.add("AutocorrelationShell")
```

Additionally, for example if you encounter any sudden issues, or
in the case you would like to contribute to the package, you can
manually choose to be on the latest (untagged) version.

```julia
Pkg.checkout("AutocorrelationShell")
```

## Examples

### 1D Autocorrelation Wavelet Transform
The following code snippet shows how to obtain the autocorrelation wavelet decomposition of a 1D signal.

```julia
using AutocorrelationShell, Wavlets, Plots

H = wavelet(WT.db2)
L = 2
Q = qfilter(H)
P = pfilter(H)

x = zeros(256)
x[128] = 1

decomposition = acwt(x, L=2, P=P, Q=Q)

wiggle(decomposition, Overlap = false)
```

### 2D Autocorrelation Wavelet Transform
```julia
acwt2D(x; L_row, L_col, P, Q)
```
The `acwt2D` function performs a forward wavelet transformation on 2D signals such as images. It returns a 4 dimensional tensor with the dimensions (num_row, num_col, levels_of_decomp_row, levels_of_decomp_col).

```julia
iacwt2D(x)
```
The `iacwt2D` function is the inverse function of `acwt2D`. It takes the output of `acwt2D`(i.e. the wavelet coefficient matrix) and reconstructs the original signal.

The following code snippet shows how to obtain the autocorrelation wavelet decomposition of an image.

```julia
H = wavelet(WT.db2)
Q = qfilter(H)
P = pfilter(H)

img = load(../test/pictures/boat.jpg)
img = Float64.(Gray.(img))

decomposition = acwt2D(img, L_row=2, L_col=2, P=P, Q=Q)

# Display the 6th row and column decomposition
acwt_heatmap(decomposition[:,:,6,6])

# Revert to original signal
reconstruct = iacwt2D(decomposition)
```

### 1D Autocorrelation Wavelet Packet Transform
```julia
acwpt(x, P, Q)
```
The `acwpt` function computes the autocorrelation wavelet packet transform for 1 dimensional signal. It returns a binary tree object where the root node contains the original signal, and each child node contains a vector of 1 dimensional autocorrelation wavelet transform coefficients.

The following code snippet shows how to obtain the autocorrelation wavelet packet transformation of a 1D signal.

```julia
using Random, Wavelets, AbstractTrees
rng = MersenneTwister(123);

X₁ = randn(rng, 4); # length 4 random signal
H = wavelet(WT.db2);
Q = qfilter(H);
P = pfilter(H);
decomp = acwpt(X₁, P, Q)

# Print the tree in the console
print_tree(decomp)

# Gather all nodes into a vector
collect(PostOrderDFS(decomp))
```

## Getting Help

To get help on specific functionality you can either look up the
information here, or if you prefer you can make use of Julia's
native doc-system.

If you encounter a bug or would like to participate in the
development of this package come find us on GitLab.

- [BoundaryValueProblems/AutocorrelationShell.jl](https://gitlab.com/BoundaryValueProblems/autocorrelation-shell)
