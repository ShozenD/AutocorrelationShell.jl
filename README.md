# AutocorrelationShell.jl

This library is an implementation of autocorrelation wavelets in Julia. It was
made by Rishi Subramanian, Christina Chang, and Shozen Dan under the supervision of Professor Naoki Saito at UC Davis.

## Dependencies
The required packages are
+ `AbstractTrees`
+ `DSP`
+ `Reexport`
+ `SpecialFunctions`
+ `StatsBase`
+ `Wavelets`

## Usage
Load the autocorrelation module
```{julia}
include("./src/AutocorrelationShell.jl")
using Main.AutocorrelationShell
```

## 1D Autocorrelation Wavelet Transform
```{julia}
# Forward Autocorrelation Wavelet Transform
fwt_ac(x,L,P,Q)

# Inverse Autocorrelation Wavelet Transform
iwt_ac(decomp)
```

### Example
Load the `Plots` package and the `wiggle` function
```{julia}
using Plots
include("./test/wiggle.jl")
```

Perform forward autocorrelation wavelet transform on the vector x
```{julia}
H = wavelet(WT.db2)
L = 2
Q = qfilter(H)
P = pfilter(H)

x = zeros(256)
x[128] = 1

decomposition = fwt_ac(x,L,P,Q)

wiggle(decomposition, Overlap = false)
```

### Result:

![Result](Presentations/2019/Overleaf/auto_decomposition.png)

## 2D Autocorrelation Wavelet Transform
```{julia}
# Forward Autocorrelation Wavelet Transform
ac2d(img,L_row,L_col,P,Q)
```
The `ac2d` function performs a forward wavelet transformation on 2D signals such as images. It returns a 4 dimensional tensor(multidimensional array) with the dimensions (num_row, num_col, levels_of_decomp_row, levels_of_decomp_col).

![AC2D transform example](Presentations/ac2d_decomp_heatmap.png)

```{julia}
# Inverse Autocorrelation Wavelet Transform
iac2d(decomp)
```
The `iac2d` function is the opposite of the `ac2d` function. It takes a transformed signal (i.e. the output of `ac2d`) and reverts it to the original signal.

### Example
```{julia}
H = wavelet(WT.db2)
L_row = 2
L_col = 2
Q = qfilter(H)
P = pfilter(H)

img = load(../test/pictures/lenna.jpg)
img = Float64.(Gray.(img))

decomposition = ac2d(img,L_row,L_col,P,Q)

# Display the 6th row and column decomposition
ac2d_heatmap(decomposition[:,:,6,6])

# Revert to original signal
reconstruct = iac2d(decomposition)
```

## Autocorrelation Wavelet Packet Transform
```{julia}
# Autocorrelation Wavelet Packets Transform
acwpt(x, P, Q)
```
The `acwpt` function computes the autocorrelation wavelet packet transform for 1 dimensional signal. It returns a binary tree object where the root node contains the original signal, and each child node contains a vector of 1 dimensional autocorrelation wavelet transform coefficients.

![AC Wavelet Packet Transform Diagram](Presentations/acwpt_diagram.png)

### Example
```{julia}
H = wavelet(WT.db2)
L_row = 2
L_col = 2
Q = qfilter(H)
P = pfilter(H)

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