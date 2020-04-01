# AutocorrelationShell.jl

This library is an implementation of autocorrelation wavelets in Julia. It was
made by Rishi Subramanian, Christina Chang, and Shozen Dan under the supervision of Professor Naoki Saito at UC Davis.

## Dependencies
The required packages are
+ `SpecialFunctions`
+ `Reexport`
+ `DSP`
+ `StatsBase`
+ `Wavelets`

Here is an example of how to install a package via the package manager and load with `using`
```{julia}
Pkg.add("PackageName")
using PackageName
```

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
iwt_ac(y)
```

### Example
Load the `Plots` package and the `wiggle` function
```{julia}
using Plots
include("./wiggle.jl")
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

![Result](Presentation/auto_decomposition.png)

## 2D Autocorrelation Wavelet Transform
```{julia}
# Forward Autocorrelation Wavelet Transform
ac2d(img,L,P,Q)
```
The `ac2d` function performs a forward wavelet transformation on 2D signals such as images. It returns a 4 dimensional tensor(multidimensional array) with the dimensions (num_row, num_col, levels_of_decomp_row, levels_of_decomp_col).

```{julia}
# Inverse Autocorrelation Wavelet Transform
iac2d(decomp)
```
The `iac2d` function is the opposite of the `ac2d` function. It takes a transformed signal (i.e. the output of `ac2d`) and reverts it to the original signal.

### Example
```{julia}
H = wavelet(WT.db2)
L = 2
Q = qfilter(H)
P = pfilter(H)

img = load(../test/pictures/lenna.jpg)
img = Float64.(Gray.(img))

decomposition = ac2d(img,L,P,Q)

# Display the 6th row and column decomposition
ac2d_heatmap(decomposition[:,:,6,6])

# Revert to original signal
reconstruct = iac2d(decomposition)
```
