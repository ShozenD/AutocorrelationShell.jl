# AutocorrelationShell.jl

This library is an implementation of autocorrelation wavelets in Julia. It was
made by Rishi Subramanian and Christina Chang under the supervision of Professor
 Naoki Saito at UC Davis.

```{julia}
using Plots
using AutocorrelationShell
include("./wiggle.jl")

H = wavelet(WT.db2)
L = 2
Q = qfilter(H)
P = pfilter(H)

x = zeros(256)
x[128] = 1

decomposition = fwt_ac(x,L,P,Q)

wiggle(decomposition, Overlap = false)
```

Result:

![Result](Presentation/auto_decomposition.png)