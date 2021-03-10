using Wavelets, Plots, LinearAlgebra, AutocorrelationShell

X = zeros(128)
X[64] = 1

wt = wavelet(WT.haar)
d = dwt(X, wt, 4)

plot(X)
plot(d)

P = qfilter(wavelet(WT.haar))
Q = pfilter(wavelet(WT.haar))

acwt(X, P=P, Q=Q)

n = 2^8
x0 = testfunction(n,"Doppler")
x = x0 + 0.05*randn(n)

acwt(x, P=P, Q=Q)
