using Test, AutocorrelationShell, Wavelets, LinearAlgebra

@test begin
    Q = qfilter(wavelet(WT.db2));
    P = pfilter(wavelet(WT.db2));

    x = zeros(256); x[128] = 1;

    decomp = acwt(x, L=2, P=P, Q=Q)

    norm(x - iacwt(decomp)) < 1e-14
end
