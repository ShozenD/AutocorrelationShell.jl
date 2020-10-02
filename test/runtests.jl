using Test, AutocorrelationShell, Wavelets, LinearAlgebra, AbstractTrees

## Autocorrelation Wavelet Transform 
@test begin
    Q = qfilter(wavelet(WT.db2));
    P = pfilter(wavelet(WT.db2));

    x = zeros(256); x[128] = 1;
    decomp = acwt(x, L=2, P=P, Q=Q)

    norm(x - iacwt(decomp)) < 1e-15
end

## Autocorrelation Wavelet Packet Transform
@test begin
    Q = qfilter(wavelet(WT.db2));
    P = pfilter(wavelet(WT.db2));

    x = zeros(256); x[128] = 1;
    X = acwt(x, L=2, P=P, Q=Q)[:,4];

    decomp = acwt(X, L=0, P=P, Q=Q);
    tree = acwpt(X, P, Q);

    norm(collect(PostOrderDFS(tree))[1].data - decomp[:,1]) < 1e-15
end
