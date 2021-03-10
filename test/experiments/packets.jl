using Wavelets, AutocorrelationShell, Plots, LinearAlgebra, AbstractTrees

## base node methods
Q = qfilter(wavelet(WT.db2));
P = pfilter(wavelet(WT.db2));

# One hot signal
x = zeros(256); x[128] = 1; # One hot signal
decomp1 = acwt(x, L=1, P=P, Q=Q);
tree1 = acwpt(x, P, Q);

plot(decomp1[:,1])
plot(collect(PostOrderDFS(tree1))[2].data)

acwptPostOrderBestBasis(tree1)

# Extract level 4 decomp
X = decomp1[:,4];
X = acwt(x, L=2, P=P, Q=Q)[:,4];

decomp2 = acwt(X, L=0, P=P, Q=Q);
tree2 = acwpt(X, P, Q);

plot(X)
plot(collect(PostOrderDFS(tree2))[1].data)
plot(decomp2[:,1])

tree2_pruned = acwptPostOrderBestBasis(tree2)

collect(PostOrderDFS(tree2_pruned))

# Check if bottom left node is the same
# norm(collect(PostOrderDFS(tree2))[1].data - decomp2[:,1]) < 1e-15
norm(acwt(x, L=0, P=P, Q=Q)[:,1] - collect(PostOrderDFS(acwpt(x, P, Q)))[1].data) == 0
