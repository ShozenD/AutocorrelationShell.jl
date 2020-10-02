using Wavelets, AutocorrelationShell, Plots, LinearAlgebra, AbstractTrees

## base node methods
Q = qfilter(wavelet(WT.db2));
P = pfilter(wavelet(WT.db2));

x = zeros(256); x[128] = 1; # One hot signal
decomp = acwt(x, L=2, P=P, Q=Q)

X = decomp[:,4];

wentropy(X, NormEntropy())

tree = acwpt_new(X, P, Q);

decomp = fwt_ac_new(X, 0, P, Q);

norm(collect(PostOrderDFS(tree))[1].data - decomp[:,1]) < 1e-15

# print_tree(tree) Too large to print
l = plot(1:256, tree.left.data, legend=false, ylims=(-0.3, 0.4))
wentropy(tree.left.data, NormEntropy())
png("/Users/shozendan/Desktop/left")

r = plot(1:256, tree.right.data, legend=false)
r = plot(1:256, tree.right.data, legend=false, ylims=(-0.3, 0.4))
wentropy(tree.right.data, NormEntropy())
png("/Users/shozendan/Desktop/right")

rl = plot(1:256, tree.right.left.data, legend=false)
rl = plot(1:256, tree.right.left.data, legend=false, ylims=(-0.3, 0.4))
wentropy(tree.right.left.data, NormEntropy())
png("/Users/shozendan/Desktop/right-left")

rr = plot(1:256, tree.right.right.data, legend=false)
rr = plot(1:256, tree.right.right.data, legend=false, ylims=(-0.3, 0.4))
wentropy(tree.right.right.data, NormEntropy())
png("/Users/shozendan/Desktop/right-right")

ll = plot(1:256, tree.left.left.data, legend=false, ylims=(-0.3, 0.4))
wentropy(tree.left.left.data, NormEntropy())
png("/Users/shozendan/Desktop/left-left")

lr = plot(1:256, tree.left.right.data, legend=false, ylims=(-0.3, 0.4))
wentropy(tree.left.right.data, NormEntropy())
png("/Users/shozendan/Desktop/left-right")
