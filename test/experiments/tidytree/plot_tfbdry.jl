using Plots, LinearAlgebra, AutocorrelationShell, Wavelets, AbstractTrees

bbpattern = zeros(8,4);

bbpattern[5:8,2] .= 1.0;
bbpattern[3:4,3] .= 1.0;
bbpattern[1:2,4] .= 1.0;


function plot_tfbdry(x, c)
    heatmap(transpose(bbpattern), yflip=true, legend=false)
    nrow, ncol = size(bbpattern)
    y_coords = [i + 0.5 for i in 1:(ncol-1)];
    hline!(y_coords, color=:red)

    for i in 1:ncol
        for j in 1:2^(i-1)
            vpos = (nrow/2^i)*(2*j-1) + 0.5
            plot!(vpos*ones(ncol-i+1), (i+0.5):(ncol+0.5), color=:red)
        end
    end
    current()
end

plot_tfbdry(bbpattern, :red)

Q = qfilter(wavelet(WT.db2));
P = pfilter(wavelet(WT.db2));

# One hot signal
x = zeros(256); x[128] = 1; # One hot signal
decomp1 = acwt(x, L=1, P=P, Q=Q);
tree1 = acwpt(x, P, Q);

# Extract level 4 decomp
X = decomp1[:,4];
X = acwt(x, L=2, P=P, Q=Q)[:,4];

decomp2 = acwt(X, L=0, P=P, Q=Q);
tree2 = acwpt(X, P, Q);

tree2 = acwptPostOrderBestBasis(tree2)

collect(PreOrderDFS(tree2))


ncol = dyadlength(tree2.data)
nrow = length(tree2.data)

bbpattern = zeros(nrow, ncol)

bbpattern[:,1]
index_hash = zeros(ncol)

for node in collect(PreOrderDFS(tree2)):
    if node.depth > 0:
        index_hash[node.depth] += length(node.data)/2^node.depth
    end
end
