using Wavelets, AutocorrelationShell, Plots

function echant(n::Integer, d::Integer, b::Integer)
    return (b + 1):(2^d):n
end

## base node methods
function leftchild(data, parent::BinaryNode)
    !isdefined(parent, :left) || error("left child is already assigned")
    node = typeof(parent)(data, parent)
    parent.left = node
end

function rightchild(data, parent::BinaryNode)
    !isdefined(parent, :right) || error("right child is already assigned")
    node = typeof(parent)(data, parent)
    parent.right = node
end

## Wavelet Packets functions
function binary_ac1d_new(x::Vector{T}, node::BinaryNode,
    P::Vector{T}, Q::Vector{T}, depth::Integer) where T<:Real

    n = length(x)
    max_depth = dyadlength(x) # maximum tree depth
    left = zeros(n)
    right = zeros(n)

    if depth <= max_depth
        for b = 0:(2^(depth-2) - 1) # depth starts from 2
            s = x[echant(n, depth-2, b)]
            h = ac_filter(s, Q)
            l = ac_filter(s, P)
            left[echant(n, depth-2, b)] = l
            right[echant(n, depth-2, b)] = h
        end

        # left
        leftchild(left, node)
        binary_ac1d_new(left, node.left, P, Q, depth+1)

        # right
        rightchild(right, node)
        binary_ac1d_new(right, node.right, P, Q, depth+1)
    end
end

function acwpt_new(x::Vector{T}, P::Vector{T}, Q::Vector{T}) where T<:Real

    root = BinaryNode(x) # original signal
    binary_ac1d_new(x, root, P, Q, 2)

    return root
end

Q = qfilter(wavelet(WT.db2));
P = pfilter(wavelet(WT.db2));

x = zeros(256); x[128] = 1; # One hot signal
decomp = acwt(x, L=2, P=P, Q=Q)

X = decomp[:,4];

wentropy(X, NormEntropy())

tree = acwpt_new(X, P, Q);

# print_tree(tree) Too large to print

l = plot(1:256, tree.left.data, legend=false, ylims=(-0.3, 0.4))
wentropy(tree.left.data, NormEntropy())
png("/Users/shozendan/Desktop/left")

#tree.right.data
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

acwptPostOrderBestBasis(tree)
