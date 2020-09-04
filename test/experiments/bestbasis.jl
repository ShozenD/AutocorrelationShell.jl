using AutocorrelationShell, Wavelets, AbstractTrees, Plots

Q = qfilter(wavelet(WT.db2));
P = pfilter(wavelet(WT.db2));

x = zeros(256); x[128] = 1; # One hot signal
plot(1:length(x), x, legend=false)
#png("/Users/shozendan/Desktop/onehot-signal")

decomp = acwt(x, L=2, P=P, Q=Q)
wiggle(decomp)
#png("/Users/shozendan/Desktop/wiggle")

X = decomp[:,4];
plot(1:length(X), X, legend=false)
png("/Users/shozendan/Desktop/coefficient-vector")
wentropy(X, NormEntropy())

tree = acwpt(X, P, Q);

# print_tree(tree) Too large to print

#tree.right.data
r = plot(1:256, tree.right.data, legend=false)
wentropy(tree.right.data, NormEntropy())
#png("/Users/shozendan/Desktop/right")

l = plot(1:256, tree.left.data, legend=false)
wentropy(tree.left.data, NormEntropy())
#png("/Users/shozendan/Desktop/left")

rr = plot(1:256, tree.right.right.data, legend=false)
wentropy(tree.right.right.data, NormEntropy())
#png("/Users/shozendan/Desktop/right-right")

rl = plot(1:256, tree.right.left.data, legend=false)
wentropy(tree.right.left.data, NormEntropy())
#png("/Users/shozendan/Desktop/right-left")

lr = plot(1:256, tree.left.right.data, legend=false)
wentropy(tree.left.right.data, NormEntropy())
#png("/Users/shozendan/Desktop/left-right")

ll = plot(1:256, tree.left.left.data, legend=false)
wentropy(tree.left.left.data, NormEntropy())
#png("/Users/shozendan/Desktop/left-left")

using LinearAlgebra

norm(tree.right.left.data - tree.left.right.data)

eroot = wentropy(tree.data, NormEntropy())

el = wentropy(tree.left.data, NormEntropy());
er = wentropy(tree.right.data, NormEntropy());
(el + er)/2 < eroot

ell = wentropy(tree.left.left.data, NormEntropy());
elr = wentropy(tree.left.right.data, NormEntropy());
(ell + elr)/2 < el

erl = wentropy(tree.right.left.data, NormEntropy());
err = wentropy(tree.right.right.data, NormEntropy());
(erl + err)/2 < er

wentropy([3.0, 3.0, 3.0], NormEntropy())

PostOrderDFS(tree)

## Post Order Best Basis
function acwptPostOrderBestBasis!(node::BinaryNode; direction::AbstractString="right", et::Wavelets.Entropy=NormEntropy())
    if isdefined(node, :left) & isdefined(node, :right)
        el = acwptPostOrderBestBasis!(node.left, direction="left")
        er = acwptPostOrderBestBasis!(node.right)
        if (el + er)/2 < wentropy(node.data, et)
            return (el + er)/2
        else isdefined(node, :parent) # handling cases where the first decomposition is unuseful
            data = node.data
            new_pruned_node = typeof(node)(data, node) # create node with no children
            if direction == "right"
                node.parent.right = new_pruned_node
            else
                node.parent.left = new_pruned_node
            end
            return wentropy(node.data, et)
        end
    else
        return wentropy(node.data, et) # return entropy value
    end
end

"""
    acwptPostOrderBestBasis(node::BinaryNode; et::Wavelets.Entropy)

Returns the best basis tree found using post order traversal. This is a democratic approach to finding the best basis tree.

# Arguments
-`tree::BinaryNode`: The entry point of the wavelet packet decomposition tree.
-`et::Wavelets.Entropy`: The type of cost function used for evaluation.
"""
function acwptPostOrderBestBasis(tree::BinaryNode; et::Wavelets.Entropy=NormEntropy())
    best_tree = deepcopy(tree)
    acwptPostOrderBestBasis!(best_tree, et=et)
    return best_tree
end

function acwptPreOrderBestBasis!(node::BinaryNode; direction::AbstractString="right", et::Wavelets.Entropy=NormEntropy())
    if isdefined(node, :left) & isdefined(node, :right)
        e0 = wentropy(node.data, et)
        e10 = wentropy(node.left.data, et)
        e11 = wentropy(node.right.data, et)
        if (e10 + e11)/2 < e0
            acwptPreOrderBestBasis!(node.left, direction="left", et=et)
            acwptPreOrderBestBasis!(node.right, et=et)
        else
            data = node.data
            new_pruned_node = typeof(node)(data, node) # create node with no children
            if isdefined(node, :parent) # handling cases where the first decomposition is unuseful
                if direction == "right"
                    node.parent.right = new_pruned_node
                else
                    node.parent.left = new_pruned_node
                end
            else print("No good tree available")
            end
        end
    end
end


"""
    acwptPreOrderBestBasis(node::BinaryNode; et::Wavelets.Entropy)

Returns the best basis tree found using pre order traversal. This is a greedy approach to finding the best basis tree.

# Arguments
-`tree::BinaryNode`: The entry point of the wavelet packet decomposition tree.
-`et::Wavelets.Entropy`: The type of cost function used for evaluation.
"""
function acwptPreOrderBestBasis(tree::BinaryNode; et::Wavelets.Entropy=NormEntropy())
    best_tree = deepcopy(tree)
    acwptPreOrderBestBasis!(best_tree)
    return best_tree
end

acwptPostOrderBestBasis(tree)
acwptPreOrderBestBasis(tree)
