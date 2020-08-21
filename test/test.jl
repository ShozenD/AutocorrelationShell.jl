include("../src/AutocorrelationShell.jl")
using Main.AutocorrelationShell
using Random, Wavelets, AbstractTrees, LinearAlgebra
rng = MersenneTwister(123);

X₁ = randn(rng, 4); # length 4 random signal
H = wavelet(WT.db2);
Q = qfilter(H);
P = pfilter(H);
decomp = acwpt(X₁, P, Q)

# Print the tree in the console
print_tree(decomp)

# Gather all nodes into a vector
collect(PostOrderDFS(decomp))

## Best Basis Tree Algorithm
function acwptBestBasisTree(node::BinaryNode; direction::AbstractString="right", et::Wavelets.Entropy=NormEntropy())
    if isdefined(node, :left) & isdefined(node, :right)
        e0 = wentropy(node.data, et)
        e10 = wentropy(node.left.data, et)
        e11 = wentropy(node.right.data, et)
        if (e10 + e11)/2 < e0
            acwptBestBasisTree(node.left, direction="left")
            acwptBestBasisTree(node.right)
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

acwptBestBasisTree(decomp, et=ShannonEntropy())
print_tree(decomp)

## Reconstruction Algorithm
function aciwpt(tree)
    if isdefined(tree, :left)
        left = aciwpt(tree.left)
    end

    if isdefined(tree, :right)
        right = aciwpt(tree.right)
    end

    if !isdefined(tree, :left) & !isdefined(tree, :right)
        return tree.data
    end

    return (left + right)/sqrt(2)
end

reconst = aciwpt(decomp)

norm(X₁ - reconst)

## Experiment
decomp = acwpt(X₁, P, Q)

decomp.left.data = zeros(4);
decomp.left.left.data = zeros(4);
decomp.left.right.data = zeros(4);

decomp.right.data = zeros(4);
decomp.right.right.data = zeros(4);
decomp.right.left.data = [0, 1, 0, 0];

reconst = aciwpt(decomp)
decomp2 = acwpt(reconst, P, Q);

acwptBestBasisTree(decomp2);

print_tree(decomp2)
