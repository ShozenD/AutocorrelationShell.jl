using Wavelets, AutocorrelationShell, Plots, LinearAlgebra, AbstractTrees

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

function fwt_ac_temp(x::Vector{T}, L::Integer, P::Vector{T}, Q::Vector{T}) where T<:Number

	n = length(x)
	J = dyadlength(x)

    # Sanity Check
    @assert L >= 0
    @assert L <= J

	D = J-L
	wp = zeros(n,D+1)

	wp[:,1] = x
    @inbounds begin
    	for d=0:(D-1)
            println(d)
    		for b=0:(2^d-1)
    		   s = wp[echant(n,d,b),1]
    		   h = ac_filter(s,Q)
    		   l = ac_filter(s,P)
    		   wp[echant(n,d,b),D+1-d] = h
    		   wp[echant(n,d,b),1] = l
    		 end
    	end
    end
	return wp
end


## Wavelet Packets functions
function acwpt_temp(x::Vector{T}, node::BinaryNode,
    P::Vector{T}, Q::Vector{T}, d::Integer) where T<:Real

    n = length(x)
    J = dyadlength(x) # maximum tree depth
    left, right = zeros(n), zeros(n)

    if d < J
        for b = 0:(2^d-1) # depth starts from 2
            s = x[echant(n, d, b)]
            h = ac_filter(s, Q)
            l = ac_filter(s, P)
            left[echant(n, d, b)] = l
            right[echant(n, d, b)] = h
        end

        # left
        leftchild(left, node)
        acwpt_temp(left, node.left, P, Q, d+1)

        # right
        rightchild(right, node)
        acwpt_temp(right, node.right, P, Q, d+1)
    end
end

function acwpt_temp(x::Vector{T}, P::Vector{T}, Q::Vector{T}) where T<:Real
    root = BinaryNode(x) # original signal
    acwpt_temp(x, root, P, Q, 0)

    return root
end


## base node methods
Q = qfilter(wavelet(WT.db2));
P = pfilter(wavelet(WT.db2));

x = zeros(256); x[128] = 1; # One hot signal
decomp = acwt(x, L=2, P=P, Q=Q)

X = decomp[:,4];

X = acwt(x, L=2, P=P, Q=Q)[:,4];
decomp = acwt(X, L=0, P=P, Q=Q);
tree = acwpt_temp(X, P, Q);

function acwptPostOrderBestBasis_temp!(node::BinaryNode; direction::AbstractString="right", et::Wavelets.Entropy=NormEntropy())
    if isdefined(node, :left) & isdefined(node, :right)
        el = acwptPostOrderBestBasis_temp!(node.left, direction="left")
        er = acwptPostOrderBestBasis_temp!(node.right)
        if (el + er)/2 < wentropy(node.data, et)
            return (el + er)/2
        elseif isdefined(node, :parent) # handling cases where the first decomposition is unuseful
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
    acwptPostOrderBestBasis(tree::BinaryNode; et::Wavelets.Entropy=NormEntropy())

Returns the best basis tree found using post order traversal. This is a democratic approach to finding the best basis tree.

# Arguments
-`tree::BinaryNode`: The entry point of the wavelet packet decomposition tree.
-`et::Wavelets.Entropy`: The type of cost function used for evaluation.
"""
function acwptPostOrderBestBasis_temp(tree::BinaryNode; et::Wavelets.Entropy=NormEntropy())
    best_tree = deepcopy(tree)
    acwptPostOrderBestBasis_temp!(best_tree, et=et)
    return best_tree
end

acwptPostOrderBestBasis_temp(tree)
norm(collect(PostOrderDFS(tree))[1].data - decomp[:,1]) < 1e-15

norm(fwt_ac_temp(x, 0, P, Q)[:,1] - collect(PostOrderDFS(acwpt_temp(x, P, Q)))[1].data) == 0
