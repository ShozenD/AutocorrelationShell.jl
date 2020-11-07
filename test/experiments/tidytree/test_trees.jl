using AbstractTrees, Plots, Wavelets

## Base node object
mutable struct BinaryNodePlus{T}
    data::T
    parent::BinaryNodePlus{T}
    left::BinaryNodePlus{T} # pointers to children
    right::BinaryNodePlus{T}

    # for graphing purposes
    ycoord::Int # tree depth
    xcoord::Float64  # xcoordinate
    offset::Int # distance to each son
    thread::Bool # whether it needs threading

    # Root constructor
    BinaryNodePlus{T}(data) where T = new{T}(data)
    # Child node constructor
    BinaryNodePlus{T}(data, parent::BinaryNodePlus{T}) where T = new{T}(data, parent)
end

BinaryNodePlus(data) = BinaryNodePlus{typeof(data)}(data)

function initBinaryNodePlus(data)
    node = BinaryNodePlus(data)
    node.ycoord = 0
    node.xcoord = 0
    node.offset = 0
    node.thread = false
    return node
end

function initBinaryNodePlus(data, parent::BinaryNodePlus)
    node = typeof(parent)(data, parent)
    node.ycoord = 0
    node.xcoord = 0
    node.offset = 0
    node.thread = false
    return node
end

function leftchildplus(data, parent::BinaryNodePlus)
    !isdefined(parent, :left) || error("left child is already assigned")
    node = initBinaryNodePlus(data, parent)
    parent.left = node
end

function rightchildplus(data, parent::BinaryNodePlus)
    !isdefined(parent, :right) || error("right child is already assigned")
    node = initBinaryNodePlus(data, parent)
    parent.right = node
end

function AbstractTrees.children(node::BinaryNodePlus)
    if isdefined(node, :left)
        if isdefined(node, :right)
            return (node.left, node.right)
        end
        return (node.left,)
    end
    isdefined(node, :right) && return (node.right,)
    return ()
end

## Initialize Trees

tree1 = initBinaryNodePlus([1])
leftchildplus([1], tree1);
rightchildplus([1], tree1);

leftchildplus([1], tree1.left);
rightchildplus([1], tree1.left);

leftchildplus([1], tree1.right);
rightchildplus([1], tree1.right);

tree1.ycoord = 0;
tree1.left.ycoord = 1;
tree1.right.ycoord = 1;

tree1.left.left.ycoord = 2;
tree1.left.right.ycoord = 2;
tree1.right.left.ycoord = 2;
tree1.right.right.ycoord = 2;
