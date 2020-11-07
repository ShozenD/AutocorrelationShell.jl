using AutocorrelationShell, AbstractTrees, Plots, LinearAlgebra

include("test_trees.jl");
## Define base node object
mutable struct BinaryNodePlus{T}
    data::T
    parent::BinaryNodePlus{T}
    left::BinaryNodePlus{T} # pointers to children
    right::BinaryNodePlus{T}

    # for graphing purposes
    ycoord::Int # tree depth
    xcoord::Int # xcoordinate
    offset = 0 # distance to each son
    thread::Bool # whether it needs threading

    # Root constructor
    BinaryNodePlus{T}(data) where T = new{T}(data)
    # Child node constructor
    BinaryNodePlus{T}(data, parent::BinaryNodePlus{T}) where T = new{T}(data, parent)
end

BinaryNodePlus(data) = BinaryNodePlus{typeof(data)}(data)

mutable struct ExtremeNode{T}
    addr::BinaryNodePlus
    off::Int=0
    lev::Int=0

    # Root constructor
    ExtremeNode{T}(node) where T = new{T}(node)
end

ExtremeNode(node) = ExtremeNode{typeof(node)}(node)

ExtremeNode(root.right)

## base node methods
function leftchildplus(data, parent::BinaryNodePlus)
    !isdefined(parent, :left) || error("left child is already assigned")
    node = typeof(parent)(data, parent)
    parent.left = node
end

function rightchildplus(data, parent::BinaryNodePlus)
    !isdefined(parent, :right) || error("right child is already assigned")
    node = typeof(parent)(data, parent)
    parent.right = node
end

root = BinaryNodePlus([1]);
leftchildplus([1], root);
rightchildplus([1], root);

root::BinaryNode # Root of subtree
level::Int # Current overall level
rmost::Int
lmost::Int # Extreme decendants

root.left # left and right sons
root.right
lr, ll, rr, rl # LR = rightmost node on lowest level of left subtree and so on...

cursep # separation on current level
rootsep # current separation at node t
loffsum, roffsum # offset from l & r to T

function tidier_tree(root::BinaryNodePlus, level::Int, minsep::Float64,
    rmost::ExtremeNode, lmost::ExtremeNode)

    if !isdefined(root)
        rmost.lev = -1;
        lmost.lev = -1;
    else
        if isdefined(root.left) & isdefined(root.right) # not a leaf
            root.ycoord = level;
            left = root.left;
            right = root.right;
            tidier_tree(left, level + 1, lr, ll);
            tidier_tree(right, level + 1, rr, rl);

            cursep = minsep;
            rootsep = minsep;
            loffsum = 0;
            roffsum = 0;

            while isdefined(left) & isdefined(right)
                if cursep < minsep # push
                    rootsep = rootsep + (minsep - cursep);
                    minsep = cursep;
                end

                # advance l and r
                if isdefined(left.right)
                    loffsum = loffsum + left.offset;
                    cursep = cursep - left.offset;
                    left = left.right
                else
                    loffsum = loffsum - left.offset;
                    cursep = cursep + left.offset;
                    left = left.left
                end
                if isdefined(right.left)
                    roffsum = roffsum - right.offset;
                    cursep = cursep - right.offset;
                    right = right.left
                else
                    roffsum = roffsum + right.offset;
                    cursep = cursep + right.offset;
                    right = right.right
                end
            end

            # set the offset in node T, and include itt in accumulated offsets for L and R
            root.offset = (rootsep + 1)/2;
            loffsum = loffsum - root.offset;
            roffsum = roffsum + root.offset;

            # Update extreeme descendants information
            if (rl.lev > ll.lev) | !isdefined(root.left)
                lmost = rl
                lmost = lmost.off + root.offset
            else
                lmost = ll;
                lmost.off = lmost.off - root.offset
            end
            if (lr.lev > rr.lev) | !isdefined(root.right)
                rmost = lr
                rmost = rmost.off - root.offset
            else
                rmost = rr;
                rmost.off = rmost.off + root.offset
            end

            # If subrees of root were of uneven heights, check to see if threading is necessary.
            # At most one thred needs to be inserted.

            if isdefined(left) & (left != root.left)
                rr.addr.thread = true;
                rr.addr.offset = abs((rr.off + root.offset) - loffsum);
                if loffsum - t.offset <= rr.off
                    rr.addr.left = left
                else
                    rr.addr.right = left
                end
            elseif isdefined(right) & (right != root.right)
                ll.addr.thread = true;
                ll.addr.offest = abs((ll.off - root.offset) - roffsum);
                if roffsum + root.offset >= ll.off
                    ll.addr.right = right
                else
                    ll.addr.left = right
                end
            end

        else # leaf
            rmost.addr = root;
            lmost.addr = root;
            rmost.lev = level;
            lmost.lev = level;
            rmost.off = 0;
            lmost.off = 0;
            root.offset = 0;
        end
    end
end
