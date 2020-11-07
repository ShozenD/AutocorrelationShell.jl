using AutocorrelationShell, AbstractTrees, Plots, Wavelets
include("test_trees.jl")

function collectchildren(node) return collect(PostOrderDFS(node)) end

function first_walk(tree)
    for node in collectchildren(tree)
        if node.ycoord > 0
            if node.parent.left == node
                node.xcoord = 0
            else
                node.xcoord = 1
            end
        else
            node.xcoord = 0
        end

        position_parent(node)
    end
end

function position_parent(node)
    if isdefined(node, :left)
        node.xcoord = (node.left.xcoord + node.right.xcoord)/2
    end
end

first_walk(tree1)

function contourRight(node)
    if isdefined(node, :right)
        return node.offset + contourRight(node.right)
    else
        return node.xcoord + node.offset
    end
end

function contourLeft(node)
    if isdefined(node, :left)
        return node.offset + contourLeft(node.left)
    else
        return node.xcoord + node.offset
    end
end

function avoid_conflict(node)
    diff = contourRight(node.left) - contourLeft(node.right)
    if diff > 0
        node.right.offset += diff + 1
    end
end

avoid_conflict(tree1)

for node in collect(PreOrderDFS(tree1))
    node.xcoord += node.offset

    if isdefined(node, :right)
        node.right.offset += node.offset
        node.left.offset += node.offset
    end
end

tree1.xcoord = (tree1.right.xcoord + tree1.left.xcoord)/2

tree1
