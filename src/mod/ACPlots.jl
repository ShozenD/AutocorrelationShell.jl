module ACPlots
export
    selectednodes_plot,
    wiggle,
    wiggle!
using ..ACWT, ..ACTransforms, ..ACUtil
using AbstractTrees, LinearAlgebra, Plots, Wavelets

"""
    nodes_to_bitmap(x::AcwptNode)

Translates the binary tree datastructure to a bitmap for visualization.
"""
function treetobitmap(x::AcwptNode)
    nrow, ncol = length(x.data), maxtransformlevels(x.data) + 1
    arr, hash = zeros(nrow, ncol), zeros(ncol - 1)

    @inbounds begin
        for n in collect(PreOrderDFS(x))
            d = n.depth
            l = length(n.data)/2^d
            _start = Int(hash[d]+1)
            _end = Int(_start + l - 1)
            if !isdefined(n, :left)
                hash[d+1:end] .= _start + l/2 - 1
                arr[_start:_end, d+1] .= 1
            end
            if !isdefined(n, :right)
                hash[d+1:end] .= _end
                arr[_start:_end, d+1] .= 1
            end
            hash[d] += l
        end
    end

    return arr
end

function treetobitmap(x::BitArray)
    M = length(x)
    ncol::Int = log2(length(x)+1) # depth of tree
    nrow::Int = 2^(ncol-1) # length of original signal
    arr = zeros(nrow, ncol)

    for (idx, flg) in enumerate(x)
        if x[idx]==1 && rightchild(idx) > M
            arr[idx,ncol] = 1
        elseif x[idx]==1 && (x[leftchild(idx)]==0 && x[rightchild(idx)]==0) #leafnode
            d = floor(Int,log2(idx)) # depth
            l = nrow >> d # length of subband
            _start = (idx-2^d)*l+1
            _end = _start+l-1
            arr[_start:_end,d+1] .= 1
        end
    end

    return arr
end

"""
    selectednodes_plot(x::BitArray, nodecolor::Symbol = :red)

Given a best basis binary tree, outputs a visual representation of the selected nodes.
"""
function selectednodes_plot(x::BitArray, nodecolor::Symbol=:red)
    bitmap = treetobitmap(x)
    nrow, ncol = size(bitmap)
    p = heatmap(
        transpose(bitmap),
        color = [:black, nodecolor],
        yflip=true,
        legend=false,
        xlims = (1, nrow+0.5),
        ylims = (0.5, ncol+0.5),
        xticks = false,
        yticks = 0:ncol
    )
    hline!([1.5:1:(ncol-0.5);], color=:white)
    @inbounds begin
        for i in 1:ncol
            for j in 1:2^(i-1)
                vpos = (nrow/2^i)*(2*j-1) + 0.5
                plot!(vpos*ones(ncol-i+1), (i+0.5):(ncol+0.5), color=:white)
            end
        end
    end
    return p
end

function selectednodes_plot(x::AcwptNode, nodecolor::Symbol=:red)
    bitmap = treetobitmap(x)
    nrow, ncol = size(bitmap)
    p = heatmap(
        transpose(bitmap),
        color = [:black, nodecolor],
        yflip=true,
        legend=false,
        xlims = (1, nrow+0.5),
        ylims = (0.5, ncol+0.5),
        xticks = false,
        yticks = 0:ncol
    )
    hline!([1.5:1:(ncol-0.5);], color=:white)
    @inbounds begin
        for i in 1:ncol
            for j in 1:2^(i-1)
                vpos = (nrow/2^i)*(2*j-1) + 0.5
                plot!(vpos*ones(ncol-i+1), (i+0.5):(ncol+0.5), color=:white)
            end
        end
    end
    return p
end

## Visualizations
"""
    wiggle(wav; taxis=1:size(wav,1), zaxis=1:size(wav,2), sc=1, EdgeColor=:black, FaceColor=:black, Orient=:across, Overlap=true, ZDir=:normal)

Plot a set of shaded wiggles

# Arguments
- `wav::Array`: matrix of waveform columns.
- `taxis::Array=1:size(wav,1)`: time axis vector
- `zaxis::Array=1:size(wav,2)`: space axis vector
- `sc::Float=1`: scale factor/magnification.
- `EdgeColor::Symbol=:black`: Sets edge of wiggles color.
- `FaceColor::Symbol=:black`: Sets shading color of wiggles.
- `Overlap::bool=true`: How signals are scaled.
        true  - Signals overlap (default);
        false - Signals are scaled so they do not overlap.
- `Orient::Symbol=:across`: Controls orientation of wiggles.
        :across - from left to right
        :down   - from top to down
- `ZDir::Symbol=:normal`: Direction of space axis.
        :normal  - First signal at bottom (default)
        :reverse - First signal at top.

Translated by Nicholas Hausch -- MATLAB file provided by Naoki Saito
The previous MATLAB version contributors are:
 Anthony K. Booer (SLB) and Bradley Marchand (NSWC-PC)
Revised by Naoki Saito, Feb. 05, 2018
"""
function wiggle(wav::AbstractArray{Float64};
                taxis=1:size(wav,1), zaxis=1:size(wav,2),
                sc::Integer=1, EdgeColor::Symbol=:black, FaceColor::Symbol=:black,
                Overlap::Bool=true, Orient::Symbol=:across, ZDir::Symbol=:normal)

    # Set axes
    (n,m) = size(wav)

    # Sanity check
    if length(taxis) != n
        error("Inconsistent taxis dimension!")
    end
    if length(zaxis) != m
        error("Inconsistent zaxis dimension!")
    end

    # For calculation purposes
    maxrow = zeros(m); minrow = zeros(m)
    for k = 1:m
      maxrow[k] = maximum(wav[:,k]); minrow[k] = minimum(wav[:,k])
    end

    # Scale the data for plotting
    wamp = deepcopy(wav)
    dt = mean(diff(taxis))
    dz = mean(diff(zaxis))
    if Overlap
      wamp *= 2 * dz * (sc/maximum(maxrow-minrow))
    else
      wmax = maximum(maxrow); wmin = minimum(minrow);
      if wmax<=0
        wmax = 0
      end
      if wmin>=0
        wmin = 0
      end
        wamp = sc*wav/(wmax-wmin)
    end

    # Set initial plot
    t0 = minimum(taxis)
    t1 = maximum(taxis)
    z0 = minimum(zaxis)
    z1 = maximum(zaxis)
    if Orient == :down
     plot(xlims=(z0-dz,z1+dz), ylims=(t0,t1), yflip=true, legend=:none)
    else
     plot(xlims=(t0,t1), ylims=(z0-dz,z1+dz), legend=:none)
    end
    if ZDir == :reverse
        wamp = flipdim(wamp,2)
    end

    # Plot each wavelet
    for k = 1:m
      sig = wamp[:,k]
      t = deepcopy(taxis)
      w_sign = sign.(sig)
      for j=1:n-1
        if (w_sign[j]!=w_sign[j+1] && w_sign[j]!=0 && w_sign[j+1]!=0)
          sig = [sig; 0]
          t = [t; t[j]-sig[j]*(t[j+1]-t[j])/(sig[j+1]-sig[j])]
        end
      end
      IX = sortperm(t)
      t = t[IX]
      sig = sig[IX]
      len = length(t)
      len1 = collect(len:-1:1)
      indperm = [1:len;len1]
      inputx = t[indperm]
      inputy = zaxis[k] .+ [sig;min.(sig[len1],0)]
        # In the plot! functions below, theoretically speaking, either
        # fillrange = zaxis[k] or fillrange=[zaxis[k], zaxis[k]+dz] should be used.
        # However, those do not generate the desired plots as of O3/19/2018.
        # Somehow, the relative value of 0, i.e., fillrange=0, works well,
        # which is used temporarily.
      if Orient == :down
        plot!(inputy, inputx, fillrange=0, fillalpha=0.75, fillcolor=FaceColor, linecolor=EdgeColor, orientation=:v)
      else
        plot!(inputx, inputy, fillrange=0, fillalpha=0.75, fillcolor=FaceColor, linecolor=EdgeColor)
      end
    end
    plot!() # flushing the display.
end

"""
    wiggle!(wav; taxis=1:size(wav,1), zaxis=1:size(wav,2), sc=1, EdgeColor=:black, FaceColor=:black, Orient=:across, Overlap=true, ZDir=:normal)

Plot a set of shaded wiggles on the current displayed graphics

# Arguments
- `wav::Array`: matrix of waveform columns.
- `taxis::Array=1:size(wav,1)`: time axis vector
- `zaxis::Array=1:size(wav,2)`: space axis vector
- `sc::Float=1`: scale factor/magnification.
- `EdgeColor::Symbol=:black`: Sets edge of wiggles color.
- `FaceColor::Symbol=:black`: Sets shading color of wiggles.
- `Overlap::bool=true`: How signals are scaled.
        true  - Signals overlap (default);
        false - Signals are scaled so they do not overlap.
- `Orient::Symbol=:across`: Controls orientation of wiggles.
        :across - from left to right
        :down   - from top to down
- `ZDir::Symbol=:normal`: Direction of space axis.
        :normal  - First signal at bottom (default)
        :reverse - First signal at top.

Translated by Nicholas Hausch -- MATLAB file provided by Naoki Saito
The previous MATLAB version contributors are:
 Anthony K. Booer (SLB) and Bradley Marchand (NSWC-PC)
Revised by Naoki Saito, Feb. 05, 2018
"""
function wiggle!(wav; taxis=1:size(wav,1), zaxis=1:size(wav,2), sc=1, EdgeColor=:black, FaceColor=:black, Overlap=true, Orient=:across, ZDir=:normal)

    # Set axes
    (n,m) = size(wav)

    # Sanity check
    if length(taxis) != n
        error("Inconsistent taxis dimension!")
    end
    if length(zaxis) != m
        error("Inconsistent zaxis dimension!")
    end

    # For calculation purposes
    maxrow = zeros(m); minrow = zeros(m)
    for k = 1:m
      maxrow[k] = maximum(wav[:,k]); minrow[k] = minimum(wav[:,k])
    end

    # Scale the data for plotting
    wamp = deepcopy(wav)
    dt = mean(diff(taxis))
    dz = mean(diff(zaxis))
        if Overlap
      wamp *= 2 * dz * (sc/maximum(maxrow-minrow))
    else
      wmax = maximum(maxrow); wmin = minimum(minrow);
      if wmax<=0
        wmax = 0
      end
      if wmin>=0
        wmin = 0
      end
        wamp = sc*wav/(wmax-wmin)
    end

    # Set initial plot
    t0 = minimum(taxis)
    t1 = maximum(taxis)
    z0 = minimum(zaxis)
    z1 = maximum(zaxis)
    if Orient == :down
     plot!(xlims=(z0-dz,z1+dz), ylims=(t0,t1), yflip=true, legend=:none)
    else
     plot!(xlims=(t0,t1), ylims=(z0-dz,z1+dz), legend=:none)
    end
    if ZDir == :reverse
        wamp = flipdim(wamp,2)
    end

    # Plot each wavelet
    for k = 1:m
      sig = wamp[:,k]
      t = deepcopy(taxis)
      w_sign = sign.(sig)
      for j=1:n-1
        if (w_sign[j]!=w_sign[j+1] && w_sign[j]!=0 && w_sign[j+1]!=0)
          sig = [sig; 0]
          t = [t; t[j]-sig[j]*(t[j+1]-t[j])/(sig[j+1]-sig[j])]
        end
      end
      IX = sortperm(t)
      t = t[IX]
      sig = sig[IX]
      len = length(t)
      len1 = collect(len:-1:1)
      indperm = [1:len;len1]
      inputx = t[indperm]
      inputy = zaxis[k] .+ [sig;min.(sig[len1],0)]
        # In the plot! functions below, theoretically speaking, either
        # fillrange = zaxis[k] or fillrange=[zaxis[k], zaxis[k]+dz] should be used.
        # However, those do not generate the desired plots as of O3/19/2018.
        # Somehow, the relative value of 0, i.e., fillrange=0, works well,
        # which is used temporarily.
      if Orient == :down
        plot!(inputy, inputx, fillrange=0, fillalpha=0.75, fillcolor=FaceColor, linecolor=EdgeColor, orientation=:v)
      else
        plot!(inputx, inputy, fillrange=0, fillalpha=0.75, fillcolor=FaceColor, linecolor=EdgeColor)
      end
    end
    plot!() # flushing the display.
end

end # End module
