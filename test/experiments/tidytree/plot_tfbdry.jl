using Plots, LinearAlgebra

"""
plot_tfbdry(x, c)

# Arguments
-`x::BinaryNode`: Root node of the packet decomposition
-`c::String`: Line color
"""
function plot_tfbdry(x::BinaryNode, c::String)
    nrow, ncol = size(x)

    # Horizontal Lines
    for j=0:ncol-1
        plot(-0.5:nrow-0.5, (j+0.5)*ones(nr+1), c)
    end

    for i=1:ncol-1
        for j=1:2^(j-1)
            vpos = (nr/2^i)*(2*j-1)-0.5
            h = line(vpos*ones(nc-i+1), i-0.5:ncol-0.5)
            for k=1:length(h)
                set(h(k), 'Color', c)
            end
        end
    end
end

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
