"""
    wiggle(wav; taxis=1:size(wav,1), zaxis=1:size(wav,2), sc=1, EdgeColor=:black, FaceColor=:black, Orient=:across, Overlap=true, ZDir=:normal)

Plot a set of shaded wiggles
# Arguments
* 'wav::Array': matrix of waveform columns.
* 'taxis::Array=1:size(wav,1)': time axis vector
* 'zaxis::Array=1:size(wav,2)': space axis vector
* 'sc::Float=1': scale factor/magnification.
* 'EdgeColor::Symbol=:black': Sets edge of wiggles color.
* 'FaceColor::Symbol=:black': Sets shading color of wiggles.
* 'Overlap::bool=true': How signals are scaled.
        true  - Signals overlap (default);
        false - Signals are scaled so they do not overlap.
* 'Orient::Symbol=:across': Controls orientation of wiggles.
        :across - from left to right
        :down   - from top to down
* 'ZDir::Symbol=:normal': Direction of space axis.
        :normal  - First signal at bottom (default)
        :reverse - First signal at top.

Translated by Nicholas Hausch -- MATLAB file provided by Naoki Saito
The previous MATLAB version contributors are:
 Anthony K. Booer (SLB) and Bradley Marchand (NSWC-PC)
Revised by Naoki Saito, Feb. 05, 2018
"""
function wiggle(wav; taxis=1:size(wav,1), zaxis=1:size(wav,2), sc=1, EdgeColor=:black, FaceColor=:black, Overlap=true, Orient=:across, ZDir=:normal)

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
* 'wav::Array': matrix of waveform columns.
* 'taxis::Array=1:size(wav,1)': time axis vector
* 'zaxis::Array=1:size(wav,2)': space axis vector
* 'sc::Float=1': scale factor/magnification.
* 'EdgeColor::Symbol=:black': Sets edge of wiggles color.
* 'FaceColor::Symbol=:black': Sets shading color of wiggles.
* 'Overlap::bool=true': How signals are scaled.
        true  - Signals overlap (default);
        false - Signals are scaled so they do not overlap.
* 'Orient::Symbol=:across': Controls orientation of wiggles.
        :across - from left to right
        :down   - from top to down
* 'ZDir::Symbol=:normal': Direction of space axis.
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
