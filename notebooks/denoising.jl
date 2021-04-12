### A Pluto.jl notebook ###
# v0.14.1

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ f9e001ae-0323-4283-9af1-1a8252b503e7
let
    import Pkg
    Pkg.activate(".")
end

# ╔═╡ 203fb8c8-4358-4908-b616-a691ce329c02
using Wavelets, AutocorrelationShell, LinearAlgebra, Plots, DataFrames, CSV, PlutoUI

# ╔═╡ 53c257e0-96ba-11eb-3615-8bfed63b2c18
md"# Denoising Experiment"

# ╔═╡ bc641876-9723-4c68-8221-ad03fb695c82
md"Let's load some test functions"

# ╔═╡ d0c98a14-ad3e-4a0b-b889-a3ea86f888f3
testdata = CSV.read("../test/data/wavelet_test_256.csv", DataFrame)

# ╔═╡ 371a96c3-e6d6-416b-8079-fe0da34d7dc8
@bind signal_name Select(
	["blocks", "bumps", "heavy_sine", "doppler", "quadchirp", "mishmash"]
)

# ╔═╡ 747ce479-6095-4d16-a2f7-28b33510ad24
x = testdata[!,signal_name]

# ╔═╡ a184ae65-7947-4ffb-b751-8b6e97a8608b
function addnoise(x::AbstractArray{<:Number,1}, s::Real=0.1)
	ϵ = randn(length(x))
	ϵ = ϵ/norm(ϵ)
	y = x + ϵ * s * norm(x)
	return y
end

# ╔═╡ 0e0fd333-471e-416b-a7b3-05b094debea8
@bind noise_size Slider(0:0.01:1)

# ╔═╡ 548bc51b-0242-4b2f-8ba3-104d0dd453b6
md"The energy of noise proportional to signal: $noise_size"

# ╔═╡ 851b04bb-e382-4a0a-98f6-7d4b983ca5ab
begin
	p1 = plot(x, ylim = (minimum(x)-1,maximum(x)+1), label = "Original signal")
	x_noisy = addnoise(x, noise_size)
	plot!(x_noisy, label = "Noisy signal");
	
	y = acwt(x_noisy, wavelet(WT.db4), 4);
	p2 = wiggle(y)
	
	plot(p1, p2, layout = (2,1))
end

# ╔═╡ 341c2551-b625-47a0-9163-3c0c0e7d4e13
histogram(vec(abs.(y)), ylim = (0,1000), legend = false)

# ╔═╡ 10452433-7123-4006-9e3d-ae4245fefdc5
threshold!(y, HardTH(), 0.2);

# ╔═╡ fd0b11ee-a0b8-4c0d-8c88-c69996b1c42d
r = iacwt(y)

# ╔═╡ 9b4ef541-9a36-4bc0-8654-10ab0a4e63b3
begin
	plot(x, label = "original")
	plot!(r, label = "denoised")
end

# ╔═╡ 18144929-3f31-42b2-9e27-df146a687ae0
norm(x - r)/length(x)

# ╔═╡ 11e63c9a-6124-4122-9a86-ceed926d25d2
md"# Data Setup
* Generate a set of original signals $\rightarrow$ `x₀`
* Add noise to original signals $\rightarrow$ `x`
"

# ╔═╡ 24bf8f8d-3880-4635-a610-93df82132a3e


# ╔═╡ 126c41e7-dd65-46c6-8c5b-2439f5624fd5
md"# Non-Redundant Transforms"

# ╔═╡ 61d745d8-5c74-479b-9698-cd50bb68b3c7


# ╔═╡ f2f949f8-772f-4787-8883-0d96137f0924
md"# Redundant Transforms"

# ╔═╡ 3895472f-0a4f-4b7a-84f6-470208b5e8cc
md"## Autocorrelation Wavelet Transforms"

# ╔═╡ 89a56a57-a5b9-4380-a618-97d8b901c01b


# ╔═╡ 0440ce41-8c23-45e6-aa67-56c6a58298a6


# ╔═╡ cb927f51-99f2-4eeb-a0e5-5f1c65464b6f
md"## Stationary Wavelet Transforms"

# ╔═╡ ab9b089f-fbb7-4436-8d6a-963db1e95670


# ╔═╡ 8774496e-a184-4c9a-9335-d2f184673cf5


# ╔═╡ Cell order:
# ╟─53c257e0-96ba-11eb-3615-8bfed63b2c18
# ╠═f9e001ae-0323-4283-9af1-1a8252b503e7
# ╠═203fb8c8-4358-4908-b616-a691ce329c02
# ╟─bc641876-9723-4c68-8221-ad03fb695c82
# ╠═d0c98a14-ad3e-4a0b-b889-a3ea86f888f3
# ╟─371a96c3-e6d6-416b-8079-fe0da34d7dc8
# ╟─747ce479-6095-4d16-a2f7-28b33510ad24
# ╠═a184ae65-7947-4ffb-b751-8b6e97a8608b
# ╟─0e0fd333-471e-416b-a7b3-05b094debea8
# ╟─548bc51b-0242-4b2f-8ba3-104d0dd453b6
# ╠═851b04bb-e382-4a0a-98f6-7d4b983ca5ab
# ╠═341c2551-b625-47a0-9163-3c0c0e7d4e13
# ╠═10452433-7123-4006-9e3d-ae4245fefdc5
# ╠═fd0b11ee-a0b8-4c0d-8c88-c69996b1c42d
# ╠═9b4ef541-9a36-4bc0-8654-10ab0a4e63b3
# ╠═18144929-3f31-42b2-9e27-df146a687ae0
# ╟─11e63c9a-6124-4122-9a86-ceed926d25d2
# ╠═24bf8f8d-3880-4635-a610-93df82132a3e
# ╟─126c41e7-dd65-46c6-8c5b-2439f5624fd5
# ╠═61d745d8-5c74-479b-9698-cd50bb68b3c7
# ╟─f2f949f8-772f-4787-8883-0d96137f0924
# ╟─3895472f-0a4f-4b7a-84f6-470208b5e8cc
# ╠═89a56a57-a5b9-4380-a618-97d8b901c01b
# ╠═0440ce41-8c23-45e6-aa67-56c6a58298a6
# ╟─cb927f51-99f2-4eeb-a0e5-5f1c65464b6f
# ╠═ab9b089f-fbb7-4436-8d6a-963db1e95670
# ╠═8774496e-a184-4c9a-9335-d2f184673cf5
