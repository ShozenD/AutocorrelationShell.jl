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
using Wavelets, AutocorrelationShell, LinearAlgebra, Plots, DataFrames, CSV, PlutoUI, WaveletsExt, Statistics, Gadfly

# ╔═╡ 53c257e0-96ba-11eb-3615-8bfed63b2c18
md"# Denoising Experiment"

# ╔═╡ ce4bf94e-edef-40c2-8ac5-b741e47a1759
md"### Activate environment"

# ╔═╡ 85c26ead-f043-42c8-8245-58c8d03a963d
md"### Import libraries"

# ╔═╡ 46e30602-a850-4175-b4bf-b4ef4b5359aa
md"## Exploratory Data Analysis"

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
@bind noise_size Slider(0:0.01:1, default=0.5, show_value=true)

# ╔═╡ 548bc51b-0242-4b2f-8ba3-104d0dd453b6
md"The energy of noise proportional to signal: $noise_size"

# ╔═╡ 851b04bb-e382-4a0a-98f6-7d4b983ca5ab
begin
	p1 = Plots.plot(x, ylim = (minimum(x)-1,maximum(x)+1), label = "Original signal")
	x_noisy = addnoise(x, noise_size)
	Plots.plot!(x_noisy, label = "Noisy signal");
	
	y = acwt(x_noisy, wavelet(WT.db4), 4);
	p2 = WaveletsExt.wiggle(y)
	
	Plots.plot(p1, p2, layout = (2,1))
end

# ╔═╡ 341c2551-b625-47a0-9163-3c0c0e7d4e13
Plots.histogram(vec(abs.(y)), ylim = (0,1000), legend = false)

# ╔═╡ 10452433-7123-4006-9e3d-ae4245fefdc5
threshold!(y, HardTH(), 0.2);

# ╔═╡ fd0b11ee-a0b8-4c0d-8c88-c69996b1c42d
r = iacwt(y)

# ╔═╡ 9b4ef541-9a36-4bc0-8654-10ab0a4e63b3
begin
	Plots.plot(x, label = "original")
	Plots.plot!(r, label = "denoised")
end

# ╔═╡ 18144929-3f31-42b2-9e27-df146a687ae0
norm(x - r)/length(x)

# ╔═╡ 11e63c9a-6124-4122-9a86-ceed926d25d2
md"# Data Setup
* Generate a set of original signals $\rightarrow$ `X₀`
* Add noise to original signals $\rightarrow$ `X`
"

# ╔═╡ 708c818b-9b35-40e2-b7f0-dd9d5ebdfa85
begin
	X₀ = generatesignals(x, 100, 2)
	X = hcat([addnoise(X₀[:,i], noise_size) for i in axes(X₀,2)]...)
end

# ╔═╡ d95b8d0a-47dc-46ac-9ee8-dd8647dd6c7f
@bind wavelet_type Select(
	["WT.haar", 
	"WT.db1", "WT.db2", "WT.db3", "WT.db4", "WT.db5", 
	"WT.db6", "WT.db7", "WT.db8", "WT.db9", "WT.db10",
	"WT.coif2", "WT.coif4", "WT.coif6", "WT.coif8",
	"WT.sym4", "WT.sym5", "WT.sym6", "WT.sym7", "WT.sym8", "WT.sym9", "WT.sym10",
	"WT.batt2", "WT.batt4", "WT.batt6"]
)

# ╔═╡ f9a7488d-12a6-45f0-9c70-e67448dfe637
@bind threshold_method Select(
	["HardTH()" => "Hard",
	"SoftTH()" => "Soft"]	
)

# ╔═╡ 6b02c425-39b9-467f-9406-3e9096873af4
begin
	wt = wavelet(eval(Meta.parse(wavelet_type)))
	th = eval(Meta.parse(threshold_method))
	dnt = VisuShrink(th, 256)
	# define variables to store results
	Y = Dict{String, AbstractArray}()
	X̂ = Dict{String, AbstractArray}()
	time = Dict{String, AbstractFloat}()
	mean_psnr = Dict{String, AbstractFloat}()
	results = DataFrame(
		transform = ["None"], 
		threshold = ["None"],
		time = 0.0,
		PSNR = mean([psnr(X[:,i], X₀[:,i]) for i in axes(X,2)])
	)
end

# ╔═╡ 126c41e7-dd65-46c6-8c5b-2439f5624fd5
md"# Non-Redundant Transforms"

# ╔═╡ 17bdc97a-4a0b-4931-a5c6-866f0c814601
md"### Discrete Wavelet Transform"

# ╔═╡ 01e43234-2194-451d-9010-176aa4799fdb
begin 
	Y["DWT"] = cat([wpt(X[:,i], wt) for i in axes(X,2)]..., dims=2)
	X̂["DWT"], time["DWT"] = @timed denoiseall(
		Y["DWT"], :dwt, wt, L=8, dnt=dnt
	)
	mean_psnr["DWT"] = mean([psnr(X̂["DWT"][:,i], X₀[:,i]) for i in axes(X,2)])
	push!(results, ["DWT", threshold_method, time["DWT"], mean_psnr["DWT"]])
end

# ╔═╡ ae8059bd-5b5b-4ff2-a6f0-5ce672bdd54d
md"### Wavelet Packet Transform - Level 4"

# ╔═╡ cf55c5cb-ead6-40b6-896a-8f7e01613a46
begin 
	Y["WPT-L4"] = cat([wpt(X[:,i], wt, 4) for i in axes(X,2)]..., dims=2)
	X̂["WPT-L4"], time["WPT-L4"] = @timed denoiseall(
		Y["WPT-L4"], :wpt, wt, tree=maketree(256,4,:full), dnt=dnt
	)
	mean_psnr["WPT-L4"] = mean([psnr(X̂["WPT-L4"][:,i], X₀[:,i]) for i in axes(X,2)])
	push!(results, ["WPT-L4", threshold_method, time["WPT-L4"], mean_psnr["WPT-L4"]])
end

# ╔═╡ c95ebbed-3d9a-4be2-943b-08c86923ad89
md"### Wavelet Packet Transform - Level 8"

# ╔═╡ 61d745d8-5c74-479b-9698-cd50bb68b3c7
begin 
	Y["WPT-L8"] = cat([wpt(X[:,i], wt) for i in axes(X,2)]..., dims=2)
	X̂["WPT-L8"], time["WPT-L8"] = @timed denoiseall(
		Y["WPT-L8"], :wpt, wt, tree=maketree(256,8,:full), dnt=dnt
	)
	mean_psnr["WPT-L8"] = mean([psnr(X̂["WPT-L8"][:,i], X₀[:,i]) for i in axes(X,2)])
	push!(results, ["WPT-L8", threshold_method, time["WPT-L8"], mean_psnr["WPT-L8"]])
end

# ╔═╡ e5738623-4ea8-4866-a1da-7e849960f4e0
md"### Wavelet Packet Transform - Best Basis"

# ╔═╡ c34c1f0c-7cfd-404b-aa59-d4bb6aa9628f
begin 
	Y["WPD"] = cat([wpd(X[:,i], wt) for i in axes(X,2)]..., dims=3)
	wpt_bt = bestbasistree(Y["WPD"], BB())
	Y["WPT-BT"] = cat([wpt(X[:,i], wt, wpt_bt[:,i]) for i in axes(X,2)]..., dims=2)
	X̂["WPT-BT"], time["WPT-BT"] = @timed hcat(
		[denoise(
			Y["WPT-BT"][:,i], :wpt, wt, tree=wpt_bt[:,i], dnt=dnt) for i in axes(X,2)
		]...
	)
	mean_psnr["WPT-BT"] = mean([psnr(X̂["WPT-BT"][:,i], X₀[:,i]) for i in axes(X,2)])
	push!(results, ["WPT-BT", threshold_method, time["WPT-BT"], mean_psnr["WPT-BT"]])
end

# ╔═╡ 069497ac-fdae-4ddf-8983-026f7d46f07a
md"### Joint Best Basis"

# ╔═╡ cd132003-1384-41a4-bfb4-91247630a24e
begin 
	jbb_tree = bestbasistree(Y["WPD"], JBB())
	Y["JBB"] = cat([wpt(X[:,i], wt, jbb_tree) for i in axes(X,2)]..., dims=2)
	X̂["JBB"], time["JBB"] = @timed denoiseall(
		Y["JBB"], :wpt, wt, tree=jbb_tree, dnt=dnt
	)
	mean_psnr["JBB"] = mean([psnr(X̂["JBB"][:,i], X₀[:,i]) for i in axes(X,2)])
	push!(results, ["JBB", threshold_method, time["JBB"], mean_psnr["JBB"]])
end

# ╔═╡ b075d6d8-228a-4ce2-8647-e2c6b962ba48
md"### Least Statistically Dependent Basis"

# ╔═╡ a7d6a82c-b143-4e8e-94ee-e8999eefc0f1
begin 
	lsdb_tree = bestbasistree(Y["WPD"], LSDB())
	Y["LSDB"] = cat([wpt(X[:,i], wt, lsdb_tree) for i in axes(X,2)]..., dims=2)
	X̂["LSDB"], time["LSDB"] = @timed denoiseall(
		Y["LSDB"], :wpt, wt, tree=lsdb_tree, dnt=dnt
	)
	mean_psnr["LSDB"] = mean([psnr(X̂["LSDB"][:,i], X₀[:,i]) for i in axes(X,2)])
	push!(results, ["LSDB", threshold_method, time["LSDB"], mean_psnr["LSDB"]])
end

# ╔═╡ f2f949f8-772f-4787-8883-0d96137f0924
md"# Redundant Transforms"

# ╔═╡ 3895472f-0a4f-4b7a-84f6-470208b5e8cc
md"## Autocorrelation Wavelet Transforms"

# ╔═╡ 89a56a57-a5b9-4380-a618-97d8b901c01b


# ╔═╡ 0440ce41-8c23-45e6-aa67-56c6a58298a6


# ╔═╡ cb927f51-99f2-4eeb-a0e5-5f1c65464b6f
md"## Stationary Wavelet Transforms"

# ╔═╡ d7da11fd-8768-4ac7-81af-82f68e794e1a
md"### Stationary Discrete Wavelet Transform"

# ╔═╡ ab9b089f-fbb7-4436-8d6a-963db1e95670
begin 
	Y["SDWT"] = cat([sdwt(X[:,i], wt) for i in axes(X,2)]..., dims=3)
	X̂["SDWT"], time["SDWT"] = @timed denoiseall(
		Y["SDWT"], :sdwt, wt, L=8, dnt=dnt
	)
	mean_psnr["SDWT"] = mean([psnr(X̂["SDWT"][:,i], X₀[:,i]) for i in axes(X,2)])
	push!(results, ["SDWT", threshold_method, time["SDWT"], mean_psnr["SDWT"]])
end

# ╔═╡ 72204688-6346-4ff5-b443-839cbd7074d8
md"### Stationary Wavelet Packet Transform - Level 4"

# ╔═╡ 6d09613a-5676-4b45-bf44-7e2c40bb71c9
begin 
	Y["SWPD"] = cat([swpd(X[:,i], wt) for i in axes(X,2)]..., dims=3)
	X̂["SWPT-L4"], time["SWPT-L4"] = @timed denoiseall(
		Y["SWPD"], :swpd, wt, tree=maketree(256,4,:full), dnt=dnt
	)
	mean_psnr["SWPT-L4"] = mean([psnr(X̂["SWPT-L4"][:,i], X₀[:,i]) for i in axes(X,2)])
	push!(
		results, 
		["SWPT-L4", threshold_method, time["SWPT-L4"], mean_psnr["SWPT-L4"]]
	)
end

# ╔═╡ 25de9867-d737-496c-8e0c-29bcdca898e8
md"### Stationary Wavelet Packet Transform - Level 8"

# ╔═╡ 3525a716-0210-4cd2-b780-3f1c27297e09
begin 
	X̂["SWPT-L8"], time["SWPT-L8"] = @timed denoiseall(
		Y["SWPD"], :swpd, wt, tree=maketree(256,8,:full), dnt=dnt
	)
	mean_psnr["SWPT-L8"] = mean([psnr(X̂["SWPT-L8"][:,i], X₀[:,i]) for i in axes(X,2)])
	push!(
		results, 
		["SWPT-L8", threshold_method, time["SWPT-L8"], mean_psnr["SWPT-L8"]]
	)
end

# ╔═╡ 8774496e-a184-4c9a-9335-d2f184673cf5
md"### Stationary Wavelet Packet Transform - Best Basis"

# ╔═╡ 609937c6-8b9f-4987-8f46-ec55ce05e861
begin
	swpt_bt = bestbasistree(Y["SWPD"], BB(stationary=true))
	X̂["SWPT-BT"], time["SWPT-BT"] = @timed hcat(
		[denoise(
			Y["SWPD"][:,:,i], :swpd, wt, tree=swpt_bt[:,i], dnt=dnt
			) for i in axes(X,2)
		]...
	)
	mean_psnr["SWPT-BT"] = mean([psnr(X̂["SWPT-BT"][:,i], X₀[:,i]) for i in axes(X,2)])
	push!(
		results, 
		["SWPT-BT", threshold_method, time["SWPT-BT"], mean_psnr["SWPT-BT"]]
	)
end

# ╔═╡ aa5fff0c-9182-4568-b657-ca48c33de141
md"### Stationary Joint Best Basis"

# ╔═╡ 7a6c247e-d765-4299-bd93-8d8271ca711f
begin 
	sjbb_tree = bestbasistree(Y["SWPD"], JBB(stationary=true))
	X̂["SJBB"], time["SJBB"] = @timed denoiseall(
		Y["SWPD"], :swpd, wt, tree=sjbb_tree, dnt=dnt
	)
	mean_psnr["SJBB"] = mean([psnr(X̂["SJBB"][:,i], X₀[:,i]) for i in axes(X,2)])
	push!(results, ["SJBB", threshold_method, time["SJBB"], mean_psnr["SJBB"]])
end

# ╔═╡ 7103af68-9ad7-4b81-8c21-c8b6d7c9f5be
md"## Results Analysis"

# ╔═╡ 56c7a1b9-c75f-48d8-a602-c219a2f432af
Gadfly.plot(
	results, 
	y=:transform, 
	x=:PSNR, 
	color=:threshold, 
    Guide.title("PSNR")
)

# ╔═╡ Cell order:
# ╟─53c257e0-96ba-11eb-3615-8bfed63b2c18
# ╟─ce4bf94e-edef-40c2-8ac5-b741e47a1759
# ╠═f9e001ae-0323-4283-9af1-1a8252b503e7
# ╟─85c26ead-f043-42c8-8245-58c8d03a963d
# ╠═203fb8c8-4358-4908-b616-a691ce329c02
# ╟─46e30602-a850-4175-b4bf-b4ef4b5359aa
# ╟─bc641876-9723-4c68-8221-ad03fb695c82
# ╠═d0c98a14-ad3e-4a0b-b889-a3ea86f888f3
# ╠═371a96c3-e6d6-416b-8079-fe0da34d7dc8
# ╟─747ce479-6095-4d16-a2f7-28b33510ad24
# ╠═a184ae65-7947-4ffb-b751-8b6e97a8608b
# ╠═0e0fd333-471e-416b-a7b3-05b094debea8
# ╟─548bc51b-0242-4b2f-8ba3-104d0dd453b6
# ╠═851b04bb-e382-4a0a-98f6-7d4b983ca5ab
# ╠═341c2551-b625-47a0-9163-3c0c0e7d4e13
# ╠═10452433-7123-4006-9e3d-ae4245fefdc5
# ╠═fd0b11ee-a0b8-4c0d-8c88-c69996b1c42d
# ╠═9b4ef541-9a36-4bc0-8654-10ab0a4e63b3
# ╠═18144929-3f31-42b2-9e27-df146a687ae0
# ╟─11e63c9a-6124-4122-9a86-ceed926d25d2
# ╠═708c818b-9b35-40e2-b7f0-dd9d5ebdfa85
# ╠═d95b8d0a-47dc-46ac-9ee8-dd8647dd6c7f
# ╠═f9a7488d-12a6-45f0-9c70-e67448dfe637
# ╠═6b02c425-39b9-467f-9406-3e9096873af4
# ╟─126c41e7-dd65-46c6-8c5b-2439f5624fd5
# ╟─17bdc97a-4a0b-4931-a5c6-866f0c814601
# ╠═01e43234-2194-451d-9010-176aa4799fdb
# ╟─ae8059bd-5b5b-4ff2-a6f0-5ce672bdd54d
# ╠═cf55c5cb-ead6-40b6-896a-8f7e01613a46
# ╟─c95ebbed-3d9a-4be2-943b-08c86923ad89
# ╠═61d745d8-5c74-479b-9698-cd50bb68b3c7
# ╟─e5738623-4ea8-4866-a1da-7e849960f4e0
# ╠═c34c1f0c-7cfd-404b-aa59-d4bb6aa9628f
# ╟─069497ac-fdae-4ddf-8983-026f7d46f07a
# ╠═cd132003-1384-41a4-bfb4-91247630a24e
# ╟─b075d6d8-228a-4ce2-8647-e2c6b962ba48
# ╠═a7d6a82c-b143-4e8e-94ee-e8999eefc0f1
# ╟─f2f949f8-772f-4787-8883-0d96137f0924
# ╟─3895472f-0a4f-4b7a-84f6-470208b5e8cc
# ╠═89a56a57-a5b9-4380-a618-97d8b901c01b
# ╠═0440ce41-8c23-45e6-aa67-56c6a58298a6
# ╟─cb927f51-99f2-4eeb-a0e5-5f1c65464b6f
# ╟─d7da11fd-8768-4ac7-81af-82f68e794e1a
# ╠═ab9b089f-fbb7-4436-8d6a-963db1e95670
# ╟─72204688-6346-4ff5-b443-839cbd7074d8
# ╠═6d09613a-5676-4b45-bf44-7e2c40bb71c9
# ╟─25de9867-d737-496c-8e0c-29bcdca898e8
# ╠═3525a716-0210-4cd2-b780-3f1c27297e09
# ╟─8774496e-a184-4c9a-9335-d2f184673cf5
# ╠═609937c6-8b9f-4987-8f46-ec55ce05e861
# ╟─aa5fff0c-9182-4568-b657-ca48c33de141
# ╠═7a6c247e-d765-4299-bd93-8d8271ca711f
# ╟─7103af68-9ad7-4b81-8c21-c8b6d7c9f5be
# ╠═56c7a1b9-c75f-48d8-a602-c219a2f432af
