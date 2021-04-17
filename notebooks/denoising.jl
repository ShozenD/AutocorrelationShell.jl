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
	Pkg.instantiate()
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
md"## I. Exploratory Data Analysis"

# ╔═╡ bc641876-9723-4c68-8221-ad03fb695c82
md"Let's load some test functions"

# ╔═╡ d0c98a14-ad3e-4a0b-b889-a3ea86f888f3
testdata = CSV.read("../test/data/wavelet_test_256.csv", DataFrame)

# ╔═╡ b0a28618-7cda-4c05-83d2-b54bbca3f9b5
md"**Select** a test function"

# ╔═╡ 7364da28-6a01-4359-9664-a3097e8bf1f1
@bind signal_name_test Select(
	["blocks", "bumps", "heavy_sine", "doppler", "quadchirp", "mishmash"],
	default = "doppler"
)

# ╔═╡ 458a5a1e-c453-4199-befe-2bf4db6825ae
md"**Adjust** the magnitude of Gaussian noise"

# ╔═╡ 97f40df0-9ccc-4e41-bebf-4e7188f33fff
@bind noise_size_test Slider(0:0.01:1, default=0.3, show_value=true)

# ╔═╡ 7c0dba17-baf3-4b9c-b1c5-486f7e4515f4
md"The `addnoise` function will add Gaussian noise that is proportional to the total energy of the signal."

# ╔═╡ a184ae65-7947-4ffb-b751-8b6e97a8608b
function addnoise(x::AbstractArray{<:Number,1}, s::Real=0.1)
	ϵ = randn(length(x))
	ϵ = ϵ/norm(ϵ)
	y = x + ϵ * s * norm(x)
	return y
end

# ╔═╡ 851b04bb-e382-4a0a-98f6-7d4b983ca5ab
begin
	x_test = testdata[!,signal_name_test]
	p1_test = Plots.plot(x_test, ylim = (minimum(x_test)-1,maximum(x_test)+1), label = "Original signal")
	x_noisy_test = addnoise(x_test, noise_size_test)
	Plots.plot!(x_noisy_test, label = "Noisy signal");
	
	y_test = acwt(x_noisy_test, wavelet(WT.db4), 4);
	p2_test = WaveletsExt.wiggle(y_test)
	
	Plots.plot(p1_test, p2_test, layout = (2,1))
end

# ╔═╡ 356d75f4-6cc1-4062-8ef6-3cc6a9c2d9a7
md"Plot a histogram of wavelet coefficients"

# ╔═╡ 341c2551-b625-47a0-9163-3c0c0e7d4e13
Plots.histogram(vec(abs.(y_test)), legend = false)

# ╔═╡ 261c2822-edaf-4c66-a032-c17fc2447627
md"Threshold the coefficients using an arbitrary threshold value"

# ╔═╡ 10452433-7123-4006-9e3d-ae4245fefdc5
threshold!(y_test, HardTH(), 0.3);

# ╔═╡ a663045c-fa0f-49fe-88c2-794450cb7806
md"Reconstruct the signal using the thresholded coefficients"

# ╔═╡ 9b4ef541-9a36-4bc0-8654-10ab0a4e63b3
begin
	r_test = iacwt(y_test)
	Plots.plot(x_test, label = "original")
	Plots.plot!(r_test, label = "denoised")
end

# ╔═╡ 4669be94-6c4c-42e2-b9d9-2dc98f1bdaea
md"Calculate the Mean Squared Error between the original signal and the denoised signal"

# ╔═╡ 18144929-3f31-42b2-9e27-df146a687ae0
norm(x_test - r_test)/length(x_test)

# ╔═╡ 11e63c9a-6124-4122-9a86-ceed926d25d2
md"# II. Data Setup"

# ╔═╡ d881753b-0432-451b-8de0-38a0b4b4382a
md"**Autorun**: Please disable autorun before updating the experiment parameters, else it will run the entire notebook, which may take a few minutes."

# ╔═╡ 7e94d13e-f84c-433c-bead-3a272c86fc9b
@bind autorun Radio(["No","Yes"], default = "No")

# ╔═╡ 8055194b-2e46-4d18-81c0-0c52bc3eb233
md"**Select** a test function"

# ╔═╡ 18b5bbe4-ecdd-4209-a764-7c8b1ecbda61
@bind signal_name Radio(
	[
		"blocks" => "Blocks", 
		"bumps" => "Bumps", 
		"heavy_sine" => "Heavy sine", 
		"doppler" => "Doppler", 
		"quadchirp" => "Quadchirp", 
		"mishmash" => "Mishmash"
	],
	default = "blocks"
)

# ╔═╡ 56ee2c61-d83c-4d76-890a-a9bd0d65cee5
md"**Adjust** the slider to add Gaussian noise to the test signal"

# ╔═╡ c50ac92e-3684-4d0a-a80d-4ee9d74ec992
@bind noise_size Slider(0:0.01:1, default=0.3, show_value=true)

# ╔═╡ e0a96592-5e77-4c29-9744-31369eea8147
md"**Select** which type of wavelet basis to use"

# ╔═╡ 53557a60-90f5-48e6-81ed-5736fc05fec0
@bind wavelet_type Select(
	["WT.haar", 
	"WT.db1", "WT.db2", "WT.db3", "WT.db4", "WT.db5", 
	"WT.db6", "WT.db7", "WT.db8", "WT.db9", "WT.db10",
	"WT.coif2", "WT.coif4", "WT.coif6", "WT.coif8",
	"WT.sym4", "WT.sym5", "WT.sym6", "WT.sym7", "WT.sym8", "WT.sym9", "WT.sym10",
	"WT.batt2", "WT.batt4", "WT.batt6"],
	default = "WT.haar"
)

# ╔═╡ c178527f-96a4-4ac7-bb0c-38b73b38c45b
md"**Select** which type of thresholding to use"

# ╔═╡ f9a7488d-12a6-45f0-9c70-e67448dfe637
@bind threshold_method Radio(
	["HardTH()" => "Hard",
	"SoftTH()" => "Soft"],
	default = "HardTH()"
)

# ╔═╡ cd9e259e-8bb3-497b-ac7f-f89a003c8032
begin
	x = testdata[!,signal_name]
	p3 = Plots.plot(x, ylim = (minimum(x)-1,maximum(x)+1), label = "Original signal")
	x_noisy = addnoise(x, noise_size)
	Plots.plot!(x_noisy, label = "Noisy signal");
end

# ╔═╡ 3246e8b5-251f-4398-b21c-397341f2542e
md"**Preparing** the data
* Generate a set of original signals $\rightarrow$ `X₀`
* Add noise to original signals $\rightarrow$ `X`
"

# ╔═╡ 82e713f8-c870-43d2-a849-e3b401b00459
begin
	X₀ = generatesignals(x, 100, 2)
	X = hcat([addnoise(X₀[:,i], noise_size) for i in axes(X₀,2)]...)
end

# ╔═╡ 6664a859-4980-4b40-8684-83cf2e7db109
md"**Baseline**: No denoising"

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
md"# III. Non-Redundant Transforms"

# ╔═╡ 17bdc97a-4a0b-4931-a5c6-866f0c814601
md"### 1. Discrete Wavelet Transform"

# ╔═╡ 01e43234-2194-451d-9010-176aa4799fdb
begin
	if autorun == "Yes"
		Y["DWT"] = cat([wpt(X[:,i], wt) for i in axes(X,2)]..., dims=2)
		X̂["DWT"], time["DWT"] = @timed denoiseall(
			Y["DWT"], :dwt, wt, L=8, dnt=dnt
		)
		mean_psnr["DWT"] = mean([psnr(X̂["DWT"][:,i], X₀[:,i]) for i in axes(X,2)])
		push!(results, ["DWT", threshold_method, time["DWT"], mean_psnr["DWT"]])
	end
end

# ╔═╡ ae8059bd-5b5b-4ff2-a6f0-5ce672bdd54d
md"### 2. Wavelet Packet Transform - Level 4"

# ╔═╡ cf55c5cb-ead6-40b6-896a-8f7e01613a46
begin
	if autorun == "Yes"
		Y["WPT-L4"] = cat([wpt(X[:,i], wt, 4) for i in axes(X,2)]..., dims=2)
		X̂["WPT-L4"], time["WPT-L4"] = @timed denoiseall(
			Y["WPT-L4"], :wpt, wt, tree=maketree(256,4,:full), dnt=dnt
		)
		mean_psnr["WPT-L4"] = mean([psnr(X̂["WPT-L4"][:,i], X₀[:,i]) for i in axes(X,2)])
		push!(results, ["WPT-L4", threshold_method, time["WPT-L4"], mean_psnr["WPT-L4"]])
	end
end

# ╔═╡ c95ebbed-3d9a-4be2-943b-08c86923ad89
md"### 3. Wavelet Packet Transform - Level 8"

# ╔═╡ 61d745d8-5c74-479b-9698-cd50bb68b3c7
begin
	if autorun == "Yes"
		Y["WPT-L8"] = cat([wpt(X[:,i], wt) for i in axes(X,2)]..., dims=2)
		X̂["WPT-L8"], time["WPT-L8"] = @timed denoiseall(
			Y["WPT-L8"], :wpt, wt, tree=maketree(256,8,:full), dnt=dnt
		)
		mean_psnr["WPT-L8"] = mean([psnr(X̂["WPT-L8"][:,i], X₀[:,i]) for i in axes(X,2)])
		push!(results, ["WPT-L8", threshold_method, time["WPT-L8"], mean_psnr["WPT-L8"]])
	end
end

# ╔═╡ e5738623-4ea8-4866-a1da-7e849960f4e0
md"### 4. Wavelet Packet Transform - Best Basis"

# ╔═╡ c34c1f0c-7cfd-404b-aa59-d4bb6aa9628f
begin
	if autorun == "Yes"
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
end

# ╔═╡ 069497ac-fdae-4ddf-8983-026f7d46f07a
md"### 5. Joint Best Basis"

# ╔═╡ cd132003-1384-41a4-bfb4-91247630a24e
begin
	if autorun == "Yes"
		jbb_tree = bestbasistree(Y["WPD"], JBB())
		Y["JBB"] = cat([wpt(X[:,i], wt, jbb_tree) for i in axes(X,2)]..., dims=2)
		X̂["JBB"], time["JBB"] = @timed denoiseall(
			Y["JBB"], :wpt, wt, tree=jbb_tree, dnt=dnt
		)
		mean_psnr["JBB"] = mean([psnr(X̂["JBB"][:,i], X₀[:,i]) for i in axes(X,2)])
		push!(results, ["JBB", threshold_method, time["JBB"], mean_psnr["JBB"]])
	end
end

# ╔═╡ b075d6d8-228a-4ce2-8647-e2c6b962ba48
md"### 6. Least Statistically Dependent Basis"

# ╔═╡ a7d6a82c-b143-4e8e-94ee-e8999eefc0f1
begin 
	if autorun == "Yes"
		lsdb_tree = bestbasistree(Y["WPD"], LSDB())
		Y["LSDB"] = cat([wpt(X[:,i], wt, lsdb_tree) for i in axes(X,2)]..., dims=2)
		X̂["LSDB"], time["LSDB"] = @timed denoiseall(
			Y["LSDB"], :wpt, wt, tree=lsdb_tree, dnt=dnt
		)
		mean_psnr["LSDB"] = mean([psnr(X̂["LSDB"][:,i], X₀[:,i]) for i in axes(X,2)])
		push!(results, ["LSDB", threshold_method, time["LSDB"], mean_psnr["LSDB"]])
	end
end

# ╔═╡ f2f949f8-772f-4787-8883-0d96137f0924
md"# IV. Redundant Transforms"

# ╔═╡ 3895472f-0a4f-4b7a-84f6-470208b5e8cc
md"## 1. Autocorrelation Wavelet Transforms"

# ╔═╡ 7c3ae1ea-887d-4af6-ba18-7fd06ea6354d
md"### 1.1 Autocorrelation Discrete Wavelet Transform"

# ╔═╡ 89a56a57-a5b9-4380-a618-97d8b901c01b
begin
	if autorun == "Yes"
		Y["ACWT"] = cat([acwt(X[:,i], wt) for i in axes(X,2)]..., dims=3)
		X̂["ACWT"], time["ACWT"] = @timed denoiseall(
			Y["ACWT"], :acwt, wt, L=8, dnt=dnt
		)
		mean_psnr["ACWT"] = mean([psnr(X̂["ACWT"][:,i], X₀[:,i]) for i in axes(X,2)])
		push!(results, ["ACWT", threshold_method, time["ACWT"], mean_psnr["ACWT"]])
	end
end

# ╔═╡ 115edde2-b1ba-4d86-9b1a-e05d76026bcf
md"### 1.2 Autocorrelation Packet Transform - Level 4"

# ╔═╡ c52bd741-00cb-4cf2-97e3-b8dbba3af9ad
begin
	if autorun == "Yes"
		Y["ACWPT"] = cat([acwpt(X[:,i], wt) for i in axes(X,2)]..., dims=3)
		X̂["ACWPT-L4"], time["ACWPT-L4"] = @timed denoiseall(
			Y["ACWPT"], :acwpt, wt, tree=maketree(256,4,:full), dnt=dnt
		)
		mean_psnr["ACWPT-L4"] = mean([psnr(X̂["ACWPT-L4"][:,i], X₀[:,i]) for i in axes(X,2)])
		push!(
			results, 
			["ACWPT-L4", threshold_method, time["ACWPT-L4"], mean_psnr["ACWPT-L4"]]
		)
	end
end

# ╔═╡ a4b3694e-2ca5-46a5-b1ce-44c2f7ec2006
md"### 1.3 Autocorrelation Packet Transform - Level 8"

# ╔═╡ 250d22dd-1d15-45d4-8aa4-3de1f37b164c
begin
	if autorun == "Yes"
		Y["ACWPT"] = cat([acwpt(X[:,i], wt) for i in axes(X,2)]..., dims=3)
		X̂["ACWPT-L8"], time["ACWPT-L8"] = @timed denoiseall(
			Y["ACWPT"], :acwpt, wt, tree=maketree(256,8,:full), dnt=dnt
		)
		mean_psnr["ACWPT-L8"] = mean([psnr(X̂["ACWPT-L8"][:,i], X₀[:,i]) for i in axes(X,2)])
		push!(
			results, 
			["ACWPT-L8", threshold_method, time["ACWPT-L8"], mean_psnr["ACWPT-L8"]]
		)
	end
end

# ╔═╡ 7d057183-c8b2-4ebd-9ee9-aa7998d9a6d5
md"### 1.4 Autocorrelation Packet Transform - Best Basis"

# ╔═╡ 991b273a-9e25-4499-b6ea-4d800ea1e6ae
begin
	if autorun == "Yes"
		acwpt_bt = bestbasistree(Y["ACWPT"], BB(stationary=true))
		X̂["ACWPT-BT"], time["ACWPT-BT"] = @timed hcat(
			[denoise(
				Y["ACWPT"][:,:,i], :acwpt, wt, tree=acwpt_bt[:,i], dnt=dnt
				) for i in axes(X,2)
			]...
		)
		mean_psnr["ACWPT-BT"] = mean([psnr(X̂["ACWPT-BT"][:,i], X₀[:,i]) for i in axes(X,2)])
		push!(
			results, 
			["ACWPT-BT", threshold_method, time["ACWPT-BT"], mean_psnr["ACWPT-BT"]]
		)
	end
end

# ╔═╡ 084decfb-5c6f-466d-a5f8-ddf8cc863d8c
md"### 1.5 Autocorrelation Joint Best Basis"

# ╔═╡ c5a90584-fc46-4f7e-8633-6866001dadf6
begin
	if autorun == "Yes"
		ajbb_tree = bestbasistree(Y["ACWPT"], JBB(stationary=true))
		X̂["ACJBB"], time["ACJBB"] = @timed denoiseall(
			Y["ACWPT"], :acwpt, wt, tree=ajbb_tree, dnt=dnt
		)
		mean_psnr["ACJBB"] = mean([psnr(X̂["ACJBB"][:,i], X₀[:,i]) for i in axes(X,2)])
		push!(results, ["ACJBB", threshold_method, time["ACJBB"], mean_psnr["ACJBB"]])
	end
end

# ╔═╡ cb927f51-99f2-4eeb-a0e5-5f1c65464b6f
md"## 2. Stationary Wavelet Transforms"

# ╔═╡ d7da11fd-8768-4ac7-81af-82f68e794e1a
md"### 2.1 Stationary Discrete Wavelet Transform"

# ╔═╡ ab9b089f-fbb7-4436-8d6a-963db1e95670
begin
	if autorun == "Yes"
		Y["SDWT"] = cat([sdwt(X[:,i], wt) for i in axes(X,2)]..., dims=3)
		X̂["SDWT"], time["SDWT"] = @timed denoiseall(
			Y["SDWT"], :sdwt, wt, L=8, dnt=dnt
		)
		mean_psnr["SDWT"] = mean([psnr(X̂["SDWT"][:,i], X₀[:,i]) for i in axes(X,2)])
		push!(results, ["SDWT", threshold_method, time["SDWT"], mean_psnr["SDWT"]])
	end
end

# ╔═╡ 72204688-6346-4ff5-b443-839cbd7074d8
md"### 2.2 Stationary Wavelet Packet Transform - Level 4"

# ╔═╡ 6d09613a-5676-4b45-bf44-7e2c40bb71c9
begin
	if autorun == "Yes"
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
end

# ╔═╡ 25de9867-d737-496c-8e0c-29bcdca898e8
md"### 2.3 Stationary Wavelet Packet Transform - Level 8"

# ╔═╡ 3525a716-0210-4cd2-b780-3f1c27297e09
begin
	if autorun == "Yes"
		X̂["SWPT-L8"], time["SWPT-L8"] = @timed denoiseall(
			Y["SWPD"], :swpd, wt, tree=maketree(256,8,:full), dnt=dnt
		)
		mean_psnr["SWPT-L8"] = mean([psnr(X̂["SWPT-L8"][:,i], X₀[:,i]) for i in axes(X,2)])
		push!(
			results, 
			["SWPT-L8", threshold_method, time["SWPT-L8"], mean_psnr["SWPT-L8"]]
		)
	end
end

# ╔═╡ 8774496e-a184-4c9a-9335-d2f184673cf5
md"### 2.4 Stationary Wavelet Packet Transform - Best Basis"

# ╔═╡ 609937c6-8b9f-4987-8f46-ec55ce05e861
begin
	if autorun == "Yes"
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
end

# ╔═╡ aa5fff0c-9182-4568-b657-ca48c33de141
md"### 2.5 Stationary Joint Best Basis"

# ╔═╡ 7a6c247e-d765-4299-bd93-8d8271ca711f
begin
	if autorun == "Yes"
		sjbb_tree = bestbasistree(Y["SWPD"], JBB(stationary=true))
		X̂["SJBB"], time["SJBB"] = @timed denoiseall(
			Y["SWPD"], :swpd, wt, tree=sjbb_tree, dnt=dnt
		)
		mean_psnr["SJBB"] = mean([psnr(X̂["SJBB"][:,i], X₀[:,i]) for i in axes(X,2)])
		push!(results, ["SJBB", threshold_method, time["SJBB"], mean_psnr["SJBB"]])
	end
end

# ╔═╡ 7103af68-9ad7-4b81-8c21-c8b6d7c9f5be
md"## V. Results Analysis"

# ╔═╡ 9818ca36-0acc-4cc3-924e-8b50813c1da1
md"**Selected Test Parameters:**

* Test signal: $signal_name

* Noise magnitude: $noise_size

* Wavelet type: $wavelet_type

* Threshold method: $threshold_method

**Note**: You might need to re-run the block below after updating a parameter to display the results>
"

# ╔═╡ 56c7a1b9-c75f-48d8-a602-c219a2f432af
if autorun == "Yes"
	Gadfly.plot(
		results, 
		y=:transform, 
		x=:PSNR, 
		color=:threshold, 
		Guide.title("PSNR")
	)
end

# ╔═╡ Cell order:
# ╟─53c257e0-96ba-11eb-3615-8bfed63b2c18
# ╟─ce4bf94e-edef-40c2-8ac5-b741e47a1759
# ╠═f9e001ae-0323-4283-9af1-1a8252b503e7
# ╟─85c26ead-f043-42c8-8245-58c8d03a963d
# ╠═203fb8c8-4358-4908-b616-a691ce329c02
# ╟─46e30602-a850-4175-b4bf-b4ef4b5359aa
# ╟─bc641876-9723-4c68-8221-ad03fb695c82
# ╠═d0c98a14-ad3e-4a0b-b889-a3ea86f888f3
# ╟─b0a28618-7cda-4c05-83d2-b54bbca3f9b5
# ╟─7364da28-6a01-4359-9664-a3097e8bf1f1
# ╟─458a5a1e-c453-4199-befe-2bf4db6825ae
# ╠═97f40df0-9ccc-4e41-bebf-4e7188f33fff
# ╟─7c0dba17-baf3-4b9c-b1c5-486f7e4515f4
# ╠═a184ae65-7947-4ffb-b751-8b6e97a8608b
# ╠═851b04bb-e382-4a0a-98f6-7d4b983ca5ab
# ╟─356d75f4-6cc1-4062-8ef6-3cc6a9c2d9a7
# ╠═341c2551-b625-47a0-9163-3c0c0e7d4e13
# ╟─261c2822-edaf-4c66-a032-c17fc2447627
# ╠═10452433-7123-4006-9e3d-ae4245fefdc5
# ╟─a663045c-fa0f-49fe-88c2-794450cb7806
# ╠═9b4ef541-9a36-4bc0-8654-10ab0a4e63b3
# ╟─4669be94-6c4c-42e2-b9d9-2dc98f1bdaea
# ╠═18144929-3f31-42b2-9e27-df146a687ae0
# ╟─11e63c9a-6124-4122-9a86-ceed926d25d2
# ╟─d881753b-0432-451b-8de0-38a0b4b4382a
# ╟─7e94d13e-f84c-433c-bead-3a272c86fc9b
# ╟─8055194b-2e46-4d18-81c0-0c52bc3eb233
# ╟─18b5bbe4-ecdd-4209-a764-7c8b1ecbda61
# ╟─56ee2c61-d83c-4d76-890a-a9bd0d65cee5
# ╠═c50ac92e-3684-4d0a-a80d-4ee9d74ec992
# ╟─e0a96592-5e77-4c29-9744-31369eea8147
# ╟─53557a60-90f5-48e6-81ed-5736fc05fec0
# ╟─c178527f-96a4-4ac7-bb0c-38b73b38c45b
# ╠═f9a7488d-12a6-45f0-9c70-e67448dfe637
# ╠═cd9e259e-8bb3-497b-ac7f-f89a003c8032
# ╟─3246e8b5-251f-4398-b21c-397341f2542e
# ╠═82e713f8-c870-43d2-a849-e3b401b00459
# ╟─6664a859-4980-4b40-8684-83cf2e7db109
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
# ╟─7c3ae1ea-887d-4af6-ba18-7fd06ea6354d
# ╠═89a56a57-a5b9-4380-a618-97d8b901c01b
# ╟─115edde2-b1ba-4d86-9b1a-e05d76026bcf
# ╠═c52bd741-00cb-4cf2-97e3-b8dbba3af9ad
# ╟─a4b3694e-2ca5-46a5-b1ce-44c2f7ec2006
# ╠═250d22dd-1d15-45d4-8aa4-3de1f37b164c
# ╟─7d057183-c8b2-4ebd-9ee9-aa7998d9a6d5
# ╠═991b273a-9e25-4499-b6ea-4d800ea1e6ae
# ╟─084decfb-5c6f-466d-a5f8-ddf8cc863d8c
# ╠═c5a90584-fc46-4f7e-8633-6866001dadf6
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
# ╟─9818ca36-0acc-4cc3-924e-8b50813c1da1
# ╠═56c7a1b9-c75f-48d8-a602-c219a2f432af
