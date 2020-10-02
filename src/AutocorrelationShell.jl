__precompile__()

module AutocorrelationShell

export
    # wavelet transforms
    acwt,
    iacwt,
    fwt_ac,
    iwt_ac,
    # filters
    pfilter,
    qfilter,
    # others
    autocorr,
    ac_filter,
    echant,
    autocorr_calc,
    dyadlength,

    # 2D forward autocorrelation wavelet transform
    ac2d,
    acwt2D,
    # inverse 2D autocorrelation wavelet transform
    iac2d,
    iacwt2D,

    # threshold
    acthreshold,

    # modify image
    make_noisy,
    # signal to noise ratio
    snr,
    get_snr,

    acwt_heatmap,
    # entropy
    wentropy,
    NormEntropy,
    ThresholdEntropy,

    # wiggle
    wiggle,
    wiggle!,

    # AC Wavelet Packets
    BinaryNode,
    acwpt,
    iacwpt,
    acwptPostOrderBestBasis,
    acwptPreOrderBestBasis

using
    AbstractTrees,
    DSP,
    Plots,
    StatsBase,
    Statistics,
    Wavelets,
    LinearAlgebra,
    Random

# the main show
include("ACW1D.jl")
include("ACW2D.jl")
include("ACWUtil.jl")
include("ACWPT1D.jl")

end # module
