# E9: noisy synthetic data from the quality-function study (PR #63) behind pchip interpolation in
# MultipleSplines (test/generators/quality_study.jl): irregular x grids, 6 points per size, a
# narrow peak and a curve with several extrema, each with 1% noise and exact error bars. Each data
# set is replicate 1 of its study configuration (not selected for a pass). The cases run the
# default quality, MultipleSplines with pchip; σmax is 3× the measured fit error, rounded to 2
# significant digits. On these replicates every fit lies within 0.6σ of the truth.
#
# `random_order_12_cubic` documents why the cubic spline is no longer the default: on 12 noisy
# random points per size it oscillates between close points, and the fit ends at x_c = 2.40,
# ν = 1.23 and p3 = 0.010 (truth 2.27, 1.0 and 0.125) with NaN errors.
using HDF5

function load_quality_study(name)
    h5open(joinpath(@__DIR__, "..", "testdata", "synthetic", "quality_study.h5"), "r") do h
        g = h[name]
        Ls = read(g, "Ls")
        xs, ys, es = ([read(g, "L$L/$v") for L in Ls] for v in ("x", "y", "e"))
        return xs, ys, es, Ls
    end
end

# p_space of the study: a grid around (not centred on) the true values
function quality_study_p_space(truth)
    xc = truth[1]
    ps = [range(round(0.8xc; digits=2), round(1.25xc; digits=2); length=46), 0.3:0.05:3.0]
    length(truth) == 3 && push!(ps, 0.0:0.05:2.5)
    return ps
end

# name => (preset, truth, σmax)
const QUALITY_STUDY_CASES = [
    "random_crossing" => (:x, [0.7, 1.2], [0.0024, 0.08]),
    "random_order_12" => (:xy, [2.27, 1.0, 0.125], [0.0033, 0.059, 0.051]),
    "random_peak" => (:xny, [2.27, 1.0, 1.75], [0.001, 0.037, 0.036]),
    "few_points_crossing" => (:x, [0.7, 1.2], [0.0019, 0.067]),
    "few_points_peak" => (:xny, [2.27, 1.0, 1.75], [0.0037, 0.12, 0.13]),
    "narrow_peak" => (:xny, [2.27, 1.0, 1.75], [0.0012, 0.087, 0.14]),
    "wiggly" => (:x, [0.7, 1.2], [0.0032, 0.086]),
]

for (variant, (preset, truth, σmax)) in QUALITY_STUDY_CASES
    push!(CASES, CollapseCase(;
        name="quality_study_$variant", kind=:exact,
        load=() -> load_quality_study(variant),
        kwargs=(sf=ScalingFunction(preset), p_space=quality_study_p_space(truth)),
        exact=Dict(i => p for (i, p) in enumerate(truth)),
        σmax=Dict(i => s for (i, s) in enumerate(σmax)),
        source="synthetic $variant data (test/generators/quality_study.jl), default quality",
    ))
end

push!(CASES, CollapseCase(;
    name="quality_study_random_order_12_cubic", kind=:exact, broken=true,
    load=() -> load_quality_study("random_order_12"),
    kwargs=(sf=ScalingFunction(:xy), p_space=quality_study_p_space([2.27, 1.0, 0.125]),
        quality=MultipleSplines(interp=:cubic),
        quality_scan=MultipleSplines(scan_mode=true, interp=:cubic)),
    exact=Dict(1 => 2.27, 2 => 1.0, 3 => 0.125),
    source="synthetic random_order_12 data (test/generators/quality_study.jl), cubic MultipleSplines",
))
