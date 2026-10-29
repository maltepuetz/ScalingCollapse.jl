# E8: synthetic finite-size-scaling data with known parameters from test/generators/synthetic.jl,
# y = L^{±p3/p2} F((x − x_c)/x_c L^{1/p2}) (1 + c L^{−ω}) + noise with x_c = 0.7, ν = 1.2, p3 = 0.3,
# run with every quality function. Exact cases have δabs = 0 (the data have no finite-size bias);
# σmax is 3× the measured fit error, rounded to 2 significant digits. `corrections` (c = 0.5,
# ω = 1) is biased by construction, so it is a regression case plus an exact `_sanity` case. The
# sanity cases need no δabs: the largest measured deviations are p3 = 0.31353 ± 0.00651 (2.1σ,
# Linear) and 0.31497 ± 0.00572 (2.6σ, MultipleSplines); all others are below 2σ.
# SingleSpline fits a single spline to all sizes, so its minimum is broad: its σ on ν is
# 0.04 to 0.16 and the fitted ν lies up to 1.7σ above the truth, which makes its exact checks weak.
using HDF5

function load_synthetic(variant)
    h5open(joinpath(@__DIR__, "..", "testdata", "synthetic", "synthetic.h5"), "r") do h
        Ls = read(h, "$variant/Ls")
        g = h[variant]
        xs = [read(g, "L$L/x") for L in Ls]
        ys = [read(g, "L$L/y") for L in Ls]
        es = haskey(g, "L$(Ls[1])/e") ? [read(g, "L$L/e") for L in Ls] : nothing
        return xs, ys, es, Ls
    end
end

# σmax per free parameter, keyed `<variant>_<quality>`
const SYNTHETIC_SIGMAMAX = Dict{String,Vector{Float64}}(
    "clean_x_Linear" => [9.4e-5, 0.006],
    "clean_x_SingleSpline" => [0.0069, 0.46],
    "clean_x_MultipleSplines" => [6.7e-5, 0.0056],
    "clean_x_Houdayer" => [0.00068, 0.05],
    "clean_xy_Linear" => [9.0e-5, 0.0063, 0.0016],
    "clean_xy_SingleSpline" => [0.007, 0.43, 0.11],
    "clean_xy_MultipleSplines" => [6.1e-5, 0.0054, 0.0014],
    "clean_xy_Houdayer" => [0.00068, 0.057, 0.013],
    "clean_xny_Linear" => [0.00018, 0.01, 0.0055],
    "clean_xny_SingleSpline" => [0.0036, 0.36, 0.2],
    "clean_xny_MultipleSplines" => [0.00013, 0.0096, 0.0046],
    "clean_xny_Houdayer" => [0.0013, 0.1, 0.045],
    "noise_err_Linear" => [0.00073, 0.025, 0.027],
    "noise_err_SingleSpline" => [0.0053, 0.11, 0.26],
    "noise_err_MultipleSplines" => [0.00049, 0.019, 0.02],
    "noise_err_Houdayer" => [0.0011, 0.038, 0.041],
    "noise_noerr_Linear" => [0.0014, 0.083, 0.021],
    "noise_noerr_SingleSpline" => [0.0072, 0.45, 0.11],
    "noise_noerr_MultipleSplines" => [0.00085, 0.078, 0.019],
    "noise_noerr_Houdayer" => [0.0012, 0.077, 0.022],
    "three_sizes_Linear" => [0.00016, 0.0084, 0.0022],
    "three_sizes_SingleSpline" => [0.0054, 0.24, 0.069],
    "three_sizes_MultipleSplines" => [6.6e-6, 0.00067, 0.00013],
    "three_sizes_Houdayer" => [0.00051, 0.028, 0.0069],
    "uneven_grid_Linear" => [0.00012, 0.0099, 0.0026],
    "uneven_grid_SingleSpline" => [0.0062, 0.49, 0.097],
    "uneven_grid_MultipleSplines" => [6.6e-5, 0.0054, 0.0014],
    "uneven_grid_Houdayer" => [0.0012, 0.11, 0.026],
    "corrections_Linear" => [0.00049, 0.017, 0.02],
    "corrections_SingleSpline" => [0.0055, 0.11, 0.26],
    "corrections_MultipleSplines" => [0.00042, 0.016, 0.017],
    "corrections_Houdayer" => [0.0011, 0.038, 0.04],
    "fixed_nu_Linear" => [8.9e-5, 0.0015],
    "fixed_nu_SingleSpline" => [0.0073, 0.18],
    "fixed_nu_MultipleSplines" => [6.5e-5, 0.0014],
    "fixed_nu_Houdayer" => [0.00085, 0.015],
)

const SYNTHETIC_QUALITIES = ["Linear" => Linear, "SingleSpline" => SingleSpline,
    "MultipleSplines" => MultipleSplines, "Houdayer" => Houdayer]

# variant => (preset, true free-parameter values)
const SYNTHETIC_VARIANTS = [
    "clean_x" => (:x, [0.7, 1.2]),
    "clean_xy" => (:xy, [0.7, 1.2, 0.3]),
    "clean_xny" => (:xny, [0.7, 1.2, 0.3]),
    "noise_err" => (:xy, [0.7, 1.2, 0.3]),
    "noise_noerr" => (:xy, [0.7, 1.2, 0.3]),
    "three_sizes" => (:xy, [0.7, 1.2, 0.3]),
    "uneven_grid" => (:xy, [0.7, 1.2, 0.3]),
    "corrections" => (:xy, [0.7, 1.2, 0.3]),
    "fixed_nu" => (:xy, [0.7, 0.3]),  # ν = 1.2 fixed: free parameters are x_c and p3
]

function synthetic_kwargs(variant, preset, quality)
    if variant == "fixed_nu"
        return (sf=ScalingFunction(:xy; p2=1.2), p_space=[0.5:0.01:0.9, 0.0:0.05:0.6],
            dx=[-1.5, 1.5], quality=quality())
    end
    p_space = [0.5:0.01:0.9, 0.6:0.1:2.0, 0.0:0.05:0.6][1:(preset === :x ? 2 : 3)]
    return (sf=ScalingFunction(preset), p_space=p_space, dx=[-1.5, 1.5], quality=quality())
end

for (variant, (preset, truth)) in SYNTHETIC_VARIANTS, (qname, quality) in SYNTHETIC_QUALITIES
    name = "synthetic_$(variant)_$(qname)"
    kwargs = synthetic_kwargs(variant, preset, quality)
    load = () -> load_synthetic(variant)
    source = "synthetic $variant data (test/generators/synthetic.jl), $qname quality"
    σmax = Dict(i => s for (i, s) in enumerate(SYNTHETIC_SIGMAMAX["$(variant)_$qname"]))
    exact = Dict(i => p for (i, p) in enumerate(truth))
    if variant == "corrections"
        push!(CASES, CollapseCase(; name, kind=:regression, load, kwargs, source))
        push!(CASES, CollapseCase(; name=name * "_sanity", kind=:exact, load, kwargs, exact, σmax, source))
    else
        push!(CASES, CollapseCase(; name, kind=:exact, load, kwargs, exact, σmax, source))
    end
end
