# E8: synthetic finite-size-scaling data with known parameters from test/generators/synthetic.jl,
# y = L^{±p3/p2} F((x − x_c)/x_c L^{1/p2}) (1 + c L^{−ω}) + noise with x_c = 0.7, ν = 1.2, p3 = 0.3,
# run with every quality function. Exact cases have δabs = 0 (the data have no finite-size bias);
# σmax is 3× the measured fit error, rounded to 2 significant digits. `corrections` (c = 0.5,
# ω = 1) is biased by construction, so it is a regression case plus an exact `_sanity` case. The
# sanity cases need no δabs: the largest measured deviation is p3 = 0.31497 ± 0.01059 (1.4σ,
# MultipleSplines); all others are below 1.1σ.
# Houdayer's ν lies 2–3% above the truth on these coarse grids; this is the grid-spacing bias
# of the method, see the `Houdayer` docstring.
# The window dx = [-1.05, 1.05] is [-1.5, 1.5] of the relative form of the generator, times x_c.
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
    "clean_x_SingleSpline" => [3.1e-5, 0.0029],
    "clean_x_MultipleSplines" => [6.6e-5, 0.0056],
    "clean_x_Houdayer" => [0.00069, 0.05],
    "clean_xy_Linear" => [0.00018, 0.0093, 0.0023],
    "clean_xy_SingleSpline" => [4.4e-5, 0.0024, 0.00066],
    "clean_xy_MultipleSplines" => [0.00012, 0.0083, 0.0022],
    "clean_xy_Houdayer" => [0.0017, 0.098, 0.021],
    "clean_xny_Linear" => [0.00018, 0.012, 0.0066],
    "clean_xny_SingleSpline" => [4.3e-5, 0.0026, 0.0015],
    "clean_xny_MultipleSplines" => [0.00013, 0.01, 0.0048],
    "clean_xny_Houdayer" => [0.0014, 0.11, 0.051],
    "noise_err_Linear" => [0.0031, 0.084, 0.054],
    "noise_err_SingleSpline" => [0.0022, 0.067, 0.053],
    "noise_err_MultipleSplines" => [0.0015, 0.043, 0.037],
    "noise_err_Houdayer" => [0.0034, 0.1, 0.079],
    "noise_noerr_Linear" => [0.0028, 0.13, 0.031],
    "noise_noerr_SingleSpline" => [0.0022, 0.11, 0.036],
    "noise_noerr_MultipleSplines" => [0.0016, 0.11, 0.03],
    "noise_noerr_Houdayer" => [0.0021, 0.1, 0.033],
    "three_sizes_Linear" => [0.00039, 0.015, 0.0038],
    "three_sizes_SingleSpline" => [8.6e-5, 0.0037, 0.0009],
    "three_sizes_MultipleSplines" => [6.6e-6, 0.00067, 0.00013],
    "three_sizes_Houdayer" => [0.0012, 0.051, 0.011],
    "uneven_grid_Linear" => [0.00028, 0.018, 0.0037],
    "uneven_grid_SingleSpline" => [6.6e-5, 0.0049, 0.0009],
    "uneven_grid_MultipleSplines" => [0.00013, 0.009, 0.002],
    "uneven_grid_Houdayer" => [0.0033, 0.25, 0.044],
    "corrections_Linear" => [0.0019, 0.045, 0.04],
    "corrections_SingleSpline" => [0.0018, 0.049, 0.048],
    "corrections_MultipleSplines" => [0.0013, 0.036, 0.032],
    "corrections_Houdayer" => [0.0024, 0.075, 0.079],
    "fixed_nu_Linear" => [0.00013, 0.0015],
    "fixed_nu_SingleSpline" => [2.2e-5, 0.00039],
    "fixed_nu_MultipleSplines" => [9.3e-5, 0.0014],
    "fixed_nu_Houdayer" => [0.0013, 0.015],
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
            dx=[-1.05, 1.05], quality=quality())
    end
    p_space = [0.5:0.01:0.9, 0.6:0.1:2.0, 0.0:0.05:0.6][1:(preset === :x ? 2 : 3)]
    return (sf=ScalingFunction(preset), p_space=p_space, dx=[-1.05, 1.05], quality=quality())
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
