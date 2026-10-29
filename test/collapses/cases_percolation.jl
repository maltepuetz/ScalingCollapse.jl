# 2D percolation: self-generated bond percolation (test/generators/percolation.jl) and the
# site-percolation order parameter shipped with autoScale.py (test/testdata/autoscale, MIT).
using HDF5

# curves `<obs>/L{L}/{x, y, e}` of test/testdata/percolation/bond_square.h5
function load_percolation_h5(obs)
    h5open(joinpath(@__DIR__, "..", "testdata", "percolation", "bond_square.h5"), "r") do h
        Ls = read(h, "Ls")
        return ([read(h, "$obs/L$L/$v") for L in Ls] for v in ("x", "y", "e"))..., Ls
    end
end

# autoScale's orderParam_L{L}.dat: columns p, y, error; `#` lines are comments
function load_autoscale(Ls)
    Ls = collect(Ls)
    cols = map(Ls) do L
        path = joinpath(@__DIR__, "..", "testdata", "autoscale", "orderParam_L$L.dat")
        rows = [parse.(Float64, split(l)) for l in eachline(path) if !startswith(l, "#") && !isempty(strip(l))]
        return [getindex.(rows, j) for j in 1:3]
    end
    return getindex.(cols, 1), getindex.(cols, 2), getindex.(cols, 3), Ls
end

const PERCOLATION_PC_SITE = 0.59274605079210

# E4: Newman–Ziff bond percolation on the open L × L square, L = 16 to 256, p ∈ [0.45, 0.55]
# (p_c = 1/2, ν = 4/3, β = 5/36). σmax ≈ 3× the measured errors.
push!(CASES, CollapseCase(;
    name="bond_spanning_exact", kind=:exact,
    load=() -> load_percolation_h5("spanning"),
    kwargs=(sf=ScalingFunction(:x), p_space=[0.48:0.002:0.52, 0.8:0.1:2.0], dx=[-1.0, 1.0]),
    exact=Dict(1 => 0.5, 2 => 4 / 3),
    # measured: p_c = 0.501317 ± 0.000502 (2.6σ), ν = 1.3058 ± 0.0515 (0.5σ)
    σmax=Dict(1 => 1.6e-3, 2 => 0.16),
    source="Newman–Ziff bond percolation, top-bottom spanning probability, L = 16 to 256",
))
push!(CASES, CollapseCase(;
    name="bond_strength_exact", kind=:exact,
    load=() -> load_percolation_h5("strength"),
    kwargs=(sf=ScalingFunction(:xy), p_space=[0.48:0.002:0.52, 0.8:0.1:2.0, 0.05:0.025:0.3], dx=[-1.0, 1.0]),
    exact=Dict(1 => 0.5, 2 => 4 / 3, 3 => 5 / 36),
    # Known finite-size bias: the largest-cluster fraction has strong corrections to scaling at
    # these small sizes, which the quality-ratio error does not cover. δabs = 1.5 × the measured
    # |p − p*| on the parameters that fail at 3σ. Measured: p_c = 0.499720 ± 0.000138 (2.0σ),
    # ν = 1.2802 ± 0.0136 (|dev| 0.0532, 3.9σ), β = 0.15530 ± 0.00306 (|dev| 0.0164, 5.4σ)
    δabs=Dict(2 => 0.080, 3 => 0.025),
    σmax=Dict(1 => 4.2e-4, 2 => 0.041, 3 => 9.2e-3),
    source="Newman–Ziff bond percolation, largest-cluster fraction, L = 16 to 256",
))

# E3: autoScale's own fit uses L = 128, 256, 512 (its inputFiles.dat) and the window
# [-1.5, 1.0] of (p − p_c)·L^{1/ν} (its README); here x is relative, (p − p1)/p1·L^{1/p2}, so the
# window is divided by p_c. autoScale itself reports p_c = 0.592708 ± 0.000075,
# 1/ν = 0.7476 ± 0.0054, β/ν = 0.1045 ± 0.0010 (its README).
push!(CASES, CollapseCase(;
    name="autoscale_site_order_exact", kind=:exact,
    load=() -> load_autoscale((128, 256, 512)),
    kwargs=(sf=ScalingFunction(:xy), p_space=[0.58:0.002:0.61, 0.8:0.1:2.0, 0.05:0.025:0.3],
        dx=[-1.5, 1.0] ./ PERCOLATION_PC_SITE),
    exact=Dict(1 => PERCOLATION_PC_SITE, 2 => 4 / 3, 3 => 5 / 36),
    # Known finite-size bias: corrections to scaling at L = 128 to 512 shift p_c below the exact
    # value for every quality function (Linear, Houdayer and MultipleSplines give 0.59261 to
    # 0.59267) and for autoScale itself. δabs = 1.5 × the measured |p − p*| on the parameters
    # that fail at 3σ. Measured: p_c = 0.5926105 ± 0.0000258 (|dev| 1.36e-4, 5.3σ),
    # ν = 1.3267 ± 0.0042 (1.6σ), β = 0.14078 ± 0.00056 (|dev| 1.89e-3, 3.4σ)
    δabs=Dict(1 => 2.1e-4, 3 => 2.9e-3),
    σmax=Dict(1 => 7.8e-5, 2 => 0.013, 3 => 1.7e-3),
    source="autoScale.py example data (O. Melchert, MIT): 2D site percolation order parameter, L = 128 to 512",
))
