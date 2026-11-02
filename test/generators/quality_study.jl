# Seven noisy synthetic data sets from the quality-function study (PR #63) that motivated pchip
# interpolation in MultipleSplines (irregular x grids, few points per size, narrow peaks, several
# extrema). Writes test/testdata/synthetic/quality_study.h5, one group per data set. Not run in CI.
#
# Model: y = L^a F(z) + noise with z = (x − x_c)/x_c · L^{1/ν} and a = 0 (:x), −p3/ν (:xy) or
# +p3/ν (:xny). The z values of every size span [−2.5, 2.5]: `npts` equidistant points shifted by
# a random offset of up to half a spacing per size (so that sizes never share scaled x values at
# the true parameters), or `npts` uniformly random points. Noise is Gaussian with σ = η|y|, and
# the error bars are the exact σ. Sizes L = 16, 32, 64, 128.
#
# Every data set is replicate 1 of the study's configuration of the same name, drawn from
# Xoshiro(hash((1, study_name, 1))). `hash` is not guaranteed to be stable across Julia versions,
# so the stored file (written with Julia 1.13) is the reference.
#
# Layout: `<name>/Ls`, `<name>/L{L}/{x, y, e}`; the group's attributes `preset`, `study_name`,
# `truth` (x_c, ν[, p3]), `npts`, `grid` and `noise` describe the data set.
#
# Usage (in an environment with HDF5):
#     julia quality_study.jl [out.h5]

using HDF5
using Random

# name => (study_name, preset, truth, F, npts, grid, noise)
const DATASETS = [
    "random_crossing" => ("grid/crossing/random", :x, [0.7, 1.2], z -> 1 + tanh(z), 30, :random, 0.01),
    "random_order_12" => ("grid/order_2d/random,n=12", :xy, [2.27, 1.0, 0.125],
        z -> (1 - tanh(z)) / 2 + 0.05, 12, :random, 0.01),
    "random_peak" => ("grid/susc_gauss/random", :xny, [2.27, 1.0, 1.75], z -> exp(-z^2) + 0.1, 30, :random, 0.01),
    "few_points_crossing" => ("npts/crossing/n=6", :x, [0.7, 1.2], z -> 1 + tanh(z), 6, :uniform, 0.01),
    "few_points_peak" => ("npts/susc_gauss/n=6", :xny, [2.27, 1.0, 1.75], z -> exp(-z^2) + 0.1, 6, :uniform, 0.01),
    "narrow_peak" => ("peaks/susc_sharp/n=20/η=0.01", :xny, [2.27, 1.0, 1.75],
        z -> exp(-z^2 / 0.08) + 0.1, 20, :uniform, 0.01),
    "wiggly" => ("peaks/wiggly/n=20/η=0.01", :x, [0.7, 1.2],
        z -> 1 + tanh(z) + 0.3sin(4z) * exp(-z^2 / 2), 20, :uniform, 0.01),
]
const LS = [16, 32, 64, 128]
const ZMAX = 2.5

function generate(study_name, preset, truth, F, npts, grid, noise)
    rng = Xoshiro(hash((1, study_name, 1)))
    xc, ν = truth[1:2]
    a = preset === :x ? 0.0 : (preset === :xy ? -1 : 1) * truth[3] / ν
    xs, ys, es = Vector{Float64}[], Vector{Float64}[], Vector{Float64}[]
    for L in LS
        Δ = 2ZMAX / (npts - 1)
        zs = grid === :uniform ? collect(range(-ZMAX, ZMAX; length=npts)) .+ (rand(rng) - 0.5) * Δ :
             sort!(-ZMAX .+ 2ZMAX .* rand(rng, npts))
        x = @. xc * (1 + zs * L^(-1 / ν))
        ztrue = @. (x - xc) / xc * L^(1 / ν)
        y = @. L^a * F(ztrue)
        σ = noise .* abs.(y)
        push!(xs, x)
        push!(ys, y .+ σ .* randn(rng, length(y)))
        push!(es, σ)
    end
    return xs, ys, es
end

function write_datasets(path)
    h5open(path, "w") do h
        for (name, (study_name, preset, truth, F, npts, grid, noise)) in DATASETS
            xs, ys, es = generate(study_name, preset, truth, F, npts, grid, noise)
            g = create_group(h, name)
            g["Ls"] = LS
            attrs(g)["preset"] = String(preset)
            attrs(g)["study_name"] = study_name
            attrs(g)["truth"] = truth
            attrs(g)["npts"] = npts
            attrs(g)["grid"] = String(grid)
            attrs(g)["noise"] = noise
            for (i, L) in enumerate(LS)
                g["L$L/x"] = xs[i]
                g["L$L/y"] = ys[i]
                g["L$L/e"] = es[i]
            end
        end
    end
end

write_datasets(isempty(ARGS) ? joinpath(@__DIR__, "..", "testdata", "synthetic", "quality_study.h5") : ARGS[1])
