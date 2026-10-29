# Synthetic finite-size-scaling data with known parameters. Writes
# test/testdata/synthetic/synthetic.h5, one group per variant. Not run in CI.
#
# Model: y = L^{s p3/p2} F((x − x_c)/x_c L^{1/p2}) (1 + c L^{−ω}) + noise, with x_c = 0.7,
# p2 = ν = 1.2, p3 = 0.3; s = −1 for :xy, +1 for :xny, and no prefactor for :x. F is the monotone
# 1 + tanh(z) or the peaked exp(−z²) + 0.1. Defaults: L = 16, 32, 64, 128, 41 equidistant
# x ∈ [0.55, 0.85], c = 0, no noise. Noise is Gaussian with σ = 1% of |y|.
#
# Layout: `<variant>/Ls`, `<variant>/L{L}/{x, y}` and, where error bars exist, `<variant>/L{L}/e`;
# the variant's attributes `preset`, `xc`, `nu`, `p3` hold the true values.
#
# Usage (in an environment with HDF5):
#     julia synthetic.jl [out.h5]

using HDF5
using Random

const XC, NU, P3 = 0.7, 1.2, 0.3
monotone(z) = 1 + tanh(z)
peaked(z) = exp(-z^2) + 0.1

function curve(L, xs, preset, F; c=0.0, ω=1.0)
    s = preset === :xy ? -1 : preset === :xny ? 1 : 0
    return [L^(s * P3 / NU) * F((x - XC) / XC * L^(1 / NU)) * (1 + c * L^(-ω)) for x in xs]
end

function variants(rng)
    Ls = [16, 32, 64, 128]
    grid = collect(range(0.55, 0.85; length=41))
    noisy(ys) = ys .+ 0.01 .* abs.(ys) .* randn(rng, length(ys))
    out = Dict{String,Any}()
    clean(preset, F, Ls=Ls) = (preset=preset, Ls=Ls, xs=[grid for _ in Ls], ys=[curve(L, grid, preset, F) for L in Ls], es=nothing)
    out["clean_x"] = clean(:x, monotone)
    out["clean_xy"] = clean(:xy, monotone)
    out["clean_xny"] = clean(:xny, peaked)

    base = clean(:xy, monotone)
    ys = [noisy(y) for y in base.ys]
    out["noise_err"] = (; base..., ys, es=[0.01 .* abs.(y) for y in ys])
    out["noise_noerr"] = (; base..., ys, es=nothing)
    out["three_sizes"] = clean(:xy, monotone, [16, 32, 64])
    xs = [sort(0.55 .+ 0.3 .* rand(rng, 41)) for _ in Ls]
    out["uneven_grid"] = (; base..., xs, ys=[curve(L, x, :xy, monotone) for (L, x) in zip(Ls, xs)])
    ys = [noisy(curve(L, grid, :xy, monotone; c=0.5, ω=1.0)) for L in Ls]
    out["corrections"] = (; base..., ys, es=[0.01 .* abs.(y) for y in ys])
    out["fixed_nu"] = base
    return out
end

function write_variants(path)
    h5open(path, "w") do h
        for (name, v) in sort!(collect(variants(Random.Xoshiro(20261029))); by=first)
            g = create_group(h, name)
            g["Ls"] = collect(v.Ls)
            attrs(g)["preset"] = String(v.preset)
            attrs(g)["xc"] = XC
            attrs(g)["nu"] = NU
            v.preset === :x || (attrs(g)["p3"] = P3)
            for (i, L) in enumerate(v.Ls)
                g["L$L/x"] = v.xs[i]
                g["L$L/y"] = v.ys[i]
                v.es === nothing || (g["L$L/e"] = v.es[i])
            end
        end
    end
end

write_variants(isempty(ARGS) ? joinpath(@__DIR__, "..", "testdata", "synthetic", "synthetic.h5") : ARGS[1])
