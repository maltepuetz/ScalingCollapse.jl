# Reduce the raw qubit-loss data of arXiv:2505.22720 to the averaged coherent information
# curves in test/testdata/qubitloss/qubitloss_Ic.h5. Not run in CI.
#
# Input: `sweep_pd.jld2` and `nishimori_line.jld2` from the `data` folder of
# NishimoriQubitLoss.zip (Zenodo, doi:10.5281/zenodo.15530142, CC-BY-4.0).
# The reduction follows code/fig3.ipynb of that repository: I_c is computed per sample from
# <ZZ> and averaged, with the standard error of the mean as uncertainty.
#
# Usage (in an environment with JLD2 and HDF5):
#     julia qubitloss_reduce.jl <data dir> <out.h5>

using JLD2, HDF5, Statistics

function Ic(kappa)
    res = -(1 + kappa) / 2 * log2(abs(1 + kappa) / 2) - (1 - kappa) / 2 * log2(abs(1 - kappa) / 2)
    if isnan(res) || isinf(res)
        b = BigFloat(kappa)
        res = -(1 + b) / 2 * log2(abs(1 + b) / 2) - (1 - b) / 2 * log2(abs(1 - b) / 2)
    end
    isnan(res) && (res = 0.0)
    return Float64(res)
end

# mean, standard error and number of samples of I_c for one parameter point
reduce_samples(zz) = (ic = Ic.(zz); (mean(ic), std(ic) / sqrt(length(ic)), length(ic)))

function write_curve(h, grp, rs, zzs)
    r = reduce_samples.(zzs)
    h["$grp/r"] = collect(rs)
    h["$grp/Ic"] = first.(r)
    h["$grp/Ic_err"] = getindex.(r, 2)
    h["$grp/nsamples"] = last.(r)
end

data, out = ARGS

f = jldopen(joinpath(data, "sweep_pd.jld2"), "r")
Ls, thetas, rs, zz = f["Ls"], f["thetas"], f["rs"], f["expect_ZZs"]
close(f)

# additional system size L = 512 on the Nishimori line (theta = 0)
g = jldopen(joinpath(data, "nishimori_line.jld2"), "r")
Ls_n, rs_n, zz_n = g["Ls"], g["rs"], g["expect_ZZs"]
close(g)

h5open(out, "w") do h
    h["Ls"] = collect(Ls)
    h["thetas"] = collect(thetas)
    for ti in eachindex(thetas), li in eachindex(Ls)
        write_curve(h, "theta$(ti)/L$(Ls[li])", rs[ti, li], zz[ti, li])
    end
    write_curve(h, "theta1/L$(Ls_n[1])", rs_n[1], zz_n[1])
end
