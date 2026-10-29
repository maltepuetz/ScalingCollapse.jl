# Converts the published data of
#   arXiv:2504.12385 (Zenodo doi:10.5281/zenodo.15227834, LearningTransitions.zip) and
#   arXiv:2608.20268 (Zenodo doi:10.5281/zenodo.21918134, LearningPotts.zip)
# from JLD2 (y-values are `Measurements.Measurement`s) to HDF5 with value and error arrays, keeping
# only the system sizes used by the collapse cases (test/collapses/cases_learning.jl).
#
# Run in a scratch environment with JLD2, Measurements and HDF5:
#   julia --project=<scratch env> convert_papers.jl <LearningTransitions/data> <LearningPotts/data>
using JLD2, Measurements, HDF5
using Measurements: value, uncertainty

lt_dir, potts_dir = ARGS
out = joinpath(@__DIR__, "..", "testdata")

function write_curve(h, group, x, y)
    h["$group/x"] = collect(Float64, x)
    h["$group/y"] = Float64.(value.(y))
    h["$group/e"] = Float64.(uncertainty.(y))
end

loadjld(f, keys...) = jldopen(f) do h
    return Tuple(h[k] for k in keys)
end

h5open(joinpath(out, "learning_ising", "betac.h5"), "w") do h
    # 2504.12385 Fig. 2 inset: I_c(γ) at β_c for L = 8 to 128; layout `Ls`, `L{L}/{x, y, e}`
    gammas, Ics, ds = loadjld(joinpath(lt_dir, "betac.jld2"), "gammas", "Ics", "ds")
    h["Ls"] = ds
    for (k, L) in enumerate(ds)
        write_curve(h, "L$L", gammas[k], Ics[k])
    end
end

h5open(joinpath(out, "learning_ising", "vertical.h5"), "w") do h
    # vertical cuts (fixed γ, scan β), L = ds[2:5]; layout `Ls`, `gammas`, `betas/L{L}` (raw β),
    # `g{i}/L{L}/{y, e}` for cut index i
    Ics, ds, betas, gammas = loadjld(joinpath(lt_dir, "sweep_pd_vertical.jld2"), "Ics", "ds", "betas", "gammas")
    ks = 2:5
    h["Ls"] = ds[ks]
    h["gammas"] = gammas
    for k in ks
        h["betas/L$(ds[k])"] = betas[k]
        for gi in eachindex(gammas)
            y = Ics[k][gi, :]
            h["g$gi/L$(ds[k])/y"] = Float64.(value.(y))
            h["g$gi/L$(ds[k])/e"] = Float64.(uncertainty.(y))
        end
    end
end

h5open(joinpath(out, "learning_ising", "horizontal.h5"), "w") do h
    # horizontal cuts (fixed β, scan γ), L = ds[2:end]; layout `Ls`, `betas`, `b{i}/L{L}/{x, y, e}`
    Ics, ds, betas, gammas = loadjld(joinpath(lt_dir, "sweep_pd_horizontal.jld2"), "Ics", "ds", "betas", "gammas")
    ks = 2:length(ds)
    h["Ls"] = ds[ks]
    h["betas"] = betas
    for bi in eachindex(betas), k in ks
        write_curve(h, "b$bi/L$(ds[k])", gammas[bi][k], Ics[bi][k])
    end
end

# 2608.20268: Potts_Ic_horizontal/3 (Fig. 4 insets, two β), Potts_Ic_horizontal/1 (40 β),
# both with x = γ; layout `Ls`, `betas`, `b{i}/L{L}/{x, y, e}`
for (src, dst) in (("horizontal/3", "insets"), ("horizontal/1", "horizontal"))
    Ics, γs, ds, βs = loadjld(joinpath(potts_dir, "Potts_Ic_$src", "averaged_results.jld2"), "Ics", "γs", "ds", "βs")
    h5open(joinpath(out, "learning_potts", dst * ".h5"), "w") do h
        h["Ls"] = ds
        h["betas"] = βs
        for bi in eachindex(βs), (k, L) in enumerate(ds)
            write_curve(h, "b$bi/L$L", γs[k][bi], Ics[k, bi])
        end
    end
end

# Potts_Ic_vertical/4 (Fig. 10): scan β at fixed γ; layout `Ls`, `gammas`, `g{i}/L{L}/{x, y, e}`
let
    Ics, γs, ds, βs = loadjld(joinpath(potts_dir, "Potts_Ic_vertical/4", "averaged_results.jld2"), "Ics", "γs", "ds", "βs")
    h5open(joinpath(out, "learning_potts", "vertical.h5"), "w") do h
        h["Ls"] = ds
        h["gammas"] = γs
        for gi in eachindex(γs), (k, L) in enumerate(ds)
            write_curve(h, "g$gi/L$L", βs[k][gi], Ics[k, gi])
        end
    end
end
