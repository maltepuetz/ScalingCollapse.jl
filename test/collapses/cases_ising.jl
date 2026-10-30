# 2D Ising Monte Carlo data (testdata/README.md)
using HDF5

function load_ising_binder()
    h5open(joinpath(@__DIR__, "..", "testdata", "IsingMonteCarlo.h5"), "r") do file
        Ts, Ls, binder = read(file, "Ts"), read(file, "Ls"), read(file, "binder")
        return [Ts for _ in eachindex(Ls)], [binder[:, i] for i in eachindex(Ls)], nothing, Ls
    end
end

# E6: exact 2D Ising T_c = 2/log(1+√2) and ν = 1
push!(CASES, CollapseCase(;
    name="ising_mc_binder_exact", kind=:exact,
    load=load_ising_binder,
    kwargs=(sf=ScalingFunction(:x), p_space=[0.1:0.1:3, 0.1:0.1:3], dx=[-1.0, 1.0]),
    exact=Dict(1 => 2 / log(1 + √2), 2 => 1.0),
    # measured: T_c = 2.26814 ± 0.00239 (|dev| 1.0e-3), ν = 0.9555 ± 0.0186 (|dev| 4.5e-2), both
    # within 3σ, so δabs = 0; σmax ≈ 3× the fitted errors
    σmax=Dict(1 => 7.2e-3, 2 => 5.6e-2),
    source="2D Ising Monte Carlo, Binder cumulant",
))
