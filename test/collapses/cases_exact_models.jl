# Exact, noise-free model data from test/generators/ising_strip.jl and curie_weiss.jl
using HDF5

# curves `prefix/L{L}/{x, y}` of test/testdata/exact/<file>
function load_exact(file, prefix="")
    h5open(joinpath(@__DIR__, "..", "testdata", "exact", file), "r") do h
        Ls = read(h, "Ls")
        p = isempty(prefix) ? "" : prefix * "/"
        return [read(h, "$(p)L$L/x") for L in Ls], [read(h, "$(p)L$L/y") for L in Ls], nothing, Ls
    end
end

# E5: ξ_L/L of the periodic Ising strip (Kaufman), T_c = 2/log(1+√2) and ν = 1. dx covers the
# whole overlap of the L = 8 and L = 64 curves.
push!(CASES, CollapseCase(;
    name="ising_strip_xi_exact", kind=:exact,
    load=() -> load_exact("ising_strip.h5"),
    kwargs=(sf=ScalingFunction(:x), p_space=[2.0:0.02:2.5, 0.5:0.1:1.6], dx=[-1.0, 1.0]),
    exact=Dict(1 => 2 / log(1 + √2), 2 => 1.0),
    # measured: T_c = 2.269210 ± 2.34e-4 (0.11σ), ν = 1.00764 ± 6.56e-3 (1.2σ)
    σmax=Dict(1 => 7.0e-4, 2 => 2.0e-2),
    source="exact periodic Ising strip correlation length, L = 8 to 64",
))

# E7: finite-N Curie–Weiss model with N in place of L, so x scales as N^{1/2} (p2 = 2), and
# m_rms ∝ N^{-1/4} (p3/p2 = 1/4). Known finite-size bias: N^{-1/2} corrections to scaling, which
# the quality-ratio error does not cover on noise-free data. The generator is checked against
# brute-force enumeration over all 2^N states (to 1e-12) and against the N → ∞ limits of U and
# ⟨M²⟩/N^{3/2}. δabs = 1.5 × the measured |p − p*| on the parameters that fail at 3σ, 0 on the
# others; σmax ≈ 3× the measured errors.
push!(CASES, CollapseCase(;
    name="curie_weiss_binder_exact", kind=:exact,
    load=() -> load_exact("curie_weiss.h5", "binder"),
    kwargs=(sf=ScalingFunction(:x), p_space=[0.8:0.02:1.2, 1.0:0.1:3.0], dx=[-2.0, 2.0]),
    exact=Dict(1 => 1.0, 2 => 2.0),
    # measured: T_c = 0.98840 ± 0.00329 (|dev| 1.16e-2, 3.5σ), ν = 1.906 ± 0.131 (|dev| 9.4e-2,
    # 0.7σ, so no δabs)
    δabs=Dict(1 => 0.018),
    σmax=Dict(1 => 9.9e-3, 2 => 3.9e-1),
    source="exact finite-N Curie–Weiss Binder ratio, N = 64 to 1024",
))
push!(CASES, CollapseCase(;
    name="curie_weiss_mrms_exact", kind=:exact,
    load=() -> load_exact("curie_weiss.h5", "m_rms"),
    kwargs=(sf=ScalingFunction(:xy), p_space=[0.8:0.02:1.2, 1.0:0.1:3.0, 0.1:0.1:1.0], dx=[-2.0, 2.0]),
    exact=Dict(1 => 1.0, 2 => 2.0, 3 => 0.5),
    # measured: T_c = 0.99756 ± 0.00362 (|dev| 2.44e-3, 0.7σ), ν = 1.847 ± 0.033 (|dev| 0.153,
    # 4.7σ), p3 = 0.4484 ± 0.0219 (|dev| 5.16e-2, 2.4σ). δabs was set from the conditional
    # errors (T_c 4.9σ, p3 14σ); with profile errors only ν needs it, the others are kept.
    δabs=Dict(1 => 0.0037, 2 => 0.24, 3 => 0.078),
    σmax=Dict(1 => 1.1e-2, 2 => 9.9e-2, 3 => 6.6e-2),
    source="exact finite-N Curie–Weiss m_rms, N = 64 to 1024",
))
