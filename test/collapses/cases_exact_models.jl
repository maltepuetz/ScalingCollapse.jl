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
    # measured: T_c = 2.269210 ± 2.05e-4 (0.12σ), ν = 1.00764 ± 5.58e-3 (1.4σ)
    σmax=Dict(1 => 6.2e-4, 2 => 1.7e-2),
    source="exact periodic Ising strip correlation length, L = 8 to 64",
))

# E7: finite-N Curie–Weiss model with N in place of L, so x scales as N^{1/2} (p2 = 2), and
# m_rms ∝ N^{-1/4} (p3/p2 = 1/4). NOT in CASES: both fail the exact criterion at 3σ with
# δabs = 0. The generator is verified (brute-force enumeration to 1e-12, N → ∞ limits of U and
# ⟨M²⟩/N^{3/2}); the misses are N^{-1/2} corrections to scaling that the quality-ratio error
# estimate does not cover on noise-free data. Measured with dx = [-2, 2]:
#   binder: T_c = 0.98840 ± 0.00304 (3.8σ), ν = 1.906 ± 0.127 (0.7σ)
#   m_rms:  T_c = 0.99756 ± 0.00050 (4.9σ), ν = 1.847 ± 0.015 (10σ), p3 = 0.4484 ± 0.0037 (14σ)
# Narrower or wider dx does not help. To be enabled as broken cases once the framework can
# mark them (Task 3); σmax is 3× the measured errors.
const CURIE_WEISS_EXACT_CASES = [
    CollapseCase(;
        name="curie_weiss_binder_exact", kind=:exact,
        load=() -> load_exact("curie_weiss.h5", "binder"),
        kwargs=(sf=ScalingFunction(:x), p_space=[0.8:0.02:1.2, 1.0:0.1:3.0], dx=[-2.0, 2.0]),
        exact=Dict(1 => 1.0, 2 => 2.0),
        σmax=Dict(1 => 9.1e-3, 2 => 3.8e-1),
        source="exact finite-N Curie–Weiss Binder ratio, N = 64 to 1024",
    ),
    CollapseCase(;
        name="curie_weiss_mrms_exact", kind=:exact,
        load=() -> load_exact("curie_weiss.h5", "m_rms"),
        kwargs=(sf=ScalingFunction(:xy), p_space=[0.8:0.02:1.2, 1.0:0.1:3.0, 0.1:0.1:1.0], dx=[-2.0, 2.0]),
        exact=Dict(1 => 1.0, 2 => 2.0, 3 => 0.5),
        σmax=Dict(1 => 1.5e-3, 2 => 4.6e-2, 3 => 1.1e-2),
        source="exact finite-N Curie–Weiss m_rms, N = 64 to 1024",
    ),
]
