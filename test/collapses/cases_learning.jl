# Learning transitions of classical Ising models and deformed toric codes (arXiv:2504.12385,
# Zenodo doi:10.5281/zenodo.15227834) and of Potts models and Z_3 toric codes (arXiv:2608.20268,
# Zenodo doi:10.5281/zenodo.21918134); data source, license and layout: ../testdata/README.md.
# The ScalingProblem calls are those of the papers. `[published]` in references.toml holds the
# stored fits of the deposits (full precision) and the numbers in the paper text.
using HDF5

const LEARNING_DATA = joinpath(@__DIR__, "..", "testdata")
const LEARNING_ISING_βC = log(1 + sqrt(2)) / 2  # critical β of the 2D Ising model

# per-size curves `<group>/L{L}/{x, y, e}`, with the x values of `xof(L)` when given
function load_learning_curves(file, group; xof=nothing)
    h5open(joinpath(LEARNING_DATA, file), "r") do h
        Ls = read(h, "Ls")
        xs = [xof === nothing ? read(h, "$group/L$L/x") : xof(h, L) for L in Ls]
        ys = [read(h, "$group/L$L/y") for L in Ls]
        es = [read(h, "$group/L$L/e") for L in Ls]
        return xs, ys, es, Ls
    end
end

# vertical cuts of the Ising phase diagram: x = β/β_c (raw β stored once per size)
load_learning_ising_vertical(gi) =
    load_learning_curves("learning_ising/vertical.h5", "g$gi"; xof=(h, L) -> read(h, "betas/L$L") ./ LEARNING_ISING_βC)

const LEARNING_ISING_VERTICAL_KWARGS = (
    sf=ScalingFunction(:xy),
    p_space=[0.6:0.01:1.3, 0.5:0.1:1.5, -0.5:0.1:0.5],
    quality=MultipleSplines(N_steps=100),
    dx=[-1.25, 1.75],
)

# E2: the vertical cut γ = 0.3 (index 16); x = β/β_c, so the 2D Ising critical point is x_c = 1,
# with ν = 1.
push!(CASES, CollapseCase(;
    name="learning_ising_vertical_g16_exact", kind=:exact,
    load=() -> load_learning_ising_vertical(16),
    kwargs=LEARNING_ISING_VERTICAL_KWARGS,
    exact=Dict(1 => 1.0, 2 => 1.0),
    # published fit 1.00058(16) is 3.6σ from 1 (measured |x_c − 1| = 5.75e-4); small-size bias
    δabs=Dict(1 => 1e-3),
    # measured: x_c = 1.000575 ± 0.000340 (1.7σ), ν = 0.9981 ± 0.0079 (0.2σ). δabs was set from
    # the conditional errors and is kept, although with profile errors x_c is now 1.7σ off.
    σmax=Dict(1 => 1.0e-3, 2 => 2.4e-2),
    source="arXiv:2504.12385, vertical cut γ = 0.3 (index 16), L = ds[2:5]",
))

# R2: Fig. 2 inset, L = 8 to 128
push!(CASES, CollapseCase(;
    name="learning_ising_betac", kind=:regression,
    load=() -> load_learning_curves("learning_ising/betac.h5", ""),
    kwargs=(sf=ScalingFunction(:xy), dx=[-2, 2],
        p_space=[0.3:0.1:0.9, 1.0:0.1:5.0, 0.0:0.01:0.1], quality=MultipleSplines()),
    source="arXiv:2504.12385, Fig. 2 inset (β_c), L = 8 to 128",
))

# R3: 31 vertical and 41 horizontal cuts
for gi in 1:31
    push!(CASES, CollapseCase(;
        name="learning_ising_vertical_g$gi", kind=:regression,
        load=() -> load_learning_ising_vertical(gi),
        kwargs=LEARNING_ISING_VERTICAL_KWARGS,
        source="arXiv:2504.12385, vertical cut $gi of 31, L = ds[2:5]",
    ))
end
for bi in 1:41
    push!(CASES, CollapseCase(;
        name="learning_ising_horizontal_b$bi", kind=:regression,
        load=() -> load_learning_curves("learning_ising/horizontal.h5", "b$bi"),
        kwargs=(sf=ScalingFunction(:xy), dx=[-1.25, 1.25],
            p_space=[0.3:0.01:1.0, 1.0:0.1:2.5, -0.4:0.1:0.4], quality=MultipleSplines(N_steps=100)),
        source="arXiv:2504.12385, horizontal cut $bi of 41, L = ds[2:end]",
    ))
end

# R4: Fig. 4 insets, N⁽¹⁾ (cut 1) and N⁽²⁾ (cut 2) of the Potts model
for (name, bi) in (("learning_potts_inset_N1", 1), ("learning_potts_inset_N2", 2))
    push!(CASES, CollapseCase(;
        name=name, kind=:regression,
        load=() -> load_learning_curves("learning_potts/insets.h5", "b$bi"),
        kwargs=(sf=ScalingFunction(:xy), dx=bi == 2 ? [-2, 2] : [-1.5, 1.5],
            p_space=[0.4:0.01:1, 0.1:0.1:3, -0.5:0.1:0.5], quality=MultipleSplines()),
        source="arXiv:2608.20268, Fig. 4 inset, N($bi)",
    ))
end

# R5: 40 horizontal and 30 vertical Potts cuts
for bi in 1:40
    # b40: no [published] entry, the stored fit (v0.2.5 + Optim 1.13.3) was hit by the Optim ≥ 1.12
    # simplex bug: quality 2.032 vs 1.866 at the true minimum; γ_c differs by 1.19σ. Checked
    # against [last] only.
    push!(CASES, CollapseCase(;
        name="learning_potts_horizontal_b$bi", kind=:regression,
        load=() -> load_learning_curves("learning_potts/horizontal.h5", "b$bi"),
        kwargs=(sf=ScalingFunction(:xy), dx=[-1.5, 1.5],
            p_space=[0.5:0.01:1, 0.1:0.1:3, -0.5:0.1:0.5], quality=MultipleSplines()),
        source="arXiv:2608.20268, horizontal cut $bi of 40 (Potts_Ic_horizontal/1)",
    ))
end
for gi in 1:30
    # g28: no [published] entry. In the profile error analysis a lower quality was found: 2.021
    # vs 2.646 at the published optimum, so the published fit was not optimal (this does not
    # show that the new fit is the global minimum). The re-optimized fit moves by 0.92σ, -0.62σ
    # and 0.85σ (published σ) in x_c, ν and p3. Checked against [last] only.
    push!(CASES, CollapseCase(;
        name="learning_potts_vertical_g$gi", kind=:regression,
        load=() -> load_learning_curves("learning_potts/vertical.h5", "g$gi"),
        kwargs=(sf=ScalingFunction(:xy), dx=[-1.5, 2.5],
            p_space=[0.99:0.001:1.01, 0.6:0.1:1, -0.5:0.1:0.5], quality=MultipleSplines(),
            optimizer=ScalingCollapse.Optim.NelderMead(;
                initial_simplex=ScalingCollapse.Optim.AffineSimplexer(; a=1e-6, b=1e-6))),
        source="arXiv:2608.20268, Fig. 10, vertical cut $gi of 30 (Potts_Ic_vertical/4)",
    ))
end
