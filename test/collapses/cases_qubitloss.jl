# Coherent information in weakly monitored quantum circuits with qubit loss
# (arXiv:2505.22720, Fig. 3; data source and license: ../testdata/README.md)
using HDF5

function load_qubitloss(theta_index, Ls)
    h5open(joinpath(@__DIR__, "..", "testdata", "qubitloss", "qubitloss_Ic.h5"), "r") do file
        xs = [read(file, "theta$theta_index/L$L/r") for L in Ls]
        ys = [read(file, "theta$theta_index/L$L/Ic") for L in Ls]
        es = [read(file, "theta$theta_index/L$L/Ic_err") for L in Ls]
        return xs, ys, es, Ls
    end
end

# the ScalingProblem call used in the paper
const QUBITLOSS_KWARGS = (
    sf=ScalingFunction(:x),
    dx=[-1.25, 1.25],
    p_space=[0.3:0.001:0.6, 1.2:0.1:1.8],
    error_threshold=1.3,
    quality=MultipleSplines(),
)

# E1: θ = π/2 is bond percolation, r_c = 1/2 and ν = 4/3
push!(CASES, CollapseCase(;
    name="qubitloss_theta_pi2_exact", kind=:exact,
    load=() -> load_qubitloss(15, [32, 64, 128, 256]),
    kwargs=QUBITLOSS_KWARGS,
    exact=Dict(1 => 0.5, 2 => 4 / 3),
    # ceilings at ~3× the error of a converged fit (4.33e-5, 5.92e-3); a scan grid point
    # returned without Nelder-Mead progress has 1.9e-4 and 8.8e-2
    σmax=Dict(1 => 1.3e-4, 2 => 1.8e-2),
    source="arXiv:2505.22720, Fig. 3, θ = π/2 (bond percolation)",
))

# R1: all 15 angles, regression against the published ν (and r_c for θ = 0)
for i in 1:15
    Ls = i == 1 ? [32, 64, 128, 256, 512] : [32, 64, 128, 256]
    push!(CASES, CollapseCase(;
        name="qubitloss_theta$i", kind=:regression,
        load=() -> load_qubitloss(i, Ls),
        kwargs=QUBITLOSS_KWARGS,
        source="arXiv:2505.22720, Fig. 3, angle index $i",
    ))
end
