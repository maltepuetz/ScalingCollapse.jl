# CDW transition in the square lattice Holstein model (data: ../testdata/cdw_holstein.jl)
include(joinpath(@__DIR__, "..", "testdata", "cdw_holstein.jl"))

# R6: the first ScalingProblem call of ../cdw.jl; dx is the window [-2, 2] of the relative form
# (x − x_c)/x_c·L^{1/ν}, times x_c ≈ 6
push!(CASES, CollapseCase(;
    name="cdw_holstein", kind=:regression,
    load=() -> ([X_L4, X_L6, X_L8, X_L10, X_L12],
        [Y_L4, Y_L6, Y_L8, Y_L10, Y_L12],
        [E_L4, E_L6, E_L8, E_L10, E_L12],
        [4, 6, 8, 10, 12]),
    kwargs=(
        sf=ScalingCollapse.ScalingFunction(:xny, p_names=["x_c", "nu", "gamma"], nu=1, gamma=1.75),
        p_space=[0.1:0.1:10],
        dx=[-12, 12],
    ),
    source="https://owenpb.github.io/FiniteSizeScaling.jl/stable/demo_1/",
))
