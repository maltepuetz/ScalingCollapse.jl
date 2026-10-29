# CDW transition in the square lattice Holstein model
# (data source: https://owenpb.github.io/FiniteSizeScaling.jl/stable/demo_1/)

include("testdata/cdw_holstein.jl")

@testset "CDW transition (Holstein model)" begin


    sp = ScalingProblem(
        [X_L4, X_L6, X_L8, X_L10, X_L12],
        [Y_L4, Y_L6, Y_L8, Y_L10, Y_L12],
        [E_L4, E_L6, E_L8, E_L10, E_L12],
        [4, 6, 8, 10, 12];
        sf=ScalingCollapse.ScalingFunction(:xny,
            p_names=["x_c", "nu", "gamma"],
            nu=1,
            gamma=1.75
        ),
        p_space=[0.1:0.1:10],
        dx=[-2, 2],
    )
    @test isapprox(sp.optimal_ps[1], 6.0; atol=0.1)
    @test try
        show(sp)
        true
    catch
        false
    end


    sp = ScalingProblem(
        [X_L4, X_L6, X_L8, X_L10, X_L12],
        [Y_L4, Y_L6, Y_L8, Y_L10, Y_L12],
        [E_L4, E_L6, E_L8, E_L10, E_L12],
        [4, 6, 8, 10, 12];
        sf=ScalingCollapse.ScalingFunction(:xny,
            p_names=["x_c", "nu", "gamma"],
            nu=1,
            #gamma=1.75
        ),
        p_space=[5:0.1:7, 1:0.1:3],
        dx=[-2, 2],
    )
    @test isapprox(sp.optimal_ps[1], 6.0; atol=0.2)
    @test isapprox(sp.optimal_ps[2], 1.75; atol=0.1)


    sp = ScalingProblem(
        [X_L4, X_L6, X_L8, X_L10, X_L12],
        [Y_L4, Y_L6, Y_L8, Y_L10, Y_L12],
        [E_L4, E_L6, E_L8, E_L10, E_L12],
        [4, 6, 8, 10, 12];
        sf=ScalingCollapse.ScalingFunction(:xny,
            p_names=["x_c", "nu", "gamma"],
            nu=1,
            #gamma=1.75
        ),
        p_space=[5:0.1:7, 1:0.1:3],
        dx=[-2, 2],
        quality=Houdayer()
    )
    @test isapprox(sp.optimal_ps[1], 6.0; atol=0.2)
    @test isapprox(sp.optimal_ps[2], 1.75; atol=0.1)

    sp = ScalingProblem(
        [X_L4, X_L6, X_L8, X_L10, X_L12],
        [Y_L4, Y_L6, Y_L8, Y_L10, Y_L12],
        [E_L4, E_L6, E_L8, E_L10, E_L12],
        [4, 6, 8, 10, 12];
        sf=ScalingCollapse.ScalingFunction(:xny,
            p_names=["x_c", "nu", "gamma"],
            nu=1,
            #gamma=1.75
        ),
        p_space=[5:0.1:7, 1:0.1:3],
        dx=[-2, 2],
        quality=MultipleSplines()
    )
    @test isapprox(sp.optimal_ps[1], 6.0; atol=0.2)
    @test isapprox(sp.optimal_ps[2], 1.75; atol=0.1)
end
