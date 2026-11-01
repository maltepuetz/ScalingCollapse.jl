@testset "Error handling" begin
    Optim = ScalingCollapse.Optim

    @testset "unknown keyword arguments" begin
        # typo `sc` instead of `sf`
        @test_throws ArgumentError ScalingProblem(Ts, binder, Ls; sc=ScalingFunction(:x))
        @test_throws ArgumentError ScalingProblem(Ts, binder, Ls;
            sf=ScalingFunction(:x), d_x=[-1.0, 1.0],
        )
        # ScalingFunction kwargs must go to ScalingFunction when `sf` is given
        @test_throws ArgumentError ScalingProblem(Ts, binder, Ls;
            sf=ScalingFunction(:x), p_names=["T_c", "nu"],
        )
        # fixing a parameter that does not exist
        @test_throws ArgumentError ScalingFunction(:x; nu=1.0)
        @test_throws ArgumentError ScalingFunction(:x; p_names=["T_c", "nu"], p2=1.0)
        @test_throws ArgumentError ScalingFunction(["a", "b", "c", "d"])

        # ScalingFunction kwargs are still forwarded to the default ScalingFunction
        sp = ScalingProblem(Ts, binder, Ls;
            p_names=["T_c", "nu"],
            dx=[-2.0, 2.0],
            starting_ps=[2.27, 1.0],
            error=false,
        )
        @test ScalingCollapse.scaled_p_names(sp.sf) == ["T_c", "nu"]
        @test isapprox(sp.optimal_ps[1], 2.269; atol=0.01)
    end

    @testset "invalid options" begin
        sf = ScalingFunction(:x)
        @test_throws ArgumentError ScalingProblem(Ts, binder, Ls; sf=sf, p_space=[0.1:0.1:3.0])
        @test_throws ArgumentError ScalingProblem(Ts, binder, Ls; sf=sf, starting_ps=[2.27])
        @test_throws ArgumentError ScalingProblem(Ts, binder, Ls; sf=sf, dx=[1.0, -1.0])
        @test_throws ArgumentError ScalingProblem(Ts, binder, Ls; sf=sf, dx=[1.0])
        @test_throws ArgumentError ScalingProblem(Ts, binder, Ls; sf=sf, error_threshold=1.0)
        @test_throws ArgumentError ScalingProblem(Ts, binder, Ls;
            sf=sf, error_threshold=[1.3, 1.3, 1.3],
        )
        @test_throws ArgumentError ScalingProblem(Ts, binder, Ls; sf=:x)
        @test_throws ArgumentError ScalingProblem(Ts, binder, Ls;
            sf=sf, p_space=[2.5:-0.1:2.0, 0.5:0.1:1.5],
        )
    end

    @testset "invalid data" begin
        sf = ScalingFunction(:x)
        @test_throws ArgumentError ScalingProblem(Ts, binder, Float64.(Ls); sf=sf)
        @test_throws ArgumentError ScalingProblem(Ts, binder[:, 1:1], Ls[1:1]; sf=sf)
        @test_throws ArgumentError ScalingProblem(Ts[1:end-1], binder, Ls; sf=sf)
        nan_binder = copy(binder)
        nan_binder[3, 2] = NaN
        @test_throws ArgumentError ScalingProblem(Ts, nan_binder, Ls; sf=sf)
        # repeated x value for one system size
        @test_throws ArgumentError ScalingProblem(
            vcat(Ts, Ts[end]), vcat(binder, binder[end:end, :]), Ls; sf=sf,
        )
        # MultipleSplines needs at least 4 points per system size
        @test_throws ArgumentError ScalingProblem(Ts[1:3], binder[1:3, :], Ls; sf=sf)
    end

    @testset "invalid custom scaling functions" begin
        p_space = [2.0:0.1:2.4, 0.5:0.1:1.5]
        # takes one parameter, but two are declared
        @test_throws ArgumentError ScalingProblem(Ts, binder, Ls;
            sf=ScalingFunction((d, a) -> d; p_names=["a", "b"]), p_space=p_space,
        )
        # does not return a Data object
        @test_throws ArgumentError ScalingProblem(Ts, binder, Ls;
            sf=ScalingFunction((d, a, b) -> (d.xs, d.ys); p_names=["a", "b"]),
            p_space=p_space,
        )
    end

    @testset "residuals without valid errors" begin
        sp = ScalingProblem(Ts, binder, Ls;
            sf=ScalingFunction(:x), dx=[-2.0, 2.0], starting_ps=[2.27, 1.0], error=false,
        )
        @test_throws ArgumentError residuals(sp; N_steps=3)
        sp.optimal_ps_error .= NaN
        @test_throws ArgumentError residuals(sp; N_steps=3)
        # an explicit p_space does not need the errors
        p_space, res = residuals(sp; p_space=[2.2:0.05:2.3, 0.9:0.05:1.0], dims=[1])
        @test size(res) == (3,)
        # only the errors of the active dimensions are needed
        sp.optimal_ps_error .= [0.01, NaN]
        p_space, res = residuals(sp; N_steps=3, dims=[1])
        @test size(res) == (3,)
    end

    @testset "no overlap in the optimization interval" begin
        sf = ScalingFunction(:x)
        @test_throws ErrorException ScalingProblem(Ts, binder, Ls;
            sf=sf, dx=[100.0, 101.0], p_space=[2.0:0.1:2.4, 0.5:0.1:1.5],
        )
        @test_throws ErrorException ScalingProblem(Ts, binder, Ls;
            sf=sf, dx=[100.0, 101.0], starting_ps=[2.27, 1.0],
        )
        # too few points for the Linear quality function
        @test_throws ErrorException ScalingProblem(Ts, binder, Ls;
            sf=sf, dx=[0.0, 1e-9], starting_ps=[2.27, 1.0],
        )
    end

    @testset "non-converged optimization warns" begin
        @test_logs (:warn, r"without converging") match_mode = :any ScalingProblem(
            Ts, binder, Ls;
            sf=ScalingFunction(:x),
            starting_ps=[2.27, 1.0],
            error=false,
            optimizer_options=Optim.Options(iterations=2),
        )
    end

    @testset "parameter scan edge cases" begin
        # quality does not depend on the second parameter -> no isolated minimum
        sp = OptimizationConfigurationHarness(
            false,
            [[-1.0, 0.0, 1.0], [-1.0, 0.0, 1.0]],
            (sp, ps; check_bounds) -> ps[1]^2,
            (sp, ps; check_bounds) -> ps[1]^2,
            zeros(2),
            Inf,
            zeros(2),
            Optim.NelderMead(),
            Optim.Options(),
        )
        err = @test_throws ErrorException ScalingCollapse.parameter_scan!(sp)
        @test occursin("along p2.", err.value.msg)

        # same for a ScalingProblem with a parameter that does not affect the scaling
        sf = ScalingFunction(
            (d, T_c, dummy) -> ScalingCollapse._power_scaling(d, T_c, 1.0);
            p_names=["T_c", "dummy"],
        )
        err = @test_throws ErrorException ScalingProblem(Ts, binder, Ls;
            sf=sf, dx=[-2.0, 2.0], p_space=[2.0:0.1:2.4, 0.5:0.1:1.5],
        )
        @test occursin("along dummy.", err.value.msg)

        # quality is infinite (or NaN) everywhere
        sp.quality_scan = (sp, ps; check_bounds) -> NaN
        @test_throws ErrorException ScalingCollapse.parameter_scan!(sp)
    end

    @testset "error analysis" begin
        harness(quality, optimal_ps, minimum) = (
            sf=ScalingFunction(:x; p_names=["a", "b"]),
            quality=quality,
            optimal_ps=optimal_ps,
            optimal_ps_error=zeros(length(optimal_ps)),
            minimum=minimum,
            error_threshold=[1.3, 1.3],
            verbose=false,
            error=:conditional,
        )

        # optimum at exactly zero (used to loop forever)
        sp = harness((sp, ps) -> 1.0 + sum(abs2, ps), [0.0, 0.0], 1.0)
        ScalingCollapse.error_analysis!(sp)
        @test all(isapprox.(sp.optimal_ps_error, sqrt(0.3); rtol=1e-6))

        # negative optimum
        sp = harness((sp, ps) -> 1.0 + sum(abs2, ps .+ 2.0), [-2.0, -2.0], 1.0)
        ScalingCollapse.error_analysis!(sp)
        @test all(isapprox.(sp.optimal_ps_error, sqrt(0.3); rtol=1e-6))

        # quality does not depend on the second parameter
        sp = harness((sp, ps) -> 1.0 + ps[1]^2, [0.0, 1.0], 1.0)
        @test_logs (:warn, r"parameter b") ScalingCollapse.error_analysis!(sp)
        @test isapprox(sp.optimal_ps_error[1], sqrt(0.3); rtol=1e-6)
        @test isnan(sp.optimal_ps_error[2])

        # quality is NaN away from the optimum
        sp = harness((sp, ps) -> iszero(ps) ? 1.0 : NaN, [0.0, 0.0], 1.0)
        @test_logs (:warn,) (:warn,) ScalingCollapse.error_analysis!(sp)
        @test all(isnan, sp.optimal_ps_error)

        # minimum of zero
        sp = harness((sp, ps) -> sum(abs2, ps), [0.0, 0.0], 0.0)
        @test_logs (:warn, r"Skipping error analysis") ScalingCollapse.error_analysis!(sp)
        @test all(isnan, sp.optimal_ps_error)
    end
end
