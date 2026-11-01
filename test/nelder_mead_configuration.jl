struct RecordingSimplexer <: ScalingCollapse.Optim.Simplexer
    calls::Base.RefValue{Int}
end

function ScalingCollapse.Optim.simplexer(simplexer::RecordingSimplexer, initial_x)
    simplexer.calls[] += 1
    return ScalingCollapse.Optim.simplexer(
        ScalingCollapse.Optim.AffineSimplexer(; a=0.01, b=0.1),
        initial_x,
    )
end

mutable struct OptimizationConfigurationHarness
    verbose::Bool
    p_space::Vector{Vector{Float64}}
    quality_scan::Function
    quality::Function
    starting_ps::Vector{Float64}
    minimum::Float64
    optimal_ps::Vector{Float64}
    optimizer
    optimizer_options
end

@testset "Nelder-Mead configuration" begin
    @testset "defaults and validation" begin
        sp = ScalingProblem(
            Ts,
            binder,
            Ls;
            sf=ScalingFunction(:x),
            starting_ps=[2.27, 1.0],
            error=false,
        )

        @test sp.optimizer isa ScalingCollapse.Optim.NelderMead
        @test sp.optimizer.initial_simplex isa ScalingCollapse.Optim.AffineSimplexer
        @test sp.optimizer.initial_simplex.a == 0.025
        @test sp.optimizer.initial_simplex.b == 0.1
        @test sp.optimizer_options isa ScalingCollapse.Optim.Options
        @test sp.optimizer_options.iterations == 1_000

        @test_throws ArgumentError ScalingProblem(
            ; optimizer=ScalingCollapse.Optim.BFGS(),
        )
        @test_throws ArgumentError ScalingProblem(; optimizer_options=(; iterations=10))
    end

    @testset "custom method and options are used for the final solve" begin
        simplex_calls = Ref(0)
        callback_calls = Ref(0)
        optimizer = ScalingCollapse.Optim.NelderMead(
            initial_simplex=RecordingSimplexer(simplex_calls),
            parameters=ScalingCollapse.Optim.FixedParameters(
                α=1.0,
                β=2.0,
                γ=0.5,
                δ=0.5,
            ),
        )
        optimizer_options = ScalingCollapse.Optim.Options(
            iterations=2,
            callback=trace -> begin
                callback_calls[] += 1
                return false
            end,
            show_warnings=false,
        )

        sp = ScalingProblem(
            Ts,
            binder,
            Ls;
            sf=ScalingFunction(:x),
            starting_ps=[2.27, 1.0],
            error=false,
            optimizer=optimizer,
            optimizer_options=optimizer_options,
        )

        @test sp.optimizer === optimizer
        @test sp.optimizer_options === optimizer_options
        @test simplex_calls[] > 0
        @test callback_calls[] > 0
        @test callback_calls[] <= optimizer_options.iterations + 1
    end

    @testset "custom method and options are used for scan refinement" begin
        simplex_calls = Ref(0)
        callback_calls = Ref(0)
        optimizer = ScalingCollapse.Optim.NelderMead(
            initial_simplex=RecordingSimplexer(simplex_calls),
        )
        optimizer_options = ScalingCollapse.Optim.Options(
            iterations=2,
            callback=trace -> begin
                callback_calls[] += 1
                return false
            end,
            show_warnings=false,
        )
        objective = (sp, ps; check_bounds) -> sum(abs2, ps)
        sp = OptimizationConfigurationHarness(
            false,
            [[-1.0, 0.0, 1.0], [-1.0, 0.0, 1.0]],
            objective,
            objective,
            [1.0, 1.0],
            Inf,
            zeros(2),
            optimizer,
            optimizer_options,
        )

        ScalingCollapse.parameter_scan!(sp)

        @test simplex_calls[] > 0
        @test callback_calls[] > 0
        @test callback_calls[] <= optimizer_options.iterations + 1
        @test all(isfinite, sp.starting_ps)
    end

    @testset "initial simplex avoids non-finite quality" begin
        # AffineSimplexer(a=0, b=0.25) puts the first vertex at 1.25 * 0.5 = 0.625, outside
        # p_space; it is reflected to 0.375, which is inside p_space but has Inf quality, and
        # must then be halved toward 0.5 (0.4375 is still Inf, 0.46875 is finite)
        f(x) = x[1] < 0.45 ? Inf : sum(abs2, x .- [0.5, 1.4])
        base = ScalingCollapse.Optim.AffineSimplexer(a=0.0, b=0.25)
        @test !all([0.3, 1.2] .<= ScalingCollapse.Optim.simplexer(base, [0.5, 1.4])[2] .<= [0.6, 1.8])
        S = ScalingCollapse.BoundedSimplexer(base, [0.3, 1.2], [0.6, 1.8], f)
        simplex = ScalingCollapse.Optim.simplexer(S, [0.5, 1.4])
        @test simplex[2] ≈ [0.46875, 1.4]
        @test all(v -> isfinite(f(v)), simplex)
        @test all(v -> all([0.3, 1.2] .<= v .<= [0.6, 1.8]), simplex)
        @test length(unique(simplex)) == 3

        # never finite away from the starting point: terminates, keeps vertices distinct
        g(x) = x == [0.5, 1.4] ? 0.0 : Inf
        S = ScalingCollapse.BoundedSimplexer(
            ScalingCollapse.Optim.AffineSimplexer(a=0.025, b=0.1), [0.3, 1.2], [0.6, 1.8], g)
        simplex = ScalingCollapse.Optim.simplexer(S, [0.5, 1.4])
        @test length(unique(simplex)) == 3
    end

    @testset "StepSimplexer" begin
        S = StepSimplexer([0.1, 0.01])
        @test ScalingCollapse.Optim.simplexer(S, [1.0, 2.0]) ≈
              [[1.0, 2.0], [1.1, 2.0], [1.0, 2.01]]
        @test_throws DimensionMismatch ScalingCollapse.Optim.simplexer(S, [1.0, 2.0, 3.0])
        @test_throws ArgumentError StepSimplexer([0.1, 0.0])
        @test_throws ArgumentError StepSimplexer([0.1, Inf])
        @test StepSimplexer([1, -2]).steps == [1.0, -2.0]

        # fits of a subset of the free parameters (profile errors) use their steps only
        @test ScalingCollapse._restrict_simplexer(StepSimplexer([0.1, 0.2, 0.3]), [1, 3]).steps ==
              [0.1, 0.3]
        affine = ScalingCollapse.Optim.AffineSimplexer()
        @test ScalingCollapse._restrict_simplexer(affine, [1]) === affine

        # one step per free parameter (fixed parameters have none)
        step_optimizer(steps) =
            ScalingCollapse.Optim.NelderMead(initial_simplex=StepSimplexer(steps))
        @test_throws ArgumentError ScalingProblem(Ts, binder, Ls; sf=ScalingFunction(:x),
            starting_ps=[2.27, 1.0], optimizer=step_optimizer([0.01, 0.05, 0.1]))
        @test_throws ArgumentError ScalingProblem(Ts, binder, Ls;
            sf=ScalingFunction(:x; nu=1.0), starting_ps=[2.27],
            optimizer=step_optimizer([0.01, 0.05]))

        # the profile fits of T_c (a valley of width ~0.002) use its small step: a brute-force
        # profile over a fine T_c grid reaches 1.3 times the minimum at Δnu ≈ -0.033 and
        # +0.027. The default simplex (step 0.25 in T_c) misses this valley and returns the
        # conditional width 0.019 for nu.
        kwargs = (sf=ScalingFunction(:x), starting_ps=[2.27, 1.0], dx=[-1.0, 1.0])
        sp_conditional = ScalingProblem(Ts, binder, Ls; kwargs..., error=:conditional)
        sp = ScalingProblem(Ts, binder, Ls; kwargs..., optimizer=step_optimizer([0.01, 0.05]))
        @test sp.error === :profile
        @test all(abs.(sp.optimal_ps .- sp_conditional.optimal_ps) .<
                  0.25 .* sp.optimal_ps_error)
        @test all(sp.optimal_ps_error .>= sp_conditional.optimal_ps_error)
        @test sp.optimal_ps_error[2] ≈ 0.033 rtol = 0.1
    end

    @testset "parameter scan builds its simplex with quality_scan" begin
        # the scan must shrink simplex vertices against quality_scan, never quality
        sp = OptimizationConfigurationHarness(
            false,
            [[-1.0, -0.5, 0.0, 0.5, 1.0], [-1.0, -0.5, 0.0, 0.5, 1.0]],
            (sp, ps; check_bounds) -> sum(abs2, ps),
            (sp, ps; check_bounds) -> error("quality must not be evaluated during the scan"),
            zeros(2),
            Inf,
            zeros(2),
            ScalingCollapse.Optim.NelderMead(),
            ScalingCollapse.Optim.Options(),
        )
        ScalingCollapse.parameter_scan!(sp)
        @test all(abs.(sp.starting_ps) .< 1e-3)
    end
end
