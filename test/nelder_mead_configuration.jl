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
        # the AffineSimplexer vertex 0.6125 leaves p_space and is reflected to 0.3875, which is
        # inside p_space but has Inf quality -> must be halved toward 0.5
        f(x) = x[1] < 0.45 ? Inf : sum(abs2, x .- [0.5, 1.4])
        S = ScalingCollapse.BoundedSimplexer(
            ScalingCollapse.Optim.AffineSimplexer(a=0.025, b=0.1), [0.3, 1.2], [0.6, 1.8], f)
        simplex = ScalingCollapse.Optim.simplexer(S, [0.5, 1.4])
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
end
