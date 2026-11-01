@testset "Several starting points" begin
    # p shifts the scaled x values of size L by (L - 1) * h(p). h has a root near p = -1
    # (global minimum of the quality) and a positive local minimum near p = 1, so a single
    # start at p = 1 gets stuck there. The data are perturbed so the minimum is positive and
    # the error analysis runs.
    h(p) = (p^2 - 1)^2 + 0.3 * (p + 1)
    sf = ScalingFunction(
        (d, p) -> ScalingCollapse.Data(d.L, d.xs .+ (d.L - 1) * h(p), d.ys, d.es);
        p_names=["p"],
    )
    Ls = [1, 2, 3]
    x = collect(-3:0.1:3)
    ys = [tanh.(x) .+ 0.01 .* sin.(7 * L .* x) for L in Ls]
    problem(starting_ps; kwargs...) = ScalingProblem([x for _ in Ls], ys, Ls;
        sf=sf, starting_ps=starting_ps, kwargs...,
    )

    @testset "two minima" begin
        local_fit = problem([1.0])
        global_fit = problem([-1.0])
        @test global_fit.minimum < local_fit.minimum / 100
        @test isapprox(global_fit.optimal_ps[1], -1.07; atol=0.02)
        @test all(isfinite, global_fit.optimal_ps_error)

        # the best start wins in either order, and the error analysis runs on its fit
        for starts in ([[1.0], [-1.0]], [[-1.0], [1.0]])
            sp = problem(starts)
            @test sp.starting_ps == [-1.0]
            @test sp.optimal_ps == global_fit.optimal_ps
            @test sp.minimum == global_fit.minimum
            @test sp.optimal_ps_error == global_fit.optimal_ps_error
        end

        # integer starts and untyped vectors are converted
        @test problem([[1], [-1]]; error=false).optimal_ps == global_fit.optimal_ps
        @test problem(Any[-1.0]; error=false).optimal_ps == global_fit.optimal_ps
        @test problem(Any[Any[1.0], [-1]]; error=false).optimal_ps == global_fit.optimal_ps
    end

    @testset "starts without a finite quality" begin
        # at p = 10 the rescaled sizes do not overlap
        sp = problem([[10.0], [-1.0]]; error=false)
        @test sp.starting_ps == [-1.0]
        @test isfinite(sp.minimum)

        err = try
            problem([[10.0], [20.0]]; error=false)
            nothing
        catch e
            e
        end
        @test err isa ErrorException
        @test occursin("[[10.0], [20.0]]", err.msg)
    end

    @testset "invalid starting points" begin
        @test_throws ArgumentError problem([[1.0], [1.0, 2.0]])
        @test_throws ArgumentError problem(Vector{Float64}[])
        @test_throws ArgumentError problem([])
        @test_throws ArgumentError problem("a")
        msg = try
            problem([[1.0], [1.0, 2.0]])
        catch e
            e.msg
        end
        @test occursin("starting_ps[2]", msg)
    end
end
