using Random
# (Logging is not a test dependency; its loggers live in Base)
using Base.CoreLogging: NullLogger, with_logger

# a quality function that always fails
struct _ThrowingQuality <: ScalingCollapse.QualityFunction end
(::_ThrowingQuality)(sp, ps; check_bounds=false) = throw(DomainError(ps, "quality fails"))

# The threaded scan, fits and error analyses must give exactly the results of threads=false.
# With a single thread both run serially, so run the tests with several threads (e.g.
# JULIA_NUM_THREADS=2) to compare them.
@testset "threads" begin
    rng = Xoshiro(17)
    Ls = [8, 16, 32, 64]
    xs = collect(range(0.4, 0.6, length=41))
    ys = [tanh.((xs .- 0.5) .* L^(1 / 1.1)) .* L^-0.1 .+ 0.01 .* randn(rng, length(xs))
          for L in Ls]
    es = [fill(0.01, length(xs)) for _ in Ls]
    kw = (; sf=ScalingFunction(:xy), p_space=[0.45:0.01:0.55, 0.8:0.1:1.5, 0.0:0.05:0.3])

    @testset "keyword" begin
        @test_throws "threads must be true or false" ScalingProblem(xs, ys, es, Ls; kw...,
            threads=1)
        @test ScalingProblem(xs, ys, es, Ls; kw..., error=false).threads
        @test !ScalingProblem(xs, ys, es, Ls; kw..., error=false, threads=false).threads
    end

    # (rng=Xoshiro(1) is created for every fit, as the bootstrap advances it)
    fit(threads; kwargs...) = with_logger(NullLogger()) do
        ScalingProblem(xs, ys, es, Ls; kw..., rng=Xoshiro(1), threads, kwargs...)
    end
    @testset "same results as threads=false: $name" for (name, kwargs) in (
        "profile" => (;),
        "conditional" => (; error=:conditional),
        "bootstrap" => (; error=:bootstrap, n_bootstrap=20),
        "several starting points" => (; starting_ps=[[0.5, 1.0, 0.1], [0.48, 1.2, 0.2]]),
        "Linear" => (; quality=Linear(), quality_scan=Linear(), error=:conditional),
    )
        serial, threaded = fit(false; kwargs...), fit(true; kwargs...)
        @test isequal(threaded.starting_ps, serial.starting_ps)
        @test isequal(threaded.optimal_ps, serial.optimal_ps)
        @test isequal(threaded.minimum, serial.minimum)
        @test isequal(threaded.optimal_ps_error, serial.optimal_ps_error)
        @test all(isfinite, serial.optimal_ps_error)
        p_space = [range(p - 0.01, p + 0.01, length=5) for p in serial.optimal_ps]
        @test isequal(residuals(threaded; p_space), residuals(serial; p_space))
        @test isequal(parameter_errors(threaded; method=:bootstrap, n_bootstrap=10,
                rng=Xoshiro(2)),
            parameter_errors(serial; method=:bootstrap, n_bootstrap=10, rng=Xoshiro(2)))
    end

    @testset "exceptions in tasks are rethrown as they are" begin
        @test_throws DomainError ScalingProblem(xs, ys, es, Ls; kw...,
            quality_scan=_ThrowingQuality())
        @test_throws DomainError ScalingProblem(xs, ys, es, Ls; kw...,
            quality=_ThrowingQuality(), starting_ps=[[0.5, 1.0, 0.1], [0.48, 1.2, 0.2]])
    end
end
