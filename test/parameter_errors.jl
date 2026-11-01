using Random

# narrow tilted valley with its minimum 1 at a = b = 0.6 (independent of the data), as in
# error_methods.jl: reported at [0, 0], the profile of a re-fits b to a lower quality
struct _Valley <: ScalingCollapse.QualityFunction end
(::_Valley)(sp, ps; check_bounds=false) = 1 + 100 * (ps[1] - ps[2])^2 + 0.01 * (ps[2] - 0.6)^2

# the data (_bootstrap_synthetic, _BOOTSTRAP_KW) come from bootstrap_errors.jl
@testset "parameter_errors" begin
    xs, ys, es, Ls_p = _bootstrap_synthetic(Xoshiro(1))
    fit(; kw...) = ScalingProblem(xs, ys, es, Ls_p; _BOOTSTRAP_KW..., kw...)

    @testset "equal to the constructor" begin
        sp = fit(error=false)
        for method in (:profile, :conditional)
            @test parameter_errors(sp; method) == fit(error=method).optimal_ps_error
        end
        @test parameter_errors(sp; method=:bootstrap, n_bootstrap=20, rng=Xoshiro(5)) ==
              fit(error=:bootstrap, n_bootstrap=20, rng=Xoshiro(5)).optimal_ps_error
        @test parameter_errors(sp; method=:profile, error_threshold=1.5) ==
              fit(error_threshold=1.5).optimal_ps_error

        # the defaults are the settings of sp, and :profile after error=false
        @test parameter_errors(sp) == fit().optimal_ps_error
        sp_c = fit(error=:conditional, error_threshold=[1.5, 2.0])
        @test parameter_errors(sp_c) == sp_c.optimal_ps_error
        sp_b = fit(error=:bootstrap, n_bootstrap=20, rng=Xoshiro(5))
        @test parameter_errors(sp_b; rng=Xoshiro(5)) == sp_b.optimal_ps_error
    end

    @testset "sp unchanged" begin
        sp = fit(error=:conditional, rng=Xoshiro(2))
        before = (copy(sp.optimal_ps), copy(sp.optimal_ps_error), sp.minimum, sp.error,
            copy(sp.error_threshold), sp.n_bootstrap, copy(sp.rng), copy(sp.starting_ps))
        data = sp.data
        ys_before = [copy(d.ys) for d in data]
        parameter_errors(sp; method=:bootstrap, n_bootstrap=20, rng=Xoshiro(3))
        parameter_errors(sp; method=:profile, error_threshold=2.0)
        @test (sp.optimal_ps, sp.optimal_ps_error, sp.minimum, sp.error, sp.error_threshold,
            sp.n_bootstrap, sp.rng, sp.starting_ps) == before
        @test sp.data === data
        @test [d.ys for d in sp.data] == ys_before
    end

    @testset "parameter_errors! stores the errors and settings" begin
        sp = fit(error=false)
        σ = parameter_errors!(sp; method=:conditional, error_threshold=1.5)
        @test σ == sp.optimal_ps_error == fit(error=:conditional, error_threshold=1.5).optimal_ps_error
        @test sp.error === :conditional
        @test sp.error_threshold == [1.5, 1.5]
        σ[1] = -1.0
        @test sp.optimal_ps_error[1] != -1.0

        rng = Xoshiro(4)
        parameter_errors!(sp; method=:bootstrap, n_bootstrap=20, rng)
        @test sp.error === :bootstrap
        @test sp.n_bootstrap == 20
        @test sp.rng === rng
    end

    @testset "keyword validation" begin
        sp = fit(error=false)
        @test_throws "method must be" parameter_errors(sp; method=false)
        @test_throws ArgumentError parameter_errors(sp; method=:bogus)
        @test_throws "error_threshold must be larger than 1" parameter_errors(sp; error_threshold=0.9)
        @test_throws "error_threshold has 3 entries" parameter_errors(sp; error_threshold=[1.3, 1.3, 1.3])
        @test_throws "n_bootstrap must be" parameter_errors(sp; method=:bootstrap, n_bootstrap=1)
        @test_throws "rng must be" parameter_errors(sp; method=:bootstrap, rng=42)
        @test_throws MethodError parameter_errors(sp; error=:profile)
        # invalid keywords leave sp unchanged
        @test_throws ArgumentError parameter_errors!(sp; method=:conditional, error_threshold=0.9)
        @test sp.error === :none

        sp_no_es = ScalingProblem(xs, ys, Ls_p; _BOOTSTRAP_KW..., error=false)
        @test_throws "error bars" parameter_errors(sp_no_es; method=:bootstrap)
    end

    @testset "lower quality during the profile" begin
        # a fit that is not at the minimum of the valley
        function reported_valley()
            sp = ScalingProblem(xs, ys, es, Ls_p; _BOOTSTRAP_KW..., quality=_Valley(), error=false)
            sp.optimal_ps = [0.0, 0.0]
            sp.minimum = _Valley()(sp, sp.optimal_ps)
            return sp
        end

        sp = reported_valley()
        @test_throws "parameter_errors!" parameter_errors(sp)
        @test sp.optimal_ps == [0.0, 0.0]
        @test sp.minimum == _Valley()(sp, [0.0, 0.0])
        @test sp.optimal_ps_error == [0.0, 0.0]

        # the conditional errors and the bootstrap do not refit
        @test all(isfinite, parameter_errors(sp; method=:conditional))
        @test sp.optimal_ps == [0.0, 0.0]

        σ = @test_logs (:warn, r"lower quality") parameter_errors!(sp)
        @test sp.minimum ≈ 1.0 atol = 1e-4
        @test sp.optimal_ps ≈ [0.6, 0.6] atol = 1e-2
        @test σ == sp.optimal_ps_error
        @test σ ≈ fill(sqrt(0.3 * 100.01), 2) rtol = 1e-2
        # at the re-optimized fit, the profile finds no lower quality
        @test parameter_errors(sp) == σ
    end
end
