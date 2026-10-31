using Random
using Statistics

# synthetic Binder-like data: y = tanh((x - x_c)/x_c * L^(1/nu)) + noise of width sigma
function _bootstrap_synthetic(rng; x_c=0.5, nu=1.0, sigma=0.02)
    Ls_b = [8, 16, 32]
    xs = collect(0.4:0.005:0.6)
    ys = [tanh.((xs .- x_c) ./ x_c .* L^(1 / nu)) .+ sigma .* randn(rng, length(xs)) for L in Ls_b]
    es = [fill(sigma, length(xs)) for _ in Ls_b]
    return xs, ys, es, Ls_b
end

# Linear quality that is Inf when `fails` returns true for the first y value of the data, to
# make chosen bootstrap refits fail
struct _FailingLinear <: ScalingCollapse.QualityFunction
    fails::Function
end
(qf::_FailingLinear)(sp, ps; check_bounds=false) =
    qf.fails(sp.data[1].ys[1]) ? Inf : Linear()(sp, ps; check_bounds)

const _BOOTSTRAP_KW = (
    sf=ScalingFunction(:x; p_names=["x_c", "nu"]),
    starting_ps=[0.5, 1.0],
    dx=[-1.0, 1.0],
)

@testset "bootstrap errors" begin
    @testset "requires error bars" begin
        xs, ys, _, Ls_b = _bootstrap_synthetic(Xoshiro(1))
        @test_throws "error bars" ScalingProblem(xs, ys, Ls_b; _BOOTSTRAP_KW..., error=:bootstrap)
    end

    @testset "keyword validation" begin
        xs, ys, es, Ls_b = _bootstrap_synthetic(Xoshiro(1))
        @test_throws "n_bootstrap must be" ScalingProblem(xs, ys, es, Ls_b; _BOOTSTRAP_KW...,
            error=:bootstrap, n_bootstrap=1)
        @test_throws "rng must be" ScalingProblem(xs, ys, es, Ls_b; _BOOTSTRAP_KW...,
            error=:bootstrap, rng=42)
    end

    @testset "reproducible with a seeded rng" begin
        xs, ys, es, Ls_b = _bootstrap_synthetic(Xoshiro(1))
        fit(seed) = ScalingProblem(xs, ys, es, Ls_b; _BOOTSTRAP_KW...,
            error=:bootstrap, n_bootstrap=20, rng=Xoshiro(seed))
        sp1, sp2 = fit(7), fit(7)
        @test sp1.error === :bootstrap
        @test sp1.optimal_ps_error == sp2.optimal_ps_error
        @test all(isfinite, sp1.optimal_ps_error) && all(>(0), sp1.optimal_ps_error)
    end

    @testset "fit of the original data unchanged" begin
        xs, ys, es, Ls_b = _bootstrap_synthetic(Xoshiro(1))
        sp_b = ScalingProblem(xs, ys, es, Ls_b; _BOOTSTRAP_KW...,
            error=:bootstrap, n_bootstrap=20, rng=Xoshiro(3))
        sp_n = ScalingProblem(xs, ys, es, Ls_b; _BOOTSTRAP_KW..., error=false)
        @test sp_b.optimal_ps == sp_n.optimal_ps
        @test sp_b.minimum == sp_n.minimum
        @test all(d -> d.ys == ys[findfirst(==(d.L), Ls_b)], sp_b.data)
    end

    @testset "matches the spread of independent noise realisations" begin
        rng = Xoshiro(2024)
        fits = map(1:100) do _
            xs, ys, es, Ls_b = _bootstrap_synthetic(rng)
            ScalingProblem(xs, ys, es, Ls_b; _BOOTSTRAP_KW..., error=false).optimal_ps
        end
        true_spread = [std(getindex.(fits, i)) for i in 1:2]

        # the bootstrap error of a single dataset scatters, so average over several datasets
        bootstrap = map(1:10) do k
            xs, ys, es, Ls_b = _bootstrap_synthetic(rng)
            ScalingProblem(xs, ys, es, Ls_b; _BOOTSTRAP_KW...,
                error=:bootstrap, n_bootstrap=50, rng=Xoshiro(k)).optimal_ps_error
        end
        ratio = mean(bootstrap) ./ true_spread
        @test all(r -> 0.7 < r < 1.5, ratio)
    end

    @testset "failed refits" begin
        xs, ys, es, Ls_b = _bootstrap_synthetic(Xoshiro(1))
        y0 = ys[1][1]

        # a few refits fail (resampled y more than 2 error bars above y0): dropped, warning
        sp = @test_logs (:warn, r"were dropped") match_mode=:any ScalingProblem(
            xs, ys, es, Ls_b; _BOOTSTRAP_KW..., quality=_FailingLinear(y -> y > y0 + 0.04),
            error=:bootstrap, n_bootstrap=100, rng=Xoshiro(1))
        @test all(isfinite, sp.optimal_ps_error)

        # every refit fails: errors are NaN, warning
        sp = @test_logs (:warn, r"Bootstrap error analysis failed") match_mode=:any ScalingProblem(
            xs, ys, es, Ls_b; _BOOTSTRAP_KW..., quality=_FailingLinear(y -> y != y0),
            error=:bootstrap, n_bootstrap=20, rng=Xoshiro(1))
        @test all(isnan, sp.optimal_ps_error)
        @test isfinite(sp.minimum)
    end

    @testset "CDW data with error bars" begin
        sp = ScalingProblem(
            [X_L4, X_L6, X_L8, X_L10, X_L12],
            [Y_L4, Y_L6, Y_L8, Y_L10, Y_L12],
            [E_L4, E_L6, E_L8, E_L10, E_L12],
            [4, 6, 8, 10, 12];
            sf=ScalingFunction(:xny, p_names=["x_c", "nu", "gamma"], nu=1),
            p_space=[5:0.1:7, 1:0.1:3],
            dx=[-2, 2],
            error=:bootstrap,
            n_bootstrap=30,
            rng=Xoshiro(11),
        )
        @test all(isfinite, sp.optimal_ps_error) && all(>(0), sp.optimal_ps_error)
    end
end
