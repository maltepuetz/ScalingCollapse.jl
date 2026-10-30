# Minimal stand-in for a ScalingProblem: the fields error_analysis! and optimize_starting_ps! read.
Base.@kwdef mutable struct ErrorHarness
    quality
    optimal_ps = [0.0, 0.0]
    optimal_ps_error = zeros(2)
    minimum = 1.0
    error_threshold = [1.3, 1.3]
    p_space = [-10.0:0.1:10.0, -10.0:0.1:10.0]
    starting_ps = [0.0, 0.0]
    verbose = false
    error = :profile
    sf = ScalingFunction(:x; p_names=["a", "b"])
    optimizer = ScalingCollapse.Optim.NelderMead()
    optimizer_options = ScalingCollapse.Optim.Options()
end

@testset "error methods" begin
    # Fixtures recorded on main before the profile-error work: the default error of that time is
    # now `error=:conditional`. Compared with rtol (not ==) because the literals come from one
    # machine and CI also runs aarch64 and 32-bit.
    @testset "conditional errors unchanged" begin
        # README Ising fit (Ts, binder, Ls come from ising.jl)
        sp = ScalingProblem(Ts, binder, Ls;
            sf=ScalingFunction(:x; p_names=["T_c", "nu"]),
            p_space=[2.0:0.1:2.4, 0.5:0.1:1.5],
            dx=[-1.0, 1.0],
            error=:conditional,
        )
        @test isapprox(sp.optimal_ps_error, [0.0016435929762187492, 0.019317580143210902]; rtol=1e-8)

        # second CDW fit (data from cdw.jl)
        sp = ScalingProblem(
            [X_L4, X_L6, X_L8, X_L10, X_L12],
            [Y_L4, Y_L6, Y_L8, Y_L10, Y_L12],
            [E_L4, E_L6, E_L8, E_L10, E_L12],
            [4, 6, 8, 10, 12];
            sf=ScalingFunction(:xny, p_names=["x_c", "nu", "gamma"], nu=1),
            p_space=[5:0.1:7, 1:0.1:3],
            dx=[-2, 2],
            error=:conditional,
        )
        @test isapprox(sp.optimal_ps_error, [0.03172257936100857, 0.03636215786455168]; rtol=1e-8)

        # bond percolation, largest-cluster fraction (collapses/cases_percolation.jl: bond_strength)
        h5open(joinpath(@__DIR__, "testdata", "percolation", "bond_square.h5"), "r") do h
            Lb = read(h, "Ls")
            xs, ys, es = ([read(h, "strength/L$L/$v") for L in Lb] for v in ("x", "y", "e"))
            sp = ScalingProblem(xs, ys, es, Lb;
                sf=ScalingFunction(:xy),
                p_space=[0.48:0.002:0.52, 0.8:0.1:2.0, 0.05:0.025:0.3],
                dx=[-1.0, 1.0],
                error=:conditional,
            )
            @test isapprox(sp.optimal_ps_error,
                [0.00010935007760675173, 0.013050655089099956, 0.003100701294555539]; rtol=1e-8)
        end

        # synthetic clean_xy data, Linear quality, no errors (collapses/cases_synthetic.jl)
        h5open(joinpath(@__DIR__, "testdata", "synthetic", "synthetic.h5"), "r") do h
            Ls_syn = read(h, "clean_xy/Ls")
            g = h["clean_xy"]
            xs = [read(g, "L$L/x") for L in Ls_syn]
            ys = [read(g, "L$L/y") for L in Ls_syn]
            sp = ScalingProblem(xs, ys, Ls_syn;
                sf=ScalingFunction(:xy),
                p_space=[0.5:0.01:0.9, 0.6:0.1:2.0, 0.0:0.05:0.6],
                dx=[-1.5, 1.5],
                quality=Linear(),
                error=:conditional,
            )
            @test isapprox(sp.optimal_ps_error,
                [3.0086057594151985e-5, 0.002085947931419029, 0.000520178097105716]; rtol=1e-8)
        end
    end

    @testset "error keyword" begin
        kw = (sf=ScalingFunction(:x), dx=[-1.0, 1.0], p_space=[2.0:0.1:2.4, 0.5:0.1:1.5])
        @test_throws ArgumentError ScalingProblem(Ts, binder, Ls; kw..., error=:bogus)

        sp = ScalingProblem(Ts, binder, Ls; kw..., error=false)
        @test sp.error === :none
        @test all(iszero, sp.optimal_ps_error)

        sp = ScalingProblem(Ts, binder, Ls; kw..., error=true)
        @test sp.error === :profile

        sp = ScalingProblem(Ts, binder, Ls; kw..., error=:conditional)
        @test sp.error === :conditional
    end

    @testset "profile width" begin
        errors!(sp) = (ScalingCollapse.error_analysis!(sp); sp.optimal_ps_error)
        correlated(ρ) = (sp, ps; check_bounds=false) -> 1 + ps[1]^2 + ps[2]^2 + 2ρ * ps[1] * ps[2]

        @testset "correlated quadratic" begin
            @test errors!(ErrorHarness(quality=correlated(0.8))) ≈
                  fill(sqrt(0.3 / (1 - 0.64)), 2) rtol = 1e-3
            @test errors!(ErrorHarness(quality=correlated(0.8), error=:conditional)) ≈
                  fill(sqrt(0.3), 2)
        end

        @testset "one free parameter" begin
            kw = (sf=ScalingFunction(:x; p_names=["T_c", "nu"], nu=1),
                p_space=[2.0:0.1:2.4], dx=[-1.0, 1.0])
            sp_profile = ScalingProblem(Ts, binder, Ls; kw..., error=:profile)
            sp_conditional = ScalingProblem(Ts, binder, Ls; kw..., error=:conditional)
            @test sp_profile.optimal_ps_error == sp_conditional.optimal_ps_error
        end

        @testset "two valleys" begin
            # a second valley around a = 1 lies beyond a barrier above the threshold
            sp = ErrorHarness(quality=(sp, ps; check_bounds=false) ->
                1 + ps[2]^2 + 25 * ps[1]^2 * (ps[1] - 1)^2)
            errors!(sp)
            @test sp.optimal_ps_error[1] ≈ 0.125224 rtol = 1e-3
        end

        @testset "invalid region" begin
            cutoff(c) = (sp, ps; check_bounds=false) -> abs(ps[1]) < c ? 1 + sum(abs2, ps) : Inf
            sp = ErrorHarness(quality=cutoff(0.5))
            @test_logs (:warn, r"undefined") ScalingCollapse.error_analysis!(sp)
            @test isnan(sp.optimal_ps_error[1])
            @test sp.optimal_ps_error[2] ≈ sqrt(0.3) rtol = 1e-3

            @test errors!(ErrorHarness(quality=cutoff(0.6)))[1] ≈ sqrt(0.3) rtol = 1e-3
        end

        @testset "nuisance bounds" begin
            # b is confined to [-0.1, 0.1]: for a > 0.125 the profile sits at b = ∓0.1
            q = correlated(0.8)
            sp = ErrorHarness(
                quality=(sp, ps; check_bounds=false) -> check_bounds ? Inf : q(sp, ps),
                p_space=[-10.0:0.1:10.0, -0.1:0.01:0.1],
            )
            # root of a² − 0.16a − 0.29 = 0
            @test errors!(sp)[1] ≈ (0.16 + sqrt(0.16^2 + 4 * 0.29)) / 2 rtol = 1e-3
        end

        @testset "per-parameter thresholds" begin
            sp = ErrorHarness(quality=(sp, ps; check_bounds=false) -> 1 + sum(abs2, ps),
                error_threshold=[1.3, 2.0])
            @test errors!(sp) ≈ [sqrt(0.3), 1.0] rtol = 1e-3
        end

        @testset "inner fit failure" begin
            # inner fits never converge in one iteration: one warning per parameter
            sp = ErrorHarness(quality=correlated(0.8),
                optimizer_options=ScalingCollapse.Optim.Options(iterations=1))
            @test_logs (:warn, r"parameter a may be underestimated") (:warn,
                r"parameter b may be underestimated") ScalingCollapse.error_analysis!(sp)
            @test all(sp.optimal_ps_error .>= sqrt(0.3) * (1 - 1e-3))
        end

        @testset "tiny conditional root" begin
            # a very steep valley (as for a minimum on a jump of the quality): the conditional
            # roots are ≈ 5e-13, which 60 outward steps of ×1.5 alone cannot carry to the
            # profile crossing; the profile 1 + a²K/(K+1) (same for b) crosses 1.3 at ≈ √0.3
            sp = ErrorHarness(quality=(sp, ps; check_bounds=false) ->
                1 + 1e24 * (ps[1] - ps[2])^2 + ps[2]^2)
            @test errors!(sp) ≈ fill(sqrt(0.3), 2) rtol = 1e-3
        end

        @testset "optimum at zero" begin
            sp = ErrorHarness(quality=(sp, ps; check_bounds=false) -> 1 + sum(abs2, ps))
            @test errors!(sp) ≈ fill(sqrt(0.3), 2) rtol = 1e-3
        end

        # profile errors are at least the conditional errors (up to the bisection precision)
        at_least_conditional(profile, conditional) =
            all(isfinite, profile) && all(profile .>= conditional .* (1 - 1e-3))

        @testset "profile ≥ conditional" begin
            kw = (sf=ScalingFunction(:x; p_names=["T_c", "nu"]),
                p_space=[2.0:0.1:2.4, 0.5:0.1:1.5], dx=[-1.0, 1.0])
            sp_profile = ScalingProblem(Ts, binder, Ls; kw..., error=:profile)
            sp_conditional = ScalingProblem(Ts, binder, Ls; kw..., error=:conditional)
            @test at_least_conditional(sp_profile.optimal_ps_error,
                sp_conditional.optimal_ps_error)
        end

        @testset "profile with starting_ps" begin
            kw = (sf=ScalingFunction(:x; p_names=["T_c", "nu"]),
                starting_ps=[2.27, 1.0], dx=[-1.0, 1.0])
            sp_profile = ScalingProblem(Ts, binder, Ls; kw..., error=:profile)
            sp_conditional = ScalingProblem(Ts, binder, Ls; kw..., error=:conditional)
            @test at_least_conditional(sp_profile.optimal_ps_error,
                sp_conditional.optimal_ps_error)
        end

        @testset "fixed parameter" begin
            # second CDW fit: nu fixed, x_c and gamma free
            fit(error) = ScalingProblem(
                [X_L4, X_L6, X_L8, X_L10, X_L12],
                [Y_L4, Y_L6, Y_L8, Y_L10, Y_L12],
                [E_L4, E_L6, E_L8, E_L10, E_L12],
                [4, 6, 8, 10, 12];
                sf=ScalingFunction(:xny, p_names=["x_c", "nu", "gamma"], nu=1),
                p_space=[5:0.1:7, 1:0.1:3],
                dx=[-2, 2],
                error,
            )
            @test at_least_conditional(fit(:profile).optimal_ps_error,
                fit(:conditional).optimal_ps_error)
        end
    end

    @testset "lower quality during the profile" begin
        # narrow tilted valley with its minimum 1 at a = b = 0.6, reported at [0, 0]: the
        # profile of a at its conditional root (≈ 0.055) re-fits b to a quality below 1.0036
        valley(sp, ps; check_bounds=false) =
            1 + 100 * (ps[1] - ps[2])^2 + 0.01 * (ps[2] - 0.6)^2
        reported_valley() = ErrorHarness(quality=valley, minimum=valley(nothing, [0.0, 0.0]))

        @testset "restart" begin
            sp = reported_valley()
            @test_logs (:warn, r"lower quality") ScalingCollapse.error_analysis!(sp)
            @test sp.minimum ≈ 1.0 atol = 1e-4
            @test sp.optimal_ps ≈ [0.6, 0.6] atol = 1e-2
            @test sp.optimal_ps_error ≈ fill(sqrt(0.3 * 100.01), 2) rtol = 1e-2
        end

        @testset "outside p_space" begin
            # the quality drops by 0.01 for |b| ≥ 0.5, which a only reaches outside p_space
            sp = ErrorHarness(
                quality=(sp, ps; check_bounds=false) ->
                    1 + (ps[1] - ps[2])^2 + (abs(ps[2]) >= 0.5 ? -0.01 : 0.0),
                p_space=[-0.12:0.01:0.12, -1.0:0.01:1.0],
            )
            @test_logs (:warn, r"outside p_space") ScalingCollapse.error_analysis!(sp)
            @test isnan(sp.optimal_ps_error[1])
            @test isfinite(sp.optimal_ps_error[2])
            @test sp.optimal_ps == [0.0, 0.0]
            @test sp.minimum == 1.0
        end

        @testset "restarts exhausted" begin
            sp = reported_valley()
            @test_logs (:warn, r"could not be confirmed") ScalingCollapse.error_analysis!(sp;
                max_restarts=0)
            @test all(isnan, sp.optimal_ps_error)
            @test sp.minimum < 1.0036
            @test sp.minimum == valley(sp, sp.optimal_ps)
            @test all(j -> first(sp.p_space[j]) <= sp.optimal_ps[j] <= last(sp.p_space[j]),
                eachindex(sp.optimal_ps))
        end
    end
end
