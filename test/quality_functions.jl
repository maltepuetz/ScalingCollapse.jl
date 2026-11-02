using Random


@testset "Quality functions" begin

    @testset "Linear with coinciding scaled x values" begin
        # grid uniform in z = (p - p_c)/p_c * L: at p1 = p_c every size has a point at
        # x = 0, and at nu = 1 all scaled x values coincide across the sizes
        Ls = [8, 16, 32]
        zs = collect(-2:0.05:2)
        xs = [0.5 .* (1 .+ zs ./ L) for L in Ls]
        ys = [tanh.(zs) for _ in Ls]
        es = [fill(0.01, length(zs)) for _ in Ls]

        for (label, args) in (("without errors", (xs, ys, Ls)), ("with errors", (xs, ys, es, Ls)))
            @testset "$label" begin
                sp = ScalingProblem(args...;
                    sf=ScalingFunction(:x), quality=Linear(), starting_ps=[0.5, 1.0], error=false,
                )
                for nu in (0.5, 1.0, 1.5)
                    @test isfinite(sp.quality(sp, [0.5, nu]))
                end
                # the data collapse exactly at the true parameters
                @test sp.quality(sp, [0.5, 1.0]) ≈ 0.0 atol = 1e-12
                @test sp.quality(sp, [0.5, 1.0]) < sp.quality(sp, [0.4875, 1.0])
                @test isapprox(sp.optimal_ps[1], 0.5; atol=1e-3)
                @test isapprox(sp.optimal_ps[2], 1.0; atol=1e-2)
            end
        end
    end

    @testset "Linear on non-tied data" begin
        # two sizes interleave in x with y = x^2; strictly inside the overlap (1, 6) lie
        # the sorted points x = 2, 3, 4, 5
        x1, x2 = [0.0, 2.0, 4.0, 6.0], [1.0, 3.0, 5.0, 7.0]
        sp = ScalingProblem([x1, x2], [x1 .^ 2, x2 .^ 2], [1, 2];
            sf=ScalingFunction((d, p1, p2) -> d; N_parameters=2), quality=Linear(),
            starting_ps=[1.0, 1.0], error=false, dx=[-10, 10],
        )
        # interior points: ybar(3) = 10, ybar(4) = 17 -> S = (1 + 1) / 2
        @test sp.quality(sp, [1.0, 1.0]) ≈ 1.0
    end

    @testset "SingleSpline" begin
        identity_sf = ScalingFunction((d, p1, p2) -> d; N_parameters=2)
        spline_problem(xs, ys, Ls) = ScalingProblem(xs, ys, Ls;
            sf=identity_sf, starting_ps=[1.0, 1.0], error=false, dx=[-10, 10],
        )

        # points outside the overlap interval do not change the fit
        x1, x2 = [0.0, 2.0, 4.0, 6.0], [1.0, 3.0, 5.0, 7.0]
        sp = spline_problem([x1, x2], [x1 .^ 2, x2 .^ 2], [1, 2])
        x2_long = [1.0, 3.0, 5.0, 7.0, 8.0, 9.0]
        sp_long = spline_problem([x1, x2_long], [x1 .^ 2, x2_long .^ 2], [1, 2])
        for q in (SingleSpline(), SingleSpline(N_knots=3))
            @test isfinite(q(sp, [1.0, 1.0]))
            @test q(sp, [1.0, 1.0]) == q(sp_long, [1.0, 1.0])
        end

        # more knots than data points are capped at one less than the 6 points in (1, 6)
        @test SingleSpline(N_knots=100)(sp, [1.0, 1.0]) == SingleSpline(N_knots=5)(sp, [1.0, 1.0])

        # p1 shifts the scaled x values of size L by (L - 1) * p1, so the problems can be
        # solved at p1 = 0 and the degenerate overlaps are reached by shifting
        shift_sf = ScalingFunction(
            (d, p1, p2) -> ScalingCollapse.Data(d.L, d.xs .+ (d.L - 1) * p1, d.ys, d.es);
            N_parameters=2,
        )
        x = [0.0, 1.0, 2.0, 3.0, 4.0]
        shift_problem(Ls) = ScalingProblem([x for _ in Ls], [x for _ in Ls], Ls;
            sf=shift_sf, starting_ps=[0.0, 1.0], error=false, dx=[-10, 10],
        )

        # fewer than 3 points in the overlap interval: at p1 = 4 only x = 4 of both sizes
        sp = shift_problem([1, 2])
        @test isfinite(SingleSpline()(sp, [0.0, 1.0]))
        @test SingleSpline()(sp, [4.0, 1.0]) == Inf

        # at p1 = 2 the overlap interval is [4, 4], whose 3 points (one per size) put all
        # knots at x = 4
        sp = shift_problem([1, 2, 3])
        @test isfinite(SingleSpline()(sp, [0.0, 1.0]))
        @test SingleSpline()(sp, [2.0, 1.0]) == Inf
    end

    @testset "MultipleSplines pchip interpolation" begin
        pchip = ScalingCollapse._pchip
        rng = Xoshiro(11)
        xs = sort!(rand(rng, 15))
        dense = collect(range(xs[1], xs[end], length=2001))
        # exact at the data points and for linear data
        ys = sin.(6xs)
        @test pchip(xs, ys, xs) ≈ ys atol = 1e-14
        @test pchip(xs, 2 .* xs .+ 1, dense) ≈ 2 .* dense .+ 1
        # monotone data give a monotone interpolant without overshoot
        ys = cumsum(rand(rng, 15))
        vals = pchip(xs, ys, dense)
        @test all(diff(vals) .>= 0)
        @test extrema(vals) == (ys[1], ys[end])
        # a peak is cut off at the highest data point, and flat stretches stay flat
        ys = exp.(-(xs .- 0.5) .^ 2 ./ 0.02)
        @test maximum(pchip(xs, ys, sort!(vcat(dense, xs)))) == maximum(ys)
        @test pchip([1.0, 2, 3, 4, 5], [0.0, 0, 0, 1, 1], collect(1:0.01:3)) == zeros(201)
        # two points: linear
        @test pchip([0.0, 2.0], [1.0, 5.0], [0.0, 0.5, 2.0]) ≈ [1.0, 2.0, 5.0]
        # converges (at second order) for smooth data
        f(x) = exp(-x^2)
        err(n) = (x = collect(range(-2.5, 2.5, length=n)); z = collect(range(-2.5, 2.5, length=1001));
            maximum(abs, pchip(x, f.(x), z) .- f.(z)))
        @test err(64) < err(16) / 10 && err(64) < 2e-3

        @test_throws ArgumentError MultipleSplines(interp=:linear)
        @test MultipleSplines().interp === :pchip
    end

    @testset "MultipleSplines pchip on an irregular grid" begin
        # noisy data at random x values: the cubic spline oscillates between close points,
        # pchip does not
        rng = Xoshiro(2)
        Ls = [16, 32, 64, 128]
        xs = [sort!(0.7 .* (1 .+ (5 .* rand(rng, 30) .- 2.5) .* L^(-1 / 1.2))) for L in Ls]
        ys = [1 .+ tanh.((x .- 0.7) ./ 0.7 .* L^(1 / 1.2)) for (x, L) in zip(xs, Ls)]
        ys = [y .* (1 .+ 0.01 .* randn(rng, length(y))) for y in ys]
        es = [0.01 .* abs.(y) for y in ys]
        sp = ScalingProblem(xs, ys, es, Ls; sf=ScalingFunction(:x),
            p_space=[0.56:0.01:0.88, 0.3:0.05:3.0])
        @test sp.quality.interp === :pchip
        @test abs(sp.optimal_ps[1] - 0.7) <= 3 * sp.optimal_ps_error[1]
        @test abs(sp.optimal_ps[2] - 1.2) <= 3 * sp.optimal_ps_error[2]
    end
end
