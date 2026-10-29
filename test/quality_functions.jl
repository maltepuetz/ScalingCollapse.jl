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
                    sf=ScalingFunction(:x), starting_ps=[0.5, 1.0], error=false,
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
            sf=ScalingFunction((d, p1, p2) -> d; N_parameters=2),
            starting_ps=[1.0, 1.0], error=false, dx=[-10, 10],
        )
        # interior points: ybar(3) = 10, ybar(4) = 17 -> S = (1 + 1) / 2
        @test sp.quality(sp, [1.0, 1.0]) ≈ 1.0
    end
end
