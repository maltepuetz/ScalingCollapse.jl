@testset "ScalingFunction" begin



    # custom scaling function
    function myfunction(d::ScalingCollapse.Data, p1, p2)

        #  initialize arrays for scaled data
        xs = zeros(length(d.xs))
        ys = zeros(length(d.ys))
        es = zeros(length(d.es))

        # scale data according to p1 and p2
        for (i, x, y, e) in zip(eachindex(d.xs), d.xs, d.ys, d.es)
            xs[i] = (x - p1) / p1 * log(d.L)
            ys[i] = y * d.L^(p2)
            es[i] = e * d.L^(p2)
        end

        # create new Data object with scaled data
        return ScalingCollapse.Data(d.L, xs, ys, es)
    end

    sf = ScalingFunction(myfunction; p_names=["myp1", "myp2"])
    @test ScalingCollapse.n_parameters(sf) == 2
    @test ScalingCollapse.scaled_p_names(sf) == ["myp1", "myp2"]
    @test length(ScalingCollapse.fixed_p_names(sf)) == 0

    sf = ScalingFunction(myfunction; N_parameters=2)
    @test ScalingCollapse.n_parameters(sf) == 2
    @test ScalingCollapse.scaled_p_names(sf) == ["p1", "p2"]
    @test length(ScalingCollapse.fixed_p_names(sf)) == 0

    try # error handling
        sf = ScalingFunction(myfunction)
        @test false
    catch
        @test true
    end


    # implicit preset via p_names
    sf = ScalingFunction(["T_c", "nu", "beta"])
    @test ScalingCollapse.n_parameters(sf) == 3
    @test ScalingCollapse.scaled_p_names(sf) == ["T_c", "nu", "beta"]


    try # error handling
        sf = ScalingCollapse.ScalingFunction(:xyxy)
        @test false
    catch
        @test true
    end


    @testset "relative keyword" begin
        d = ScalingCollapse.Data(16, [-0.3, 0.0, 0.2, 0.7], [1.0, 2.0, 3.0, 4.0], [0.1, 0.2, 0.3, 0.4])
        pc, nu, b = 0.4, 0.8, 0.125
        s = d.L^(1 / nu)
        x_rel = (d.xs .- pc) ./ pc .* s
        x_abs = (d.xs .- pc) .* s

        # scaled data match the formula, with and without fixed parameters
        for (kw, x) in (((; relative=true), x_rel), ((; relative=false), x_abs), ((;), x_abs))
            for (sf, ps) in (
                (ScalingFunction(:x; kw...), (pc, nu)),
                (ScalingFunction(:x; p_names=["p_c", "nu"], nu=nu, kw...), (pc,)),
                (ScalingFunction(:x; p_names=["p_c", "nu"], p_c=pc, kw...), (nu,)),
            )
                sd = sf.f(d, ps...)
                @test sd.xs ≈ x
                @test sd.ys == d.ys
                @test sd.es == d.es
            end
            for (preset, sign) in ((:xy, 1), (:xny, -1))
                for (sf, ps) in (
                    (ScalingFunction(preset; kw...), (pc, nu, b)),
                    (ScalingFunction(preset; p_names=["p_c", "nu", "b"], nu=nu, kw...), (pc, b)),
                    (ScalingFunction(preset; p_names=["p_c", "nu", "b"], p_c=pc, b=b, kw...), (nu,)),
                )
                    sd = sf.f(d, ps...)
                    @test sd.xs ≈ x
                    @test sd.ys ≈ d.ys .* d.L^(sign * b / nu)
                    @test sd.es ≈ d.es .* d.L^(sign * b / nu)
                end
            end
        end

        # p_c = 0 is well defined
        @test ScalingFunction(:x).f(d, 0.0, nu).xs ≈ d.xs .* s

        # matching x_scale
        for preset in (:x, :xy, :xny)
            p_names = ["p_c", "nu", "b"][1:(preset == :x ? 2 : 3)]
            @test ScalingFunction(preset; p_names).x_scale == "(x - p_c) * L^(1/nu)"
            @test ScalingFunction(preset; p_names, relative=true).x_scale ==
                  "(x - p_c)/p_c * L^(1/nu)"
        end
        @test ScalingFunction(:x; relative=true, x_scale="custom").x_scale == "custom"

        # error handling
        @test_throws ArgumentError ScalingFunction(:x; relative=0)
        @test_throws ArgumentError ScalingFunction(:x; f=myfunction, relative=true)
        @test_throws ArgumentError ScalingFunction(myfunction; N_parameters=2, relative=true)
    end


    @testset "collapses with p_c = 0 and p_c < 0" begin
        Ls_syn = [8, 16, 32, 64]
        xs_syn = collect(-1.0:0.02:1.0)
        # small deterministic noise, so that the quality minimum is not zero
        noise(x, L) = 1e-3 * sin(997 * x + 13 * L)

        # :x, p_c = 0
        ys_syn = [[tanh(x * L^(1 / 1.0)) + noise(x, L) for x in xs_syn] for L in Ls_syn]
        sp = ScalingProblem(xs_syn, ys_syn, Ls_syn;
            sf=ScalingFunction(:x; p_names=["p_c", "nu"]),
            p_space=[-0.5:0.1:0.5, 0.5:0.1:2.0],
            dx=[-2.0, 2.0],
        )
        @test isapprox(sp.optimal_ps[1], 0.0; atol=1e-3)
        @test isapprox(sp.optimal_ps[2], 1.0; atol=2e-2)
        @test all(isfinite, sp.optimal_ps_error)

        # :xy, p_c < 0
        pc, nu, b = -0.3, 0.8, 0.125
        g(u) = 1 / (1 + exp(u))
        ys_syn = [[L^(-b / nu) * g((x - pc) * L^(1 / nu)) + noise(x, L) for x in xs_syn]
                  for L in Ls_syn]
        sp = ScalingProblem(xs_syn, ys_syn, Ls_syn;
            sf=ScalingFunction(:xy; p_names=["p_c", "nu", "b"]),
            p_space=[-0.5:0.1:0.0, 0.5:0.1:1.5, 0.0:0.05:0.3],
            dx=[-2.0, 2.0],
        )
        @test isapprox(sp.optimal_ps[1], pc; atol=1e-3)
        @test isapprox(sp.optimal_ps[2], nu; atol=2e-2)
        @test isapprox(sp.optimal_ps[3], b; atol=1e-2)
        @test all(isfinite, sp.optimal_ps_error)

        # relative is passed on to the default ScalingFunction: with p_c ∈ [0.2, 0.4] (data
        # mirrored to p_c = 0.3) and the window divided by p_c, the relative form finds the same fit
        sp_abs = ScalingProblem(-xs_syn, ys_syn, Ls_syn;
            p_names=["p_c", "nu", "b"],
            p_space=[0.2:0.1:0.4, 0.5:0.1:1.5, 0.0:0.05:0.3],
            dx=[-2.0, 2.0],
            error=false,
        )
        sp_rel = ScalingProblem(-xs_syn, ys_syn, Ls_syn;
            p_names=["p_c", "nu", "b"], relative=true,
            p_space=[0.2:0.1:0.4, 0.5:0.1:1.5, 0.0:0.05:0.3],
            dx=[-2.0, 2.0] ./ 0.3,
            error=false,
        )
        @test sp_abs.optimal_ps ≈ sp_rel.optimal_ps rtol = 1e-3
        @test sp_rel.sf.x_scale == "(x - p_c)/p_c * L^(1/nu)"
    end

end
