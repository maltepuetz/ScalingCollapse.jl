@info "Loading testdata..."
file = h5open("testdata/IsingMonteCarlo.h5", "r")
Ts = read(file, "Ts")
Ls = read(file, "Ls")
binder = read(file, "binder")
M_abs_mean = read(file, "M_abs_mean")
susceptibility = read(file, "susceptibility")
close(file)
@info "Done. Starting tests..."

@testset "Ising model 2D" begin
    @testset "Binder cumulant" begin

        # create different input variants
        _Ts = [Ts for _ in eachindex(Ls)]
        _binder = [binder[:, i] for i in eachindex(Ls)]

        sp = ScalingCollapse.ScalingProblem(_Ts, _binder, Ls;
            sf=ScalingCollapse.ScalingFunction(:x),
            p_space=[0.1:0.1:3, 0.1:0.1:3],
            dx=[-1.0, 1.0],
        )
        @test isapprox(sp.optimal_ps[1], 2.269; atol=0.01)
        @test isapprox(sp.optimal_ps[2], 1.0; atol=0.05)


        # test with SingleSpline quality function
        sp = ScalingCollapse.ScalingProblem(_Ts, _binder, Ls;
            sf=ScalingCollapse.ScalingFunction(:x),
            p_space=[0.1:0.1:3, 0.1:0.1:3],
            dx=[-1.0, 1.0],
            quality=SingleSpline(),
        )
        @test isapprox(sp.optimal_ps[1], 2.269; atol=0.01)
        @test isapprox(sp.optimal_ps[2], 1.0; atol=0.05)


        # create different input variants
        _Ts = zeros(size(binder))
        for i in axes(_Ts, 2), j in axes(_Ts, 1)
            _Ts[j, i] = Ts[j]
        end
        _binder = binder

        sp = ScalingCollapse.ScalingProblem(_Ts, _binder, Ls;
            sf=ScalingCollapse.ScalingFunction(:x, nu=1, p_names=["T_c", "nu"]),
            dx=[-1.0, 1.0],
            error=false,
        )
        @test isapprox(sp.optimal_ps[1], 2.269; atol=0.001)
        @test try
            show(sp)
            true
        catch
            false
        end

        sp = ScalingCollapse.ScalingProblem(Ts, binder, Ls;
            sf=ScalingCollapse.ScalingFunction(:x; p_names=["T_c", "nu"]),
            dx=[-1.0, 1.0],
            starting_ps=[2.27, 1.0]
        )
        @test isapprox(sp.optimal_ps[1], 2.269; atol=0.01)
        @test isapprox(sp.optimal_ps[2], 1.0; atol=0.05)

        sp = ScalingCollapse.ScalingProblem(Ts, binder, Ls;
            sf=ScalingCollapse.ScalingFunction(:x; p_names=["T_c", "nu"]),
            dx=[-1.0, 1.0],
            starting_ps=[2.27, 1.0],
            error=true
        )
        @test isapprox(sp.optimal_ps[1], 2.269; atol=0.01)
        @test isapprox(sp.optimal_ps[2], 1.0; atol=0.05)

        @test try
            show(sp)
            true
        catch
            false
        end
    end

    @testset "Magnetization" begin

        sp = ScalingCollapse.ScalingProblem(Ts, M_abs_mean, Ls;
            dx=[-1.0, 1.0],
        )
        @test isapprox(sp.optimal_ps[1], 2.269; atol=0.01)
        @test isapprox(sp.optimal_ps[2], 1.0; atol=0.1)
        @test isapprox(sp.optimal_ps[3], 0.125; atol=0.015)

        sp = ScalingCollapse.ScalingProblem(Ts, M_abs_mean, Ls;
            sf=ScalingCollapse.ScalingFunction(:xy, beta=0.125, p_names=["T_c", "nu", "beta"]),
            dx=[-1.0, 1.0],
        )
        @test isapprox(sp.optimal_ps[1], 2.269; atol=0.01)
        @test isapprox(sp.optimal_ps[2], 1.0; atol=0.1)

        sp = ScalingCollapse.ScalingProblem(Ts, M_abs_mean, Ls;
            sf=ScalingCollapse.ScalingFunction(:xy, beta=0.125, p_names=["T_c", "nu", "beta"]),
            dx=[-1.0, 1.0],
            starting_ps=[2.26, 1.0],
            qualtiy=MultipleSplines()
        )
        @test isapprox(sp.optimal_ps[1], 2.269; atol=0.01)
        @test isapprox(sp.optimal_ps[2], 1.0; atol=0.1)

        # test with SingleSpline quality function
        sp = ScalingCollapse.ScalingProblem(Ts, M_abs_mean, Ls;
            dx=[-1.0, 1.0],
            quality=SingleSpline()
        )
        @test isapprox(sp.optimal_ps[1], 2.269; atol=0.01)
        @test isapprox(sp.optimal_ps[2], 1.0; atol=0.1)
        @test isapprox(sp.optimal_ps[3], 0.125; atol=0.015)
    end

    @testset "Susceptibility" begin

        sp = ScalingCollapse.ScalingProblem(Ts, susceptibility, Ls;
            sf=ScalingCollapse.ScalingFunction(:xny),
            dx=[-1.0, 1.0],
        )
        @test isapprox(sp.optimal_ps[1], 2.269; atol=0.01)
        @test isapprox(sp.optimal_ps[2], 1.0; atol=0.1)
        @test isapprox(sp.optimal_ps[3], 1.75; atol=0.2)

        sp = ScalingCollapse.ScalingProblem(Ts, susceptibility, Ls;
            sf=ScalingCollapse.ScalingFunction(:xny,
                gamma=1.75,
                p_names=["T_c", "nu", "gamma"]
            ),
            dx=[-1.0, 1.0],
        )
        @test isapprox(sp.optimal_ps[1], 2.269; atol=0.02)
        @test isapprox(sp.optimal_ps[2], 1.0; atol=0.02)

        sp = ScalingCollapse.ScalingProblem(Ts, susceptibility, Ls;
            sf=ScalingCollapse.ScalingFunction(:xny,
                gamma=1.75,
                p_names=["T_c", "nu", "gamma"]
            ),
            dx=[-1.0, 1.0],
            starting_ps=[2.27, 1.0]
        )
        @test isapprox(sp.optimal_ps[1], 2.269; atol=0.02)
        @test isapprox(sp.optimal_ps[2], 1.0; atol=0.02)

        # test with SingleSpline quality function
        sp = ScalingCollapse.ScalingProblem(Ts, susceptibility, Ls;
            sf=ScalingCollapse.ScalingFunction(:xny),
            dx=[-1.0, 1.0],
            quality=SingleSpline()
        )
        @test isapprox(sp.optimal_ps[1], 2.269; atol=0.01)
        @test isapprox(sp.optimal_ps[2], 1.0; atol=0.1)
        @test isapprox(sp.optimal_ps[3], 1.75; atol=0.2)

        # residual landscape
        sp = ScalingCollapse.ScalingProblem(Ts, susceptibility, Ls;
            sf=ScalingCollapse.ScalingFunction(:xny, p_names=["T_c", "nu", "gamma"]),
            p_space=[2:0.1:2.5, 0.5:0.1:1.5, 1.5:0.1:2.5],
            dx=[-1.0, 2.2],
        )
        p_space, residuals = ScalingCollapse.residuals(sp; dims=[1, 2, 3], N_steps=5)
        ind = argmin(residuals)
        @test p_space[1][ind[1]] == sp.optimal_ps[1]
        @test p_space[2][ind[2]] == sp.optimal_ps[2]
        @test p_space[3][ind[3]] == sp.optimal_ps[3]

        p_space, residuals = ScalingCollapse.residuals(sp; dims=[1, 2], N_steps=5)
        ind = argmin(residuals)
        @test p_space[1][ind[1]] == sp.optimal_ps[1]
        @test p_space[2][ind[2]] == sp.optimal_ps[2]

        p_space, residuals = ScalingCollapse.residuals(sp; dims=[1], N_steps=5)
        ind = argmin(residuals)
        @test p_space[1][ind[1]] == sp.optimal_ps[1]

        try # test error handling
            p_space, residuals = ScalingCollapse.residuals(sp; dims=[1, 2, 4], N_steps=5)
            @test false
        catch
            @test true
        end
        try # test error handling
            p_space, residuals = ScalingCollapse.residuals(sp; dims=[1, 2, 2], N_steps=5)
            @test false
        catch
            @test true
        end
        try # test error handling
            p_space, residuals = ScalingCollapse.residuals(sp; dims=[], N_steps=5)
            @test false
        catch
            @test true
        end


        sp = ScalingCollapse.ScalingProblem(Ts, susceptibility, Ls;
            sf=ScalingCollapse.ScalingFunction(:xny, p_names=["T_c", "nu", "gamma"]),
            p_space=[2:0.1:2.5, 0.5:0.1:1.5, 1.5:0.1:2.5],
            dx=[-1.0, 2.2],
            quality=MultipleSplines()
        )
        # export scaled data
        sx, sy, se, sL = ScalingCollapse.scaled_data(sp)
        @test length(sx) == length(sL)
        @test length(sy) == length(sL)
        @test length(se) == length(sL)
        @test sL == Ls

        # export scaled data with splines
        sx, sy, se, sL, xspl, yspl, espl = ScalingCollapse.scaled_data(sp; splines=true)
        @test length(sx) == length(sL)
        @test length(sy) == length(sL)
        @test length(se) == length(sL)
        @test sL == Ls
        @test length(xspl) == length(sL)
        @test length(yspl) == length(sL)
        @test length(espl) == length(sL)

    end
end
