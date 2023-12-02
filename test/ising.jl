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

        sp = Scaling.ScalingProblem(_Ts, _binder, Ls;
            sf=Scaling.ScalingFunction(:x),
            p_space=[0.1:0.1:3, 0.1:0.1:3],
            dx=[-1.0, 1.0],
        )
        @test isapprox(sp.optimal_ps[1], 2.269; atol=0.01)
        @test isapprox(sp.optimal_ps[2], 1.0; atol=0.05)


        # create different input variants
        _Ts = zeros(size(binder))
        for i in axes(_Ts, 2), j in axes(_Ts, 1)
            _Ts[j, i] = Ts[j]
        end
        _binder = binder

        sp = Scaling.ScalingProblem(_Ts, _binder, Ls;
            sf=Scaling.ScalingFunction(:x, nu=1, p_names=["T_c", "nu"]),
            dx=[-1.0, 1.0],
        )
        @test isapprox(sp.optimal_ps[1], 2.269; atol=0.001)

        sp = Scaling.ScalingProblem(Ts, binder, Ls;
            sf=Scaling.ScalingFunction(:x; p_names=["T_c", "nu"]),
            dx=[-1.0, 1.0],
            starting_ps=[2.27, 1.0]
        )
        @test isapprox(sp.optimal_ps[1], 2.269; atol=0.01)
        @test isapprox(sp.optimal_ps[2], 1.0; atol=0.05)

        sp = Scaling.ScalingProblem(Ts, binder, Ls;
            sf=Scaling.ScalingFunction(:x; p_names=["T_c", "nu"]),
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

        sp = Scaling.ScalingProblem(Ts, M_abs_mean, Ls;
            dx=[-1.0, 1.0],
        )
        @test isapprox(sp.optimal_ps[1], 2.269; atol=0.01)
        @test isapprox(sp.optimal_ps[2], 1.0; atol=0.1)
        @test isapprox(sp.optimal_ps[3], 0.125; atol=0.015)

        sp = Scaling.ScalingProblem(Ts, M_abs_mean, Ls;
            sf=Scaling.ScalingFunction(:xy, beta=0.125, p_names=["T_c", "nu", "beta"]),
            dx=[-1.0, 1.0],
        )
        @test isapprox(sp.optimal_ps[1], 2.269; atol=0.01)
        @test isapprox(sp.optimal_ps[2], 1.0; atol=0.1)

        sp = Scaling.ScalingProblem(Ts, M_abs_mean, Ls;
            sf=Scaling.ScalingFunction(:xy, beta=0.125, p_names=["T_c", "nu", "beta"]),
            dx=[-1.0, 1.0],
            starting_ps=[2.2, 1.0]
        )
        @test isapprox(sp.optimal_ps[1], 2.269; atol=0.01)
        @test isapprox(sp.optimal_ps[2], 1.0; atol=0.1)
    end

    @testset "Susceptibility" begin

        sp = Scaling.ScalingProblem(Ts, susceptibility, Ls;
            sf=Scaling.ScalingFunction(:xny),
            dx=[-1.0, 1.0],
        )
        @test isapprox(sp.optimal_ps[1], 2.269; atol=0.01)
        @test isapprox(sp.optimal_ps[2], 1.0; atol=0.1)
        @test isapprox(sp.optimal_ps[3], 1.75; atol=0.2)

        sp = Scaling.ScalingProblem(Ts, susceptibility, Ls;
            sf=Scaling.ScalingFunction(:xny, gamma=1.75, p_names=["T_c", "nu", "gamma"]),
            dx=[-1.0, 1.0],
        )
        @test isapprox(sp.optimal_ps[1], 2.269; atol=0.01)
        @test isapprox(sp.optimal_ps[2], 1.0; atol=0.01)

        sp = Scaling.ScalingProblem(Ts, susceptibility, Ls;
            sf=Scaling.ScalingFunction(:xny, gamma=1.75, p_names=["T_c", "nu", "gamma"]),
            dx=[-1.0, 1.0],
            starting_ps=[2.27, 1.0]
        )
        @test isapprox(sp.optimal_ps[1], 2.269; atol=0.01)
        @test isapprox(sp.optimal_ps[2], 1.0; atol=0.01)
    end
end
