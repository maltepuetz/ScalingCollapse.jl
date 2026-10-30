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
end
