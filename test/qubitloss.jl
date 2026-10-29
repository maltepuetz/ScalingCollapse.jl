# Coherent information in weakly monitored quantum circuits with qubit loss
# (arXiv:2505.22720, Fig. 3; data source and license: testdata/README.md)

function load_qubitloss(theta_index, Ls)
    h5open("testdata/qubitloss/qubitloss_Ic.h5", "r") do file
        xs = [read(file, "theta$theta_index/L$L/r") for L in Ls]
        ys = [read(file, "theta$theta_index/L$L/Ic") for L in Ls]
        es = [read(file, "theta$theta_index/L$L/Ic_err") for L in Ls]
        return xs, ys, es, Ls
    end
end

# the ScalingProblem call used in the paper
function fit_qubitloss(theta_index, Ls)
    xs, ys, es, Ls = load_qubitloss(theta_index, Ls)
    return ScalingProblem(xs, ys, es, Ls;
        sf=ScalingFunction(:x),
        dx=[-1.25, 1.25],
        p_space=[0.3:0.001:0.6, 1.2:0.1:1.8],
        error_threshold=1.3,
        quality=MultipleSplines(),
    )
end

@testset "Qubit loss (arXiv:2505.22720, Fig. 3)" begin
    @testset "θ = π/2: bond percolation (exact)" begin
        sp = fit_qubitloss(15, [32, 64, 128, 256])
        (rc, nu), (σrc, σnu) = sp.optimal_ps, sp.optimal_ps_error
        @test abs(rc - 0.5) <= 3σrc
        @test abs(nu - 4 / 3) <= 3σnu
        # ceilings at ~3× the error of a converged fit (4.33e-5, 5.92e-3); a scan grid point
        # returned without Nelder-Mead progress has 1.9e-4 and 8.8e-2
        @test σrc <= 1.3e-4
        @test σnu <= 1.8e-2
    end
    @testset "θ = 0: Nishimori (published, v0.2.2 + Optim 1.11)" begin
        sp = fit_qubitloss(1, [32, 64, 128, 256, 512])
        ref = [0.33681399873460816, 1.5316126320454644]
        σref = [2.0879977065466418e-5, 0.004006908264828457]
        @test all(abs.(sp.optimal_ps .- ref) .<= 0.25 .* σref)
        @test all(0.5 .<= sp.optimal_ps_error ./ σref .<= 2)
    end
end
