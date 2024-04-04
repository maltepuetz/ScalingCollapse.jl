# CDW transition in the square lattice Holstein model
# (data source: https://owenpb.github.io/FiniteSizeScaling.jl/stable/demo_1/)

# define data
X_L4 = [1.973, 2.989, 3.978, 4.513, 4.754, 4.968, 5.182, 5.476, 5.743, 5.957, 6.225, 6.492, 7.000, 7.989, 8.979]
Y_L4 = [1.250, 1.806, 2.222, 3.056, 3.333, 3.750, 4.306, 4.583, 4.861, 5.278, 5.694, 5.972, 6.250, 6.806, 7.083]
E_L4 = [0.114, 0.123, 0.132, 0.128, 0.135, 0.121, 0.130, 0.101, 0.108, 0.139, 0.133, 0.138, 0.128, 0.199, 0.137]

X_L6 = [2.000, 2.989, 4.005, 4.513, 4.754, 4.995, 5.262, 5.476, 5.743, 6.011, 6.492, 7.000, 7.455, 7.989, 9.005, 9.995, 11.010, 12.000]
Y_L6 = [1.245, 1.667, 2.778, 4.444, 5.694, 6.667, 8.056, 9.306, 10.139, 10.972, 12.222, 13.333, 13.750, 14.167, 14.444, 14.583, 14.722, 14.722]
E_L6 = [0.164, 0.153, 0.152, 0.148, 0.165, 0.151, 0.150, 0.161, 0.138, 0.169, 0.153, 0.148, 0.138, 0.149, 0.167, 0.155, 0.147, 0.160]

X_L8 = [4.487, 4.727, 4.968, 5.235, 5.503, 5.743, 6.011, 6.225, 6.492, 6.733, 6.973, 7.267, 7.481, 7.775, 8.043, 8.524, 8.979, 10.021, 10.983, 12.000]
Y_L8 = [4.861, 5.694, 7.639, 10.278, 13.472, 15.417, 17.917, 19.583, 20.694, 21.528, 22.083, 22.361, 22.778, 23.194, 23.612, 23.889, 24.028, 24.583, 24.722, 24.722]
E_L8 = [0.252, 0.289, 0.243, 0.324, 0.489, 0.524, 0.578, 0.609, 0.575, 0.589, 0.625, 0.621, 0.650, 0.623, 0.621, 0.602, 0.658, 0.665, 0.667, 0.689]

X_L10 = [4.489, 4.757, 5.000, 5.244, 5.487, 5.755, 5.971, 6.293, 6.509, 6.750, 7.018, 7.259, 7.500, 7.768, 8.035, 10.041, 12.020]
Y_L10 = [4.583, 6.111, 9.306, 15.556, 19.861, 21.806, 25.833, 27.639, 31.111, 32.361, 33.333, 33.750, 34.722, 35.972, 35.694, 37.222, 37.500]
E_L10 = [0.254, 0.289, 0.501, 1.189, 1.184, 0.790, 0.804, 1.045, 0.600, 0.584, 0.598, 0.655, 0.680, 0.640, 0.601, 0.575, 0.478]

X_L12 = [4.493, 4.734, 4.981, 5.280, 5.528, 5.637, 5.777, 6.023, 6.268, 6.539, 6.809, 7.024, 7.292]
Y_L12 = [5.432, 6.964, 11.838, 16.713, 23.120, 25.627, 31.058, 36.072, 40.669, 43.454, 45.822, 47.075, 48.050]
E_L12 = [0.321, 0.398, 0.401, 0.556, 0.531, 0.598, 0.601, 0.592, 0.625, 0.630, 0.587, 0.658, 0.688]

@testset "CDW transition (Holstein model)" begin


    sp = ScalingProblem(
        [X_L4, X_L6, X_L8, X_L10, X_L12],
        [Y_L4, Y_L6, Y_L8, Y_L10, Y_L12],
        [E_L4, E_L6, E_L8, E_L10, E_L12],
        [4, 6, 8, 10, 12];
        sf=ScalingCollapse.ScalingFunction(:xny,
            p_names=["x_c", "nu", "gamma"],
            nu=1,
            gamma=1.75
        ),
        p_space=[0.1:0.1:10],
        dx=[-2, 2],
    )
    @test isapprox(sp.optimal_ps[1], 6.0; atol=0.1)
    @test try
        show(sp)
        true
    catch
        false
    end


    sp = ScalingProblem(
        [X_L4, X_L6, X_L8, X_L10, X_L12],
        [Y_L4, Y_L6, Y_L8, Y_L10, Y_L12],
        [E_L4, E_L6, E_L8, E_L10, E_L12],
        [4, 6, 8, 10, 12];
        sf=ScalingCollapse.ScalingFunction(:xny,
            p_names=["x_c", "nu", "gamma"],
            nu=1,
            #gamma=1.75
        ),
        p_space=[5:0.1:7, 1:0.1:3],
        dx=[-2, 2],
    )
    @test isapprox(sp.optimal_ps[1], 6.0; atol=0.2)
    @test isapprox(sp.optimal_ps[2], 1.75; atol=0.1)


    sp = ScalingProblem(
        [X_L4, X_L6, X_L8, X_L10, X_L12],
        [Y_L4, Y_L6, Y_L8, Y_L10, Y_L12],
        [E_L4, E_L6, E_L8, E_L10, E_L12],
        [4, 6, 8, 10, 12];
        sf=ScalingCollapse.ScalingFunction(:xny,
            p_names=["x_c", "nu", "gamma"],
            nu=1,
            #gamma=1.75
        ),
        p_space=[5:0.1:7, 1:0.1:3],
        dx=[-2, 2],
        quality=Houdayer()
    )
    @test isapprox(sp.optimal_ps[1], 6.0; atol=0.2)
    @test isapprox(sp.optimal_ps[2], 1.75; atol=0.1)

    sp = ScalingProblem(
        [X_L4, X_L6, X_L8, X_L10, X_L12],
        [Y_L4, Y_L6, Y_L8, Y_L10, Y_L12],
        [E_L4, E_L6, E_L8, E_L10, E_L12],
        [4, 6, 8, 10, 12];
        sf=ScalingCollapse.ScalingFunction(:xny,
            p_names=["x_c", "nu", "gamma"],
            nu=1,
            #gamma=1.75
        ),
        p_space=[5:0.1:7, 1:0.1:3],
        dx=[-2, 2],
        quality=Spline()
    )
    @test isapprox(sp.optimal_ps[1], 6.0; atol=0.2)
    @test isapprox(sp.optimal_ps[2], 1.75; atol=0.1)

    sp = ScalingProblem(
        [X_L4, X_L6, X_L8, X_L10, X_L12],
        [Y_L4, Y_L6, Y_L8, Y_L10, Y_L12],
        [E_L4, E_L6, E_L8, E_L10, E_L12],
        [4, 6, 8, 10, 12];
        sf=ScalingCollapse.ScalingFunction(:xny,
            p_names=["x_c", "nu", "gamma"],
            nu=1,
            #gamma=1.75
        ),
        p_space=[5:0.1:7, 1:0.1:3],
        dx=[-2, 2],
        quality=SingleSpline()
    )
    @test isapprox(sp.optimal_ps[1], 6.0; atol=0.2)
    @test isapprox(sp.optimal_ps[2], 1.75; atol=0.1)
end
