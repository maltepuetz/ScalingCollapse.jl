@testset "Utility" begin
    p_space = [1:2]
    res = [[i] for i in 1.0:2.0]
    @test res == ScalingCollapse._parameter_combinations(p_space)

    p_space = [1:2, 10:11]
    res = [[i, j] for i in 1.0:2.0, j in 10.0:11.0]
    @test res == ScalingCollapse._parameter_combinations(p_space)

    p_space = [1:2, 10:11, 100:101]
    res = [[i, j, k] for i in 1.0:2.0, j in 10.0:11.0, k in 100.0:101.0]
    @test res == ScalingCollapse._parameter_combinations(p_space)

    p_space = [1:2, 10:11, 100:101, 1000:1001]
    res = [[i, j, k, l] for i in 1.0:2.0, j in 10.0:11.0, k in 100.0:101.0, l in 1000.0:1001.0]
    @test res == ScalingCollapse._parameter_combinations(p_space)

    p_space = [1:2, 10:11, 100:101, 1000:1001, 10000:10001]
    res = [[i, j, k, l, m] for i in 1.0:2.0, j in 10.0:11.0, k in 100.0:101.0, l in 1000.0:1001.0, m in 10000.0:10001.0]
    @test res == ScalingCollapse._parameter_combinations(p_space)
end
