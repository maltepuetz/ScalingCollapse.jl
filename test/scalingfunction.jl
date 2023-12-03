@testset "ScalingFunction" begin



    # custom scaling function
    function myfunction(d::Scaling.Data, p1, p2)

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
        return Scaling.Data(d.L, xs, ys, es)
    end

    sf = ScalingFunction(myfunction; p_names=["myp1", "myp2"])
    @test Scaling.n_parameters(sf) == 2
    @test Scaling.scaled_p_names(sf) == ["myp1", "myp2"]
    @test length(Scaling.fixed_p_names(sf)) == 0

    sf = ScalingFunction(myfunction; N_parameters=2)
    @test Scaling.n_parameters(sf) == 2
    @test Scaling.scaled_p_names(sf) == ["p1", "p2"]
    @test length(Scaling.fixed_p_names(sf)) == 0

    try # error handling
        sf = ScalingFunction(myfunction)
        @test false
    catch
        @test true
    end


    # implicit preset via p_names
    sf = ScalingFunction(["T_c", "nu", "beta"])
    @test Scaling.n_parameters(sf) == 3
    @test Scaling.scaled_p_names(sf) == ["T_c", "nu", "beta"]


    try # error handling
        sf = Scaling.ScalingFunction(:xyxy)
        @test false
    catch
        @test true
    end

end
