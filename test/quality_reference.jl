# The quality functions were optimized for speed without changing their results. These tests
# compare them bit for bit with the previous implementations and guard the type stability that
# the speed depends on.
using Random
using Statistics

const _Dierckx = ScalingCollapse.Dierckx
const _Measurements = ScalingCollapse.Measurements

# previous implementation of MultipleSplines without error bars
function _S_unweighted_reference(scaled_data, interval, N_steps)
    y_splines = [
        _Dierckx.Spline1D(scaled_data[i].xs, scaled_data[i].ys, k=3) for i in eachindex(scaled_data)
    ]
    xvals = range(interval[1], interval[2], length=N_steps)
    yvals = zeros(N_steps, length(scaled_data))
    for l in eachindex(scaled_data)
        for (i, x) in enumerate(xvals)
            yvals[i, l] = y_splines[l](x)
        end
    end
    offset = minimum(yvals)
    for i in eachindex(yvals)
        yvals[i] -= offset
    end
    norm_factor = maximum(yvals)
    for i in eachindex(yvals)
        yvals[i] /= norm_factor
    end
    S = 0.0
    for i in axes(yvals, 1)
        S += var(yvals[i, :])
    end
    S /= N_steps
    return S
end

# previous implementation of MultipleSplines with error bars, which propagated the
# uncertainties with Measurements.jl
function _S_weighted_reference(scaled_data, interval, N_steps, weight_density)
    y_splines = [
        _Dierckx.Spline1D(scaled_data[i].xs, scaled_data[i].ys, k=3) for i in eachindex(scaled_data)
    ]
    e_splines = [
        _Dierckx.Spline1D(scaled_data[i].xs, scaled_data[i].es, k=1) for i in eachindex(scaled_data)
    ]
    dp_densities = zeros(length(scaled_data))
    for l in eachindex(scaled_data)
        i_l = searchsortedlast(scaled_data[l].xs, interval[1])
        i_r = searchsortedfirst(scaled_data[l].xs, interval[2])
        x_l = scaled_data[l].xs[i_l]
        x_r = scaled_data[l].xs[i_r]
        dp_densities[l] = (i_r - i_l - 1) / (x_r - x_l)
    end
    !weight_density && (dp_densities .= 1.0)
    xvals = range(interval[1], interval[2], length=N_steps)
    yvals = zeros(_Measurements.Measurement, length(scaled_data), N_steps)
    for (i, x) in enumerate(xvals)
        for l in eachindex(scaled_data)
            y = y_splines[l](x)
            e = e_splines[l](x) + 1e-6
            yvals[l, i] = _Measurements.measurement(y, e)
        end
    end
    offset = minimum(yvals).val
    for i in eachindex(yvals)
        yvals[i] -= offset
    end
    norm_factor = maximum(yvals).val
    for i in eachindex(yvals)
        yvals[i] /= norm_factor
    end
    S = 0.0
    for i in axes(yvals, 2)
        for l in axes(yvals, 1)
            y_mean = mean(yvals[:, i])
            S += dp_densities[l] * (yvals[l, i].val - y_mean.val)^2 / (yvals[l, i].err^2 + y_mean.err^2)
        end
    end
    S /= N_steps
    S /= sum(dp_densities)
    return S
end

# random scaled data of n_sizes sizes on x in about [-3, 3] with overlap interval [-1.5, 1.5]
function _random_scaled_data(rng, n_sizes; n_points=40, zero_errors=false)
    map(1:n_sizes) do l
        xs = sort!(6 .* rand(rng, n_points) .- 3)
        ys = tanh.(xs .+ 0.1 .* randn(rng)) .* (1 + 0.1 * randn(rng)) .+ 0.01 .* randn(rng, n_points)
        es = zero_errors ? zeros(n_points) : 0.05 .* rand(rng, n_points)
        ScalingCollapse.Data(8 * l, xs, ys, es)
    end
end

# number of allocations of q(sp, ps) after compilation
function _allocations(q, sp, ps)
    q(sp, ps)
    return @allocations q(sp, ps)
end

@testset "Quality functions match the previous implementations" begin

    @testset "MultipleSplines bit for bit" begin
        rng = Xoshiro(20261031)
        interval = [-1.5, 1.5]
        for n_sizes in 2:12, trial in 1:20
            scaled_data = _random_scaled_data(rng, n_sizes)
            N_steps = rand(rng, (10, 57, 100))
            @test isequal(
                ScalingCollapse._S_unweighted(scaled_data, interval, N_steps),
                _S_unweighted_reference(scaled_data, interval, N_steps),
            )
            for weight_density in (true, false)
                @test isequal(
                    ScalingCollapse._S_weighted(scaled_data, interval, N_steps, weight_density),
                    _S_weighted_reference(scaled_data, interval, N_steps, weight_density),
                )
            end
        end
        # zero error bars (only the 1e-6 regularization is left) and identical sizes (S = 0)
        for scaled_data in (
            _random_scaled_data(rng, 4; zero_errors=true),
            fill(first(_random_scaled_data(rng, 1)), 3),
        )
            @test isequal(
                ScalingCollapse._S_weighted(scaled_data, interval, 100, true),
                _S_weighted_reference(scaled_data, interval, 100, true),
            )
            @test isequal(
                ScalingCollapse._S_unweighted(scaled_data, interval, 100),
                _S_unweighted_reference(scaled_data, interval, 100),
            )
        end
    end

    @testset "MultipleSplines through a ScalingProblem" begin
        rng = Xoshiro(7)
        Ls = [8, 16, 32, 64]
        xs = [collect(range(0.3, 0.7, length=60)) for _ in Ls]
        ys = [tanh.((x .- 0.5) ./ 0.5 .* L) .* L^0.1 .+ 0.01 .* randn(rng, 60) for (x, L) in zip(xs, Ls)]
        es = [fill(0.01, 60) for _ in Ls]
        sp = ScalingProblem(xs, ys, es, Ls; sf=ScalingFunction(:xy), quality=MultipleSplines(),
            starting_ps=[0.5, 1.0, 0.1], error=false)
        for ps in ([0.5, 1.0, 0.1], [0.49, 1.1, 0.12], [0.51, 0.9, 0.05])
            scaled = [sp.sf.f(d, ps...) for d in sp.data]
            interval = [maximum(d.xs[1] for d in scaled), minimum(d.xs[end] for d in scaled)]
            @test isequal(MultipleSplines()(sp, ps), _S_weighted_reference(scaled, interval, 100, true))
            @test isequal(MultipleSplines(scan_mode=true)(sp, ps),
                _S_unweighted_reference(scaled, interval, 100))
        end
    end

    @testset "Data sorting" begin
        rng = Xoshiro(3)
        xs = sort!(rand(rng, 50))
        ys, es = randn(rng, 50), randn(rng, 50)
        d = ScalingCollapse.Data(4, xs, ys, es)
        @test d.xs == xs && d.ys == ys && d.es == abs.(es)
        # the same result from shuffled input
        perm = randperm(rng, 50)
        @test isequal(ScalingCollapse.Data(4, xs[perm], ys[perm], es[perm]), d)
        # the data are copied, also when the input is already sorted
        xs[1] = -1.0
        ys[1] = 42.0
        @test d.xs[1] != -1.0 && d.ys[1] != 42.0
        # equal x values keep their order (stable sort), sorted or not
        d = ScalingCollapse.Data(4, [1.0, 2.0, 2.0, 3.0], [1.0, 2.0, 3.0, 4.0], zeros(4))
        @test d.ys == [1.0, 2.0, 3.0, 4.0]
        d = ScalingCollapse.Data(4, [3.0, 2.0, 2.0, 1.0], [4.0, 2.0, 3.0, 1.0], zeros(4))
        @test d.xs == [1.0, 2.0, 2.0, 3.0] && d.ys == [1.0, 2.0, 3.0, 4.0]
    end

    @testset "Quality functions are type stable" begin
        # a type instability in the loops over the rescaled data allocates for every data
        # point, so the number of allocations would grow with the number of points
        function problem(n_points)
            Ls = [8, 16, 32]
            xs = [collect(range(0.3, 0.7, length=n_points)) for _ in Ls]
            ys = [tanh.((x .- 0.5) ./ 0.5 .* L) for (x, L) in zip(xs, Ls)]
            es = [fill(0.01, n_points) for _ in Ls]
            return ScalingProblem(xs, ys, es, Ls; sf=ScalingFunction(:x),
                starting_ps=[0.5, 1.0], error=false)
        end
        small, large = problem(50), problem(500)
        # (a fixed number of knots: the spline fit itself allocates more for more knots)
        for q in (Linear(), SingleSpline(N_knots=10), MultipleSplines(),
            MultipleSplines(scan_mode=true))
            n_small = _allocations(q, small, [0.5, 1.0])
            n_large = _allocations(q, large, [0.5, 1.0])
            @test n_large < 2 * n_small
        end
        # with error bars, MultipleSplines propagated the uncertainties with Measurements.jl,
        # which allocated about 13000 times here (now about 200)
        @test _allocations(MultipleSplines(), small, [0.5, 1.0]) < 2000
    end
end
