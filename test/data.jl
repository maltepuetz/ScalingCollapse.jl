@testset "Data read-in" begin

    # with errors and "correct ordered"
    Ls = rand(1:1000, 3)
    xs = randn(10)
    ys = randn(10, length(Ls))
    es = randn(10, length(Ls))

    d1 = ScalingCollapse.unzip_data(xs, ys, es, Ls)
    d2 = ScalingCollapse.unzip_data([xs xs xs], ys, es, Ls)
    d3 = ScalingCollapse.unzip_data(
        [xs for _ in eachindex(Ls)],
        [ys[:, l] for l in eachindex(Ls)],
        [es[:, l] for l in eachindex(Ls)],
        Ls
    )
    d4 = ScalingCollapse.unzip_data(xs, ScalingCollapse.measurement.(ys, es), Ls)
    d5 = ScalingCollapse.unzip_data([xs xs xs], ScalingCollapse.measurement.(ys, es), Ls)
    d6 = ScalingCollapse.unzip_data(
        [xs for _ in eachindex(Ls)],
        [ScalingCollapse.measurement.(ys[:, l], es[:, l]) for l in eachindex(Ls)],
        Ls
    )
    @test isequal(d1, d2)
    @test isequal(d1, d3)
    @test isequal(d1, d4)
    @test isequal(d1, d5)
    @test isequal(d1, d6)

    # with errors and "not-correct ordered"
    Ls = rand(1:1000, 3)
    xs = randn(10)
    ys = randn(length(Ls), 10)
    es = randn(length(Ls), 10)

    d1 = ScalingCollapse.unzip_data(xs, ys, es, Ls)
    d2 = ScalingCollapse.unzip_data([xs[j] for i in axes(ys, 1), j in axes(ys, 2)], ys, es, Ls)
    d3 = ScalingCollapse.unzip_data(
        [xs for _ in eachindex(Ls)],
        [ys[l, :] for l in eachindex(Ls)],
        [es[l, :] for l in eachindex(Ls)],
        Ls
    )
    d4 = ScalingCollapse.unzip_data(xs, ScalingCollapse.measurement.(ys, es), Ls)
    d5 = ScalingCollapse.unzip_data(
        [xs[j] for i in axes(ys, 1), j in axes(ys, 2)],
        ScalingCollapse.measurement.(ys, es),
        Ls
    )
    d6 = ScalingCollapse.unzip_data(
        [xs for _ in eachindex(Ls)],
        [ScalingCollapse.measurement.(ys[l, :], es[l, :]) for l in eachindex(Ls)],
        Ls
    )
    @test isequal(d1, d2)
    @test isequal(d1, d3)
    @test isequal(d1, d4)
    @test isequal(d1, d5)
    @test isequal(d1, d6)


    # without errors and "not-correct ordered"
    Ls = rand(1:1000, 3)
    xs = randn(10)
    ys = randn(length(Ls), 10)

    d1 = ScalingCollapse.unzip_data(xs, ys, Ls)
    d2 = ScalingCollapse.unzip_data([xs[j] for i in axes(ys, 1), j in axes(ys, 2)], ys, Ls)
    d3 = ScalingCollapse.unzip_data(
        [xs for _ in eachindex(Ls)],
        [ys[l, :] for l in eachindex(Ls)],
        Ls
    )
    @test isequal(d1, d2)
    @test isequal(d1, d3)
end
