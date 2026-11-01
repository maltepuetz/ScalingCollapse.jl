"""
    SingleSpline()

Quality function that interpolates all scaled data points to a single spline. The quality
is then calculated as the sum of the squared differences between the data points and the
interpolated values. The spline is fitted using the KissSmoothing.jl package, only to the
data points in the overlap interval of all system sizes (intersected with `dx`), i.e. the
same points the quality is evaluated on. The knots are placed at evenly spaced quantiles of
these points, so the spline keeps its flexibility when the scaling stretches or compresses
the x axis. The number of knots can be specified using the `N_knots` keyword argument. If
`N_knots` is not specified, the number of knots is set to
`number_of_data_points_in_interval / 5`.

Points entering or leaving the interval change the fitted spline abruptly, so on uniformly
gridded data the quality is rugged like that of `Linear` (see its docstring).
`MultipleSplines` gives a smooth quality.
"""
struct SingleSpline <: QualityFunction
    N_knots::Int
    function SingleSpline(; N_knots::Int=0)
        new(N_knots)
    end
end


##### Single Spline Quality Function #####
function (sqf::SingleSpline)(sp, parameters; check_bounds=false)

    # check if parameters are in bounds -- if not return Inf
    if check_bounds
        for (i, p) in enumerate(parameters)
            if p < sp.p_space[i][1] || p > sp.p_space[i][end]
                return Inf
            end
        end
    end

    scaled_data = Data[sp.sf.f(sp.data[i], parameters...) for i in eachindex(sp.data)]

    # set interval in which we want to optimize
    interval = Vector{Float64}(undef, 2)
    interval[1] = max(
        maximum(scaled_data[i].xs[1] for i in eachindex(scaled_data)),
        sp.dx[1]
    )
    interval[2] = min(
        minimum(scaled_data[i].xs[end] for i in eachindex(scaled_data)),
        sp.dx[2]
    )

    # check whether the constructed optimization interval is valid
    interval[1] > interval[2] && return Inf

    # collect all data points in the interval -- fitting the spline to points outside of it
    # would make its flexibility inside the interval depend on the scaling parameters
    xs = Float64[]
    ys = Float64[]
    es = Float64[]
    for d in scaled_data
        for i in eachindex(d.xs)
            interval[1] <= d.xs[i] <= interval[2] || continue
            push!(xs, d.xs[i])
            push!(ys, d.ys[i])
            push!(es, d.es[i])
        end
    end

    # the spline needs at least 2 knots and more data points than knots
    length(xs) < 3 && return Inf

    # sort data points
    perm = sortperm(xs)
    xs = xs[perm]
    ys = ys[perm]
    es = es[perm]


    # set N_knots
    N_knots = sqf.N_knots != 0 ? sqf.N_knots : div(length(xs), 5, RoundUp)
    N_knots = clamp(N_knots, 2, length(xs) - 1)


    # fit spline with knots at evenly spaced quantiles of the data
    # (broadcast instead of a generator: a closure capturing the reassigned xs would box it)
    knot_idxs = round.(Int, LinRange(0, 1, N_knots) .* (length(xs) - 1)) .+ 1
    knots = unique(xs[knot_idxs])
    length(knots) < 2 && return Inf
    spline = fit_nspline(xs, ys, knots)

    if sp.errors_defined
        es = es .+ 1e-9  # add small number to avoid division by zero
        S = sum(1 ./ (es .^ 2) .* (spline.(xs) .- ys) .^ 2) / length(xs)
        return S
    else
        S = sum((spline.(xs) .- ys) .^ 2) / length(xs)
        return S
    end
end
