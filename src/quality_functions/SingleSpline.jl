"""
    SingleSpline()

Quality function that interpolates all scaled data points to a single spline. The quality
is then calculated as the sum of the squared differences between the data points and the
interpolated values. The spline is fitted using the KissSmoothing.jl package. The number of
knots can be specified using the `N_knots` keyword argument. If `N_knots` is not specified,
the number of knots is set to `total_number_of_data_points / 10`.
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

    scaled_data = [sp.sf.f(sp.data[i], parameters...) for i in eachindex(sp.data)]

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

    # collect all data points in a single array
    xs = Float64[]
    ys = Float64[]
    es = Float64[]
    for d in scaled_data
        for i in eachindex(d.xs)
            push!(xs, d.xs[i])
            push!(ys, d.ys[i])
            push!(es, d.es[i])
        end
    end

    # sort data points
    perm = sortperm(xs)
    xs = xs[perm]
    ys = ys[perm]
    es = es[perm]


    # set N_knots
    N_knots = sqf.N_knots != 0 ? sqf.N_knots : div(length(xs), 10, RoundUp)


    # fit spline
    knots = LinRange(xs[1], xs[end], N_knots)
    spline = fit_nspline(xs, ys, knots)

    # filter for data points in interval
    mask = findall(x -> x < interval[2] && x > interval[1], xs)

    if sp.errors_defined
        es = es .+ 1e-9  # add small number to avoid division by zero
        S = sum(1 ./ (es[mask] .^ 2) .* (spline.(xs[mask]) .- ys[mask]) .^ 2) / length(mask)
        return S
    else
        S = sum((spline.(xs[mask]) .- ys[mask]) .^ 2) / length(mask)
        return S
    end
end
