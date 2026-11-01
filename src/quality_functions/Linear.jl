
"""
    Linear()

Quality function that interpolates all scaled data points to a single master curve. For
every data point, we take the closest data points, to the left and to the right to
interpolate a linear function. The quality is then calculated as the sum of the squared
differences between the data points and the interpolated values.

# Limitations
The quality is piecewise: it jumps whenever scaled points of different sizes swap order in x
(the neighbour of a third point changes from one size to the other) and whenever a point
enters or leaves the interval. It also dips where the scaled x grids of the sizes line up:
near-coincident points of different sizes are compared with each other without interpolation
error. On uniformly gridded data this gives false minima. For example:
- percolation data on a grid uniform in `(p - 1/2) L^{0.7}`: the quality changes abruptly at
  exactly `p_c = 1/2`, where every size has a point at `x = 0`, and has its lowest profile
  minimum near the commensurate `ν = 1/0.7` instead of near the true `ν = 4/3`;
- 2D Ising Binder cumulant on a uniform temperature grid (`L = 8` to `64`): minimized over
  `T_c`, the quality oscillates by about ±50% every 0.02 in `ν`, with 19 local minima for
  `ν` in `[0.8, 1.1]`.
`Houdayer` and `SingleSpline` show the same kind of ruggedness. For uniformly gridded data,
prefer `MultipleSplines`, which compares smooth interpolations of every size and gives a
smooth quality, and check a fit by scanning the quality around it with [`residuals`](@ref).
"""
struct Linear <: QualityFunction end



##### Linear Quality Function #####
function (::Linear)(sp, parameters; check_bounds=false)

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

    # collect all data points in the interval
    xs = Float64[]
    ys = Float64[]
    es = Float64[]
    for d in scaled_data
        i_l = max(1, searchsortedlast(d.xs, interval[1]) + 1)
        i_r = min(length(d.xs), searchsortedfirst(d.xs, interval[2]) - 1)
        for i in i_l:i_r
            push!(xs, d.xs[i])
            push!(ys, d.ys[i])
            push!(es, d.es[i])
        end
    end

    # the linear interpolation needs at least 3 data points in the interval
    length(xs) < 3 && return Inf

    # sort data points
    perm = sortperm(xs)
    xs = xs[perm]
    ys = ys[perm]
    es = es[perm]

    # calculate quality of the scaling
    S = 0.0
    if sp.errors_defined
        for i in eachindex(xs)[2:end-1]
            # interpolate linearly
            w_l, w_r = interpolation_weights(xs, i)
            ybar = w_l * ys[i-1] + w_r * ys[i+1]
            expected_var = es[i]^2 + w_l^2 * es[i-1]^2 + w_r^2 * es[i+1]^2
            S += (ys[i] - ybar)^2 / expected_var
        end
    else
        for i in eachindex(xs)[2:end-1]
            # interpolate linearly
            w_l, w_r = interpolation_weights(xs, i)
            ybar = w_l * ys[i-1] + w_r * ys[i+1]
            S += (ys[i] - ybar)^2
        end
    end
    S /= (length(xs) - 2)
    return S
end

# weights of the left and right neighbour for the linear interpolation at xs[i]; if both
# neighbours coincide with xs[i] (e.g. several sizes share a scaled x value), the
# interpolation is undefined (0/0) and we take their average instead
function interpolation_weights(xs, i)
    dx = xs[i+1] - xs[i-1]
    dx == 0 && return 0.5, 0.5
    return (xs[i+1] - xs[i]) / dx, (xs[i] - xs[i-1]) / dx
end
