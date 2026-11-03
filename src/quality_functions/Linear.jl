
"""
    Linear()

Quality function that measures the distance of every data point to the line between its
neighbours, with the points of all system sizes merged.

# How it works
All rescaled data points in the window where all sizes overlap (intersected with `dx`) are
sorted by x. Every point `i` except the first and the last is compared with the linear
interpolation `ȳ_i` between its left and right neighbour, which usually belong to other
system sizes. The quality is the mean squared deviation,

    S = 1/(n - 2) Σ_i (y_i - ȳ_i)²,

or, with error bars, the mean of `(y_i - ȳ_i)² / (e_i² + w_l² e_{i-1}² + w_r² e_{i+1}²)`,
with the interpolation weights `w_l` and `w_r` of the neighbours. This was the default
quality function up to version 0.2.

# Limitations
The quality is piecewise: it jumps whenever scaled points of different sizes swap order in x
(the neighbour of a third point changes from one size to the other) and whenever a point
enters or leaves the window. It also dips where the scaled x grids of the sizes line up:
near-coincident points of different sizes are compared with each other without interpolation
error. On uniformly gridded data this gives false minima. For example:
- percolation data on a grid uniform in `(p - 1/2) L^{0.7}`: the quality changes abruptly at
  exactly `p_c = 1/2`, where every size has a point at `x = 0`, and has its lowest profile
  minimum near the commensurate `ν = 1/0.7` instead of near the true `ν = 4/3`;
- 2D Ising Binder cumulant on a uniform temperature grid (`L = 8` to `64`): minimized over
  `T_c`, the quality oscillates by about ±50% every 0.02 in `ν`, with 19 local minima for
  `ν` in `[0.8, 1.1]`.
`Houdayer` and `SingleSpline` show the same kind of ruggedness. For uniformly gridded data,
prefer [`MultipleSplines`](@ref), which compares smooth interpolations of every size and
gives a smooth quality, and check a fit by scanning the quality around it with
[`residuals`](@ref).
"""
struct Linear <: QualityFunction end



##### Linear Quality Function #####
function (qf::Linear)(sp, parameters; check_bounds=false)
    return qf(Workspace(sp), sp, parameters; check_bounds)
end

function (qf::Linear)(ws::Workspace, sp, parameters; check_bounds=false)

    # check if parameters are in bounds -- if not return Inf
    check_bounds && _out_of_bounds(ws, parameters) && return Inf

    scaled_data = _scale_all!(ws, sp.data, parameters)

    # set interval in which we want to optimize
    interval = _overlap_interval(scaled_data, sp.dx)

    # check whether the constructed optimization interval is valid (also rejects NaN)
    interval[1] <= interval[2] || return Inf

    # all data points strictly inside the interval, sorted by x
    ranges = ws.ranges
    for (l, d) in enumerate(scaled_data)
        i_l = max(1, searchsortedlast(d.xs, interval[1]) + 1)
        i_r = min(length(d.xs), searchsortedfirst(d.xs, interval[2]) - 1)
        ranges[l] = i_l:i_r
    end

    # the linear interpolation needs at least 3 data points in the interval
    sum(length, ranges) < 3 && return Inf

    xs, ys, es = ws.xs, ws.ys, ws.es
    n = _merge_sorted!(xs, ys, es, ws.next, scaled_data, ranges)

    # calculate quality of the scaling
    S = 0.0
    if sp.errors_defined
        for i in 2:n-1
            # interpolate linearly
            w_l, w_r = interpolation_weights(xs, i)
            ybar = w_l * ys[i-1] + w_r * ys[i+1]
            expected_var = es[i]^2 + w_l^2 * es[i-1]^2 + w_r^2 * es[i+1]^2
            S += (ys[i] - ybar)^2 / expected_var
        end
    else
        for i in 2:n-1
            # interpolate linearly
            w_l, w_r = interpolation_weights(xs, i)
            ybar = w_l * ys[i-1] + w_r * ys[i+1]
            S += (ys[i] - ybar)^2
        end
    end
    S /= (n - 2)
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
