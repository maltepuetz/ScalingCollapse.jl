"""
    SingleSpline(; N_knots=0)

Quality function that fits one smoothing spline to the rescaled data of all system sizes and
measures the deviation of the data points from it.

# How it works
All rescaled data points in the window where all sizes overlap (intersected with `dx`) are
fitted with a natural cubic least-squares spline (KissSmoothing.jl). Its knots lie at evenly
spaced quantiles of these points, so the spline keeps its flexibility when the scaling
stretches or compresses the x axis. The quality is the mean squared deviation of the points
from the spline, divided by their squared errors if error bars are given.

`N_knots` sets the number of knots; the default `0` uses one knot per 5 data points in the
window (at least 2).

# Limitations
Points entering or leaving the window change the fitted spline abruptly, so on uniformly
gridded data the quality is rugged like that of [`Linear`](@ref).
[`MultipleSplines`](@ref) gives a smooth quality.
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
    check_bounds && _out_of_bounds(sp, parameters) && return Inf

    scaled_data = _scale_all(sp, parameters)

    # set interval in which we want to optimize
    interval = _overlap_interval(scaled_data, sp.dx)

    # check whether the constructed optimization interval is valid (also rejects NaN)
    interval[1] <= interval[2] || return Inf

    # collect all data points in the interval -- fitting the spline to points outside of it
    # would make its flexibility inside the interval depend on the scaling parameters
    ranges = map(scaled_data) do d
        searchsortedfirst(d.xs, interval[1]):searchsortedlast(d.xs, interval[2])
    end

    # the spline needs at least 2 knots and more data points than knots
    sum(length, ranges) < 3 && return Inf

    # sort data points
    xs, ys, es = _merge_sorted(scaled_data, ranges)


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
