"""
    MultipleSplines(; N_steps::Int=100, weight_density::Bool=true, interp::Symbol=:pchip)

Quality function that interpolates the scaled data of every system size separately and
compares the interpolated curves.

# How does it work?
TO BE DOCUMENTED

# Keyword arguments
- `N_steps::Int=100`: Number of points used in the interval to calculate the quality.
- `weight_density::Bool=true`: If `true`, the quality is weighted by the density of data
  points in the interval.
- `interp::Symbol=:pchip`: How the data of each system size are interpolated.
    - `:pchip` (default): monotone piecewise-cubic Hermite interpolation (Fritsch–Carlson).
      It never overshoots between data points, so closely spaced noisy points (e.g. on an
      irregular x grid) do not make the interpolant oscillate. A local extremum of the
      interpolant always lies on a data point, so peaks are cut off at the highest data point.
    - `:cubic`: cubic interpolating spline (Dierckx). It is more accurate for noise-free
      data, especially for master curves with peaks or several extrema, but it oscillates
      for noisy data on irregular grids. It needs at least 4 data points per system size.
"""
struct MultipleSplines <: QualityFunction
    N_steps::Int
    weight_density::Bool
    scan_mode::Bool
    interp::Symbol
    function MultipleSplines(; N_steps::Int=100, weight_density=true, scan_mode=false,
        interp::Symbol=:pchip)
        interp in (:pchip, :cubic) || throw(ArgumentError(
            "interp must be :pchip or :cubic, got :$(interp)."
        ))
        new(N_steps, weight_density, scan_mode, interp)
    end
end

# minimum number of data points per system size
min_points(sqf::MultipleSplines) = sqf.interp === :cubic ? 4 : 2


##### MultipleSplines Quality Function #####
function (sqf::MultipleSplines)(sp, parameters; check_bounds=false)

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

    # check whether the constructed optimization interval is valid (also rejects NaN)
    interval[1] <= interval[2] || return Inf


    # calculate S
    S = 0.0
    if sp.errors_defined && !sqf.scan_mode
        S = _S_weighted(scaled_data, interval, sqf.N_steps, sqf.weight_density, sqf.interp)
    else
        S = _S_unweighted(scaled_data, interval, sqf.N_steps, sqf.interp)
    end

    return S
end

# interpolated y values of one system size at the sorted points xs_eval
function _interpolate(interp::Symbol, xs, ys, xs_eval)
    interp === :cubic && return evaluate(Spline1D(xs, ys, k=3), xs_eval)
    return _pchip(xs, ys, xs_eval)
end

# Monotone piecewise-cubic Hermite interpolation (Fritsch & Carlson 1980), with the slopes of
# Fritsch & Butland (1984) as in Moler's pchip: at an interior point the slope is a weighted
# harmonic mean of the two adjacent secant slopes, or zero if they differ in sign (the point
# is a local extremum of the data). The end slopes are one-sided three-point estimates,
# limited to keep the shape. Outside [xs[1], xs[end]] the end cubics are extrapolated.
function _pchip(xs, ys, xs_eval)
    n = length(xs)
    d = Vector{Float64}(undef, n)
    if n == 2
        d[1] = d[2] = (ys[2] - ys[1]) / (xs[2] - xs[1])
    else
        h_prev = xs[2] - xs[1]
        δ_prev = (ys[2] - ys[1]) / h_prev
        for i in 2:n-1
            h = xs[i+1] - xs[i]
            δ = (ys[i+1] - ys[i]) / h
            if δ_prev * δ > 0
                w1, w2 = 2h + h_prev, h + 2h_prev
                d[i] = (w1 + w2) / (w1 / δ_prev + w2 / δ)
            else
                d[i] = 0.0
            end
            h_prev, δ_prev = h, δ
        end
        h1, h2 = xs[2] - xs[1], xs[3] - xs[2]
        d[1] = _pchip_end_slope(h1, h2, (ys[2] - ys[1]) / h1, (ys[3] - ys[2]) / h2)
        h1, h2 = xs[n] - xs[n-1], xs[n-1] - xs[n-2]
        d[n] = _pchip_end_slope(h1, h2, (ys[n] - ys[n-1]) / h1, (ys[n-1] - ys[n-2]) / h2)
    end

    vals = similar(xs_eval, Float64)
    i = 1
    for (k, x) in enumerate(xs_eval)
        # the interval [xs[i], xs[i+1]] that contains x (xs_eval is sorted)
        while i < n - 1 && x >= xs[i+1]
            i += 1
        end
        h = xs[i+1] - xs[i]
        t = (x - xs[i]) / h
        vals[k] = (2t^3 - 3t^2 + 1) * ys[i] + (t^3 - 2t^2 + t) * h * d[i] +
                  (-2t^3 + 3t^2) * ys[i+1] + (t^3 - t^2) * h * d[i+1]
    end
    return vals
end

function _pchip_end_slope(h1, h2, δ1, δ2)
    s = ((2h1 + h2) * δ1 - h1 * δ2) / (h1 + h2)
    sign(s) != sign(δ1) && return 0.0
    (sign(δ1) != sign(δ2) && abs(s) > abs(3δ1)) && return 3δ1
    return s
end

function _S_unweighted(scaled_data, interval, N_steps, interp)

    xvals = range(interval[1], interval[2], length=N_steps)
    xs_eval = collect(xvals)
    yvals = zeros(N_steps, length(scaled_data))
    for l in eachindex(scaled_data)
        yvals[:, l] = _interpolate(interp, scaled_data[l].xs, scaled_data[l].ys, xs_eval)
    end

    # normalize y, s.t. minimum(y) == 0, maximum(y) == 1
    offset = minimum(yvals)
    for i in eachindex(yvals)
        yvals[i] -= offset
    end
    norm_factor = maximum(yvals)
    for i in eachindex(yvals)
        yvals[i] /= norm_factor
    end

    # calculate quality of the scaling
    S = 0.0
    row = zeros(size(yvals, 2))
    for i in axes(yvals, 1)
        copyto!(row, view(yvals, i, :))
        S += var(row)
    end
    S /= N_steps

    return S
end

function _S_weighted(scaled_data, interval, N_steps, weight_density, interp)

    # linear interpolation of the error bars of each system size
    e_splines = [
        Spline1D(scaled_data[i].xs, scaled_data[i].es, k=1) for i in eachindex(scaled_data)
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
    xs_eval = collect(xvals)
    n = length(scaled_data)
    yvals = Matrix{Float64}(undef, n, N_steps)
    evals = Matrix{Float64}(undef, n, N_steps)  # uncertainties of the spline values
    for l in eachindex(scaled_data)
        yvals[l, :] = _interpolate(interp, scaled_data[l].xs, scaled_data[l].ys, xs_eval)
        # add small number to avoid division by zero
        evals[l, :] = evaluate(e_splines[l], xs_eval) .+ 1e-6
    end

    # normalize y, s.t. minimum(y) == 0, maximum(y) == 1
    offset = minimum(yvals)
    yvals .-= offset
    norm_factor = maximum(yvals)
    yvals ./= norm_factor
    inv_norm = 1 / norm_factor

    # The uncertainties were propagated with Measurements.jl before. The uncertainty of the
    # mean sums the squared uncertainties in Measurements.jl's order, so that the quality
    # does not change in the last digits.
    ord = [1]
    for k in 2:n
        ord = [reverse(ord); k]
    end

    # calculate quality of the scaling
    S = 0.0
    for i in axes(yvals, 2)
        y_sum = yvals[1, i]
        for l in 2:n
            y_sum += yvals[l, i]
        end
        y_mean = y_sum / n
        e2_sum = 0.0
        for l in ord
            e2_sum += abs2(evals[l, i] * inv_norm)
        end
        e_mean = abs(1 / n * sqrt(e2_sum))
        for l in axes(yvals, 1)
            e_l = abs(inv_norm * evals[l, i])
            S += dp_densities[l] * (yvals[l, i] - y_mean)^2 / (e_l^2 + e_mean^2)
        end
    end
    S /= N_steps
    S /= sum(dp_densities)

    return S
end
