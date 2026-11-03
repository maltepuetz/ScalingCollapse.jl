"""
    MultipleSplines(; interp=:pchip, N_steps=100, weight_density=true, scan_mode=false)

The default quality function. It interpolates the rescaled data of every system size
separately and measures how much the interpolated curves differ.

# How it works
The curves are compared at `N_steps` evenly spaced points `x_k` in the window where all
rescaled sizes overlap (intersected with `dx`). With the interpolated values `y_l(x_k)` of
system size `l`, normalized so that all of them lie between 0 and 1 (and the errors
scaled alike), the quality is
- without error bars: the variance across system sizes, averaged over the points,

      S = 1/N_steps Σ_k var_l(y_l(x_k));

- with error bars: the error-weighted squared deviation from the mean over system sizes,

      S = Σ_k Σ_l ρ_l (y_l(x_k) - ȳ(x_k))² / (e_l(x_k)² + ē(x_k)²) / (N_steps Σ_l ρ_l),

    where `e_l` is the linearly interpolated error of size `l`, `ē = sqrt(Σ_l e_l²) / n`
    the error of the mean `ȳ` of the `n` sizes, and `ρ_l` the number of data points of size
    `l` per unit of rescaled x in the window (or 1 with `weight_density=false`).

Since the curves are compared on a fixed set of points, the quality changes smoothly with
the parameters, also for data on uniform grids (see `Linear` for the alternative). The
normalization makes the quality independent of the scale of y, which matters when y is
rescaled with a power of `L`.

# Keyword arguments
- `interp::Symbol=:pchip`: the interpolation of every system size.
    - `:pchip` (default): monotone piecewise-cubic Hermite interpolation (Fritsch–Carlson).
        It never overshoots between data points, so noisy points on irregular x grids do
        not make the curve oscillate. A local extremum of the curve always lies on a data
        point, so a peak is cut off at its highest data point.
    - `:cubic`: cubic interpolating spline (Dierckx.jl). More accurate for noise-free data,
        especially for master curves with peaks or several extrema, but it oscillates for
        noisy data on irregular grids. Needs at least 4 data points per system size.
- `N_steps::Int=100`: the number of points at which the curves are compared.
- `weight_density::Bool=true`: with error bars, weight every system size by its density
    of data points in the window, so that sizes with many points count more.
- `scan_mode::Bool=false`: ignore the error bars (always use the unweighted form). This
    is the setting of the default `quality_scan` of [`ScalingProblem`](@ref).
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
    check_bounds && _out_of_bounds(sp, parameters) && return Inf

    scaled_data = _scale_all(sp, parameters)

    # set interval in which we want to optimize
    interval = _overlap_interval(scaled_data, sp.dx)

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

# write the interpolated y values of one system size at the sorted points xs_eval into vals;
# slopes is a buffer of at least length(xs) for pchip
function _interpolate!(vals, slopes, interp::Symbol, xs, ys, xs_eval)
    if interp === :cubic
        copyto!(vals, evaluate(Spline1D(xs, ys, k=3), collect(xs_eval)))
    else
        _pchip!(vals, slopes, xs, ys, xs_eval)
    end
    return vals
end

# Monotone piecewise-cubic Hermite interpolation (Fritsch & Carlson 1980), with the slopes of
# Fritsch & Butland (1984) as in Moler's pchip: at an interior point the slope is a weighted
# harmonic mean of the two adjacent secant slopes, or zero if they differ in sign (the point
# is a local extremum of the data). The end slopes are one-sided three-point estimates,
# limited to keep the shape. Outside [xs[1], xs[end]] the end cubics are extrapolated.
function _pchip(xs, ys, xs_eval)
    vals = Vector{Float64}(undef, length(xs_eval))
    return _pchip!(vals, Vector{Float64}(undef, length(xs)), xs, ys, xs_eval)
end

# _pchip writing into vals, with the slopes d in a buffer of at least length(xs)
function _pchip!(vals, d, xs, ys, xs_eval)
    n = length(xs)
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

    i = 1
    for (k, x) in zip(eachindex(vals), xs_eval)
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

# write the linear interpolation of (xs, ys) at the sorted points xs_eval, which lie in
# [xs[1], xs[end]], into vals. The arithmetic is that of Dierckx's evaluation of a linear
# B-spline (fpbspl), which this replaces, so the values agree bit for bit.
function _linear_interpolation!(vals, xs, ys, xs_eval)
    n = length(xs)
    i = 1
    for (k, x) in zip(eachindex(vals), xs_eval)
        while i < n - 1 && x >= xs[i+1]
            i += 1
        end
        f = 1 / (xs[i+1] - xs[i])
        vals[k] = ys[i] * (f * (xs[i+1] - x)) + ys[i+1] * (f * (x - xs[i]))
    end
    return vals
end

# a buffer for the pchip slopes of every size of scaled_data
_slopes_buffer(scaled_data) = Vector{Float64}(undef, maximum(d -> length(d.xs), scaled_data))

function _S_unweighted(scaled_data, interval, N_steps, interp)

    xs_eval = range(interval[1], interval[2], length=N_steps)
    yvals = Matrix{Float64}(undef, N_steps, length(scaled_data))
    slopes = _slopes_buffer(scaled_data)
    for l in eachindex(scaled_data)
        d = scaled_data[l]
        _interpolate!(view(yvals, :, l), slopes, interp, d.xs, d.ys, xs_eval)
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
    for i in axes(yvals, 1)
        S += var(view(yvals, i, :))
    end
    S /= N_steps

    return S
end

function _S_weighted(scaled_data, interval, N_steps, weight_density, interp)

    dp_densities = zeros(length(scaled_data))
    for l in eachindex(scaled_data)
        i_l = searchsortedlast(scaled_data[l].xs, interval[1])
        i_r = searchsortedfirst(scaled_data[l].xs, interval[2])
        x_l = scaled_data[l].xs[i_l]
        x_r = scaled_data[l].xs[i_r]
        dp_densities[l] = (i_r - i_l - 1) / (x_r - x_l)
    end
    !weight_density && (dp_densities .= 1.0)

    xs_eval = range(interval[1], interval[2], length=N_steps)
    n = length(scaled_data)
    yvals = Matrix{Float64}(undef, n, N_steps)
    evals = Matrix{Float64}(undef, n, N_steps)  # uncertainties of the spline values
    slopes = _slopes_buffer(scaled_data)
    for l in eachindex(scaled_data)
        d = scaled_data[l]
        _interpolate!(view(yvals, l, :), slopes, interp, d.xs, d.ys, xs_eval)
        # linear interpolation of the error bars
        _linear_interpolation!(view(evals, l, :), d.xs, d.es, xs_eval)
    end
    # add small number to avoid division by zero
    evals .+= 1e-6

    # normalize y, s.t. minimum(y) == 0, maximum(y) == 1
    offset = minimum(yvals)
    yvals .-= offset
    norm_factor = maximum(yvals)
    yvals ./= norm_factor
    inv_norm = 1 / norm_factor

    # calculate quality of the scaling
    S = 0.0
    for i in axes(yvals, 2)
        y_sum = yvals[1, i]
        for l in 2:n
            y_sum += yvals[l, i]
        end
        y_mean = y_sum / n
        # The uncertainties were propagated with Measurements.jl before. The uncertainty of
        # the mean sums the squared uncertainties in Measurements.jl's order (n-1, n-3, ...
        # down to 1 or 2, then up to n in the other parity), so that the quality does not
        # change in the last digits.
        e2_sum = 0.0
        for l in n-1:-2:1
            e2_sum += abs2(evals[l, i] * inv_norm)
        end
        for l in (iseven(n) ? 2 : 1):2:n
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
