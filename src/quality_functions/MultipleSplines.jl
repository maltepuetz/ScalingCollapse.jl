"""
    MultipleSplines(; N_steps::Int=100, weight_density::Bool=true)

Quality function that uses a spline interpolation to calculate the quality of a scaling.

# How does it work?
TO BE DOCUMENTED

# Keyword arguments
- `N_steps::Int=100`: Number of points used in the interval to calculate the quality.
- `weight_density::Bool=true`: If `true`, the quality is weighted by the density of data
  points in the interval.
"""
struct MultipleSplines <: QualityFunction
    N_steps::Int
    weight_density::Bool
    scan_mode::Bool
    function MultipleSplines(; N_steps::Int=100, weight_density=true, scan_mode=false)
        new(N_steps, weight_density, scan_mode)
    end
end


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

    # check whether the constructed optimization interval is valid
    interval[1] > interval[2] && return Inf


    # calculate S
    S = 0.0
    if sp.errors_defined && !sqf.scan_mode
        S = _S_weighted(scaled_data, interval, sqf.N_steps, sqf.weight_density)
    else
        S = _S_unweighted(scaled_data, interval, sqf.N_steps)
    end

    return S
end

function _S_unweighted(scaled_data, interval, N_steps)

    # create spline for each system size
    y_splines = [
        Spline1D(scaled_data[i].xs, scaled_data[i].ys, k=3) for i in eachindex(scaled_data)
    ]

    xvals = range(interval[1], interval[2], length=N_steps)
    xs_eval = collect(xvals)
    yvals = zeros(N_steps, length(scaled_data))
    for l in eachindex(scaled_data)
        yvals[:, l] = evaluate(y_splines[l], xs_eval)
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

function _S_weighted(scaled_data, interval, N_steps, weight_density)

    # create spline for each system size
    y_splines = [
        Spline1D(scaled_data[i].xs, scaled_data[i].ys, k=3) for i in eachindex(scaled_data)
    ]
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
        yvals[l, :] = evaluate(y_splines[l], xs_eval)
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
