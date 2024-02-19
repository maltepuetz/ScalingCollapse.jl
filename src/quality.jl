abstract type QualityFunction end

function Base.show(io::IO, qf::QualityFunction)
    print(io, typeof(qf), " <: QualityFunction")
end
"""
    SingleMasterCurve()

Quality function that interpolates all scaled data points to a single master curve. For
every data point, we take the closest data points, to the left and to the right to
interpolate a linear function. The quality is then calculated as the sum of the squared
differences between the data points and the interpolated values.
"""
struct SingleMasterCurve <: QualityFunction end


"""
    Spline(; N_steps::Int=100, weight_density::Bool=true)

Quality function that uses a spline interpolation to calculate the quality of a scaling.

# How does it work?
TO BE DOCUMENTED

# Keyword arguments
- `N_steps::Int=100`: Number of points used in the interval to calculate the quality.
- `weight_density::Bool=true`: If `true`, the quality is weighted by the density of data
  points in the interval.
"""
struct Spline <: QualityFunction
    N_steps::Int
    weight_density::Bool
    scan_mode::Bool
    function Spline(; N_steps::Int=100, weight_density=true, scan_mode=false)
        new(N_steps, weight_density, scan_mode)
    end
end

"""
    Houdayer()

Quality function that uses the Houdayer & Hartmann method to calculate the quality of a
  scaling. It is less stable than the spline method.

# How does it work?
TO BE DOCUMENTED
"""
struct Houdayer <: QualityFunction end



##### Single Master Curve Quality Function #####
function (::SingleMasterCurve)(sp, parameters; check_bounds=false)

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
            ybar = ((xs[i+1] - xs[i]) * ys[i-1] - (xs[i-1] - xs[i]) * ys[i+1]) / (xs[i+1] - xs[i-1])
            expected_var = es[i]^2 + ((xs[i+1] - xs[i]) / (xs[i+1] - xs[i-1]))^2 * es[i-1]^2 + ((xs[i-1] - xs[i]) / (xs[i+1] - xs[i-1]))^2 * es[i+1]^2
            S += (ys[i] - ybar)^2 / expected_var
        end
    else
        for i in eachindex(xs)[2:end-1]
            # interpolate linearly
            ybar = ((xs[i+1] - xs[i]) * ys[i-1] - (xs[i-1] - xs[i]) * ys[i+1]) / (xs[i+1] - xs[i-1])
            S += (ys[i] - ybar)^2
        end
    end
    S /= (length(xs) - 2)
    return S
end



##### Spline Quality Function #####
function (sqf::Spline)(sp, parameters; check_bounds=false)

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
    yvals = zeros(N_steps, length(scaled_data))
    for l in eachindex(scaled_data)
        for (i, x) in enumerate(xvals)
            yvals[i, l] = y_splines[l](x)
        end
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
        S += var(yvals[i, :])
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
    yvals = zeros(Measurement, length(scaled_data), N_steps)

    for (i, x) in enumerate(xvals)
        for l in eachindex(scaled_data)
            y = y_splines[l](x)
            e = e_splines[l](x) + 1e-6  # add small number to avoid division by zero
            yvals[l, i] = measurement(y, e)
        end
    end

    # normalize y, s.t. minimum(y) == 0, maximum(y) == 1
    offset = minimum(yvals).val
    for i in eachindex(yvals)
        yvals[i] -= offset
    end
    norm_factor = maximum(yvals).val
    for i in eachindex(yvals)
        yvals[i] /= norm_factor
    end

    # calculate quality of the scaling
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


##### Houdayer Quality Function #####
function (::Houdayer)(sp, parameters; check_bounds=false)

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

    N = 0  # counts number of terms in the sum
    S = 0.0  # quality of the scaling

    for (i, d) in enumerate(scaled_data)

        j = searchsortedfirst(d.xs, interval[1])  # index of first element in interval

        while (j <= length(d.xs)) && (d.xs[j] <= interval[2])

            Y_ij, dY2_ij = fit_mastercurve(scaled_data, i, j)
            Y_ij == Inf && (j += 1; continue)  # no contribution

            S += (d.ys[j] - Y_ij)^2 / (d.es[j]^2 + dY2_ij)

            j += 1
            N += 1
        end
    end
    S /= N

    return S
end

function fit_mastercurve(scaled_data, i, j)

    # set arrays to store set of data used for fitting
    xs = Float64[]
    ys = Float64[]
    ws = Float64[]

    for (i_prime, d_prime) in enumerate(scaled_data)
        i_prime == i && continue
        j_prime = searchsortedlast(d_prime.xs, scaled_data[i].xs[j])
        ((j_prime == 0) || (j_prime == length(d_prime.xs))) && continue  # no contribution
        push!(xs, d_prime.xs[j_prime])
        push!(xs, d_prime.xs[j_prime+1])
        push!(ys, d_prime.ys[j_prime])
        push!(ys, d_prime.ys[j_prime+1])
        push!(ws, d_prime.es[j_prime])
        push!(ws, d_prime.es[j_prime+1])
    end

    length(xs) == 0 && return Inf, Inf  # no contribution

    # check whether all elements in ws are zeros
    all(ws .== 0.0) && (ws .= 1.0)
    ws .= 1 ./ ws .^ 2

    # if there are single zero errors, set weight to max(ws \ Inf)
    for i in eachindex(ws)
        isinf(ws[i]) && (ws[i] = maximum(ws[isfinite.(ws)]))
    end

    K = sum(ws)
    Kx = sum(ws .* xs)
    Ky = sum(ws .* ys)
    Kxx = sum(ws .* xs .^ 2)
    Kxy = sum(ws .* xs .* ys)
    Delta = K * Kxx - Kx^2

    a = (Kxx * Ky - Kx * Kxy) / Delta
    b = (K * Kxy - Kx * Ky) / Delta

    Y_ij = a + scaled_data[i].xs[j] * b
    dY2_ij = 1 / Delta * (Kxx - 2 * scaled_data[i].xs[j] * Kx + scaled_data[i].xs[j]^2 * K)

    return Y_ij, dY2_ij
end
