function quality_hybrid(sp, parameters; check_bounds=false, scanmode=false)

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
    N_steps = 100
    S = 0.0
    if sp.errors_defined && !scanmode
        S = _S_weighted(scaled_data, interval, N_steps)
    else
        S = _S_unweighted(scaled_data, interval, N_steps)
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

function _S_weighted(scaled_data, interval, N_steps)

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
    #S *= (interval[2] - interval[1])

    return S
end





#=
function quality_hybrid(sp, parameters; check_bounds=false)

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

    N_steps = 100
    S = 0.0
    sum_dp_densities = 0.0

    for l in eachindex(sp.data)

        if sp.errors_defined
            S_l, dp_density_l = _S_weighted(scaled_data, l, interval, N_steps)
            S += S_l * dp_density_l
            sum_dp_densities += dp_density_l
        else
            S_l, dp_density_l = _S_unweighted(scaled_data, l, interval, N_steps)
            S += S_l * dp_density_l
            sum_dp_densities += dp_density_l

            #@info "l = $l, S_l = $S_l, dp_density_l = $dp_density_l"
        end
    end
    S /= sum_dp_densities
    S /= N_steps

    return S
end

function _S_weighted(scaled_data, l, interval, N_steps)

    res = 0.0
    x_vals = range(interval[1], interval[2], length=N_steps)
    y_spline = Spline1D(scaled_data[l].xs, scaled_data[l].ys, k=3)  #! maybe k=1?
    e_spline = Spline1D(scaled_data[l].xs, scaled_data[l].es, k=1)

    # set the data point density in the interval (used for weighting)
    i_l = searchsortedlast(scaled_data[l].xs, interval[1])
    i_r = searchsortedfirst(scaled_data[l].xs, interval[2])
    dp_density = (i_r - i_l - 1) / (scaled_data[l].xs[i_r] - scaled_data[l].xs[i_l])
    dp_density == 0 && (dp_density = 1 / (scaled_data[l].xs[i_r] - scaled_data[l].xs[i_l]))

    # set structure to store the left and right element of the other data sets
    N_data_sets = length(scaled_data)
    ints = _IntervalsW(N_data_sets - 1)

    # find the elemenets from the other data sets that are closest xvals[1]
    for (i, j) in enumerate(setdiff(eachindex(scaled_data), l))

        # find indices of left element from interval[1]
        ind = searchsortedlast(scaled_data[j].xs, interval[1])
        ind == 0 && error("No element found in interval.")
        ind == length(scaled_data[j].xs) && error("No element found in interval.")

        ints.indices[i] = ind
        ints.xs_l[i] = scaled_data[j].xs[ind]
        ints.ys_l[i] = scaled_data[j].ys[ind]
        ints.es_l[i] = scaled_data[j].es[ind]
        ints.xs_r[i] = scaled_data[j].xs[ind+1]
        ints.ys_r[i] = scaled_data[j].ys[ind+1]
        ints.es_r[i] = scaled_data[j].es[ind+1]
    end

    # initialize fit
    fit = _Fit(ints)

    for x in x_vals

        # update right (and left) elements if necessary
        changes = update!(ints, scaled_data, l, x)

        if changes
            update!(fit, ints)
        end

        Y, dY2 = fit(x)
        res += (y_spline(x) - Y)^2 / (e_spline(x)^2 + dY2)
    end
    return res, dp_density
end

function _S_unweighted(scaled_data, l, interval, N_steps)

    res = 0.0
    x_vals = range(interval[1], interval[2], length=N_steps)
    y_spline = Spline1D(scaled_data[l].xs, scaled_data[l].ys, k=3)  #! maybe k=1?

    # set the data point density in the interval (used for weighting)
    i_l = searchsortedlast(scaled_data[l].xs, interval[1])
    i_r = searchsortedfirst(scaled_data[l].xs, interval[2])
    dp_density = (i_r - i_l - 1) / (scaled_data[l].xs[i_r] - scaled_data[l].xs[i_l])
    dp_density == 0 && (dp_density = 1 / (scaled_data[l].xs[i_r] - scaled_data[l].xs[i_l]))

    # set structure to store the left and right element of the other data sets
    N_data_sets = length(scaled_data)
    ints = _IntervalsUW(N_data_sets - 1)

    # find the elemenets from the other data sets that are closest xvals[1]
    for (i, j) in enumerate(setdiff(eachindex(scaled_data), l))

        # find indices of left element from interval[1]
        ind = searchsortedlast(scaled_data[j].xs, interval[1])
        ind == 0 && error("No element found in interval.")
        ind == length(scaled_data[j].xs) && error("No element found in interval.")

        ints.indices[i] = ind
        ints.xs_l[i] = scaled_data[j].xs[ind]
        ints.ys_l[i] = scaled_data[j].ys[ind]
        ints.xs_r[i] = scaled_data[j].xs[ind+1]
        ints.ys_r[i] = scaled_data[j].ys[ind+1]
    end

    # initialize fit
    fit = _Fit(ints)

    for x in x_vals

        # update right (and left) elements if necessary
        changes = update!(ints, scaled_data, l, x)

        if changes
            update!(fit, ints)
        end

        Y, dY2 = fit(x)
        res += (y_spline(x) - Y)^2 / dY2
    end
    return res, dp_density
end

mutable struct _IntervalsW
    indices::Vector{Int}
    xs_l::Vector{Float64}
    xs_r::Vector{Float64}
    ys_l::Vector{Float64}
    ys_r::Vector{Float64}
    es_l::Vector{Float64}
    es_r::Vector{Float64}

    function _IntervalsW(N)
        new(
            Vector{Int}(undef, N),
            Vector{Float64}(undef, N),
            Vector{Float64}(undef, N),
            Vector{Float64}(undef, N),
            Vector{Float64}(undef, N),
            Vector{Float64}(undef, N),
            Vector{Float64}(undef, N),
        )
    end
end

mutable struct _IntervalsUW
    indices::Vector{Int}
    xs_l::Vector{Float64}
    xs_r::Vector{Float64}
    ys_l::Vector{Float64}
    ys_r::Vector{Float64}

    function _IntervalsUW(N)
        new(
            Vector{Int}(undef, N),
            Vector{Float64}(undef, N),
            Vector{Float64}(undef, N),
            Vector{Float64}(undef, N),
            Vector{Float64}(undef, N),
        )
    end
end

# this function updates the right (and left) elements if x > xs_r[i]
function update!(ints::_IntervalsW, scaled_data, l, x)

    changes = false
    for (i, j) in enumerate(setdiff(eachindex(scaled_data), l))
        if x > ints.xs_r[i]
            ints.indices[i] += 1
            ints.xs_l[i] = ints.xs_r[i]
            ints.ys_l[i] = ints.ys_r[i]
            ints.es_l[i] = ints.es_r[i]
            ints.xs_r[i] = scaled_data[j].xs[ints.indices[i]+1]
            ints.ys_r[i] = scaled_data[j].ys[ints.indices[i]+1]
            ints.es_r[i] = scaled_data[j].es[ints.indices[i]+1]
            changes = true
        end
    end

    return changes
end

function update!(ints::_IntervalsUW, scaled_data, l, x)

    changes = false
    for (i, j) in enumerate(setdiff(eachindex(scaled_data), l))
        if x > ints.xs_r[i]
            ints.indices[i] += 1
            ints.xs_l[i] = ints.xs_r[i]
            ints.ys_l[i] = ints.ys_r[i]
            ints.xs_r[i] = scaled_data[j].xs[ints.indices[i]+1]
            ints.ys_r[i] = scaled_data[j].ys[ints.indices[i]+1]
            changes = true
        end
    end

    return changes
end

struct _Fit
    K::Float64
    Kx::Float64
    Ky::Float64
    Kxx::Float64
    Kxy::Float64
    Delta::Float64
    a::Float64
    b::Float64

    function _Fit(ints::_IntervalsW)

        ws_l = 1 ./ ints.es_l .^ 2
        ws_r = 1 ./ ints.es_r .^ 2
        max_noninf = max(maximum(ws_l[isfinite.(ws_l)]), maximum(ws_r[isfinite.(ws_r)]))

        # if there are single zero errors, set weigth to max(ws \ Inf)
        for i in eachindex(ws_l)
            isinf(ws_l[i]) && (ws_l[i] = max_noninf)
            ising(ws_r[i]) && (ws_r[i] = max_noninf)
        end

        K = sum(ws_l) + sum(ws_r)
        Kx = sum(ws_l .* ints.xs_l) + sum(ws_r .* ints.xs_r)
        Ky = sum(ws_l .* ints.ys_l) + sum(ws_r .* ints.ys_r)
        Kxx = sum(ws_l .* ints.xs_l .^ 2) + sum(ws_r .* ints.xs_r .^ 2)
        Kxy = sum(ws_l .* ints.xs_l .* ints.ys_l) + sum(ws_r .* ints.xs_r .* ints.ys_r)
        Delta = K * Kxx - Kx^2
        a = (Kxx * Ky - Kx * Kxy) / Delta
        b = (K * Kxy - Kx * Ky) / Delta




        new(K, Kx, Ky, Kxx, Kxy, Delta, a, b)
    end

    function _Fit(ints::_IntervalsUW)

        # for the unweighted fit we can just set all weights to 1 / N
        N = length(ints.indices) * 2
        K = 1.0
        Kx = (sum(ints.xs_l) + sum(ints.xs_r)) / N
        Ky = (sum(ints.ys_l) + sum(ints.ys_r)) / N
        Kxx = (sum(ints.xs_l .^ 2) + sum(ints.xs_r .^ 2)) / N
        Kxy = (sum(ints.xs_l .* ints.ys_l) + sum(ints.xs_r .* ints.ys_r)) / N
        Delta = K * Kxx - Kx^2
        a = (Kxx * Ky - Kx * Kxy) / Delta
        b = (K * Kxy - Kx * Ky) / Delta

        if isnan(a) || !isfinite(a)
            @warn "a = ", a
        elseif isnan(b) || !isfinite(b)
            @warn "b = ", b
        elseif isnan(Delta) || !isfinite(Delta)
            @warn "Delta = ", Delta
        elseif isnan(K) || !isfinite(K)
            @warn "K = ", K
        elseif isnan(Kx) || !isfinite(Kx)
            @warn "Kx = ", Kx
        elseif isnan(Ky) || !isfinite(Ky)
            @warn "Ky = ", Ky
        elseif isnan(Kxx) || !isfinite(Kxx)
            @warn "Kxx = ", Kxx
        elseif isnan(Kxy) || !isfinite(Kxy)
            @warn "Kxy = ", Kxy
        end




        new(K, Kx, Ky, Kxx, Kxy, Delta, a, b)
    end
end

function update!(fit::_Fit, ints)
    fit = _Fit(ints)
    return nothing
end

function (f::_Fit)(x)
    Y = f.a + x * f.b
    dY2 = 1 / f.Delta * (f.Kxx - 2 * x * f.Kx + x^2 * f.K)
    return Y, dY2
end

=#
