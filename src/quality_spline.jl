function quality_spline(data::Vector{Data}, scaling_function, parameters; kwargs...)

    # check if parameters are in bounds -- if not return Inf
    if get(kwargs, :check_bounds, false)
        p_space = get(
            kwargs, :p_space, [0.1:0.1:3.0 for _ in 1:3]
        )
        for (i, p) in enumerate(parameters)
            if p < p_space[i][1] || p > p_space[i][end]
                return Inf
            end
        end
    end

    scaled_data = [scaling_function(data[i], parameters...) for i in eachindex(data)]

    # set interval in which we want to optimize
    dx = get(kwargs, :dx, [-Inf, Inf])
    interval = Vector{Float64}(undef, 2)
    interval[1] = max(
        maximum(scaled_data[i].xs[1] for i in eachindex(scaled_data)),
        dx[1]
    )
    interval[2] = min(
        minimum(scaled_data[i].xs[end] for i in eachindex(scaled_data)),
        dx[2]
    )

    # check whether the constructed optimization interval is valid
    interval[1] > interval[2] && return Inf

    # create spline for each system size
    splines = [
        Spline1D(scaled_data[i].xs, scaled_data[i].ys, k=3) for i in eachindex(scaled_data)
    ]

    N_steps = 100
    xvals = range(interval[1], interval[2], length=N_steps)
    yvals = zeros(N_steps, length(scaled_data))
    for (l, spline) in enumerate(splines)
        for (i, x) in enumerate(xvals)
            yvals[i, l] = spline(x)
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
