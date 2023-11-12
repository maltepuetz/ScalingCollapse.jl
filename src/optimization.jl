

function optimize_parameters(
    data::Vector{Data},
    sf::ScalingFunction,
    p_space,
    verbose;
    kwargs...
)

    # check if starting parameters have been provided
    starting_ps = get(kwargs, :starting_ps, nothing)
    if starting_ps !== nothing
        return optimize_parameters(
            data,
            sf,
            p_space,
            starting_ps,
            verbose;
            kwargs...
        )
    end

    # do a scan of the parameter space to find good starting points
    verbose && @info "Scanning parameter space..."
    p_combos = _parameter_combinations(p_space)
    ssr_p_space = zeros(size(p_combos))
    for (i, ps) in enumerate(p_combos)
        ssr_p_space[i] = squared_sum_residuals(
            data, sf.f, ps; kwargs...
        )
    end
    loc_starting_ps = local_minima(ssr_p_space, p_combos)
    verbose && @info "Found $(length(loc_starting_ps)) local minima."
    verbose && @info "Optimizing each starting point..."

    # optimize each starting point
    loc_minima = zeros(length(loc_starting_ps))
    loc_optimal_ps = Vector{Vector{Float64}}(undef, length(loc_starting_ps))
    for (i, l_ps) in enumerate(loc_starting_ps)
        result = optimize(
            ps -> squared_sum_residuals(
                data, sf.f, ps;
                check_bounds=true,
                p_space=p_space,
                kwargs...),
            l_ps
        )
        loc_minima[i] = result.minimum
        loc_optimal_ps[i] = result.minimizer
    end

    # find global minimum
    indexmin = argmin(loc_minima)
    minimum = loc_minima[indexmin]
    optimal_ps = loc_optimal_ps[indexmin]

    return optimal_ps, minimum
end

function optimize_parameters(
    data::Vector{Data},
    sf::ScalingFunction,
    p_space,
    starting_ps,
    verbose;
    kwargs...
)
    verbose && @info "Starting optimzation of provided startign parameters..."
    # optimize starting parameters
    result = optimize(
        ps -> squared_sum_residuals(
            data, sf.f, ps;
            check_bounds=true,
            p_space=p_space,
            kwargs...),
        starting_ps
    )

    minimum = result.minimum
    optimal_ps = result.minimizer

    return optimal_ps, minimum
end

















"""
    squared_sum_residuals(data::Vector{Data}, scaling_function, parameters; kwargs...)

Calculate the squared sum of residuals for a given scaling function and parameters.

# Arguments
- `data::Vector{Data}`: data to be scaled
- `scaling_function::Function`: scaling function
- `parameters`: parameters for the scaling function
"""
function squared_sum_residuals(data::Vector{Data}, scaling_function, parameters; kwargs...)

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

    data_scaled = [scaling_function(data[i], parameters...) for i in eachindex(data)]

    # set interval in which we want to optimize
    dx = get(kwargs, :dx, [-Inf, Inf])
    interval = Vector{Float64}(undef, 2)
    interval[1] = max(
        maximum(data_scaled[i].xs[1] for i in eachindex(data_scaled)),
        dx[1]
    )
    interval[2] = min(
        minimum(data_scaled[i].xs[end] for i in eachindex(data_scaled)),
        dx[2]
    )

    # check whether the constructed optimization interval is valid
    interval[1] > interval[2] && return Inf

    # create spline for each system size
    splines = [
        Spline1D(data_scaled[i].xs, data_scaled[i].ys, k=3) for i in eachindex(data_scaled)
    ]

    N_steps = 100
    xvals = range(interval[1], interval[2], length=N_steps)
    yvals = zeros(N_steps, length(data_scaled))
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

    # calculate squared sum of residuals
    ssr = 0.0
    for i in axes(yvals, 1)
        ssr += var(yvals[i, :])
    end
    ssr /= N_steps

    return ssr
end






##### Little Helper Functions #####
# small helper function to span the parameter space
function _parameter_combinations(p_space)
    X = p_space  # for convenience
    if length(X) == 1
        return [[X[1][i]] for i in eachindex(X[1])]
    elseif length(X) == 2
        return [
            [X[1][i], X[2][j]] for i in eachindex(X[1]), j in eachindex(X[2])
        ]
    elseif length(X) == 3
        return [
            [X[1][i], X[2][j], X[3][k]]
            for i in eachindex(X[1]), j in eachindex(X[2]), k in eachindex(X[3])
        ]
    end
    error(
        "Not implemented for more than 3 parameters - " *
        "feel free to create an issue on GitHub! :)"
    )
end

# small helper function to find local minimal parameter combinations
function local_minima(arr, p_combos)
    indices = CartesianIndices(arr)
    dims = size(arr)
    local_minima = Vector{Vector{Float64}}(undef, 0)
    for i in indices
        neighbors = get_neighbors(i, dims)
        if all(j -> arr[i] < arr[j], neighbors)
            push!(local_minima, p_combos[i])
        end
    end

    return local_minima
end

# small helper function to get the neighbors of an index in an array
function get_neighbors(i::CartesianIndex, dims)
    neighbors = Vector{CartesianIndex{length(dims)}}(undef, 0)

    d = 1
    inds = Tuple(i)
    while d <= length(dims)
        i[d] > 1 && push!(
            neighbors,
            CartesianIndex(inds[1:d-1]..., i[d] - 1, inds[d+1:end]...)
        )
        i[d] < dims[d] && push!(
            neighbors,
            CartesianIndex(inds[1:d-1]..., i[d] + 1, inds[d+1:end]...)
        )
        d += 1
    end

    return neighbors
end
