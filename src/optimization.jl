

function optimize_parameters(
    data::Vector{Data},
    sf::ScalingFunction,
    p_space,
    verbose;
    kwargs...
)

    # check if algorithm is specified
    algorithm = get(kwargs, :algorithm, :hybrid)

    # check if starting parameters have been provided
    starting_ps = get(kwargs, :starting_ps, nothing)
    if starting_ps !== nothing

        # set correct quality function
        quality = quality_houdayer
        algorithm == :spline && (quality = quality_spline)

        return optimize_parameters(
            data,
            sf,
            quality,
            p_space,
            starting_ps,
            verbose;
            kwargs...
        )
    end

    # set correct quality function
    quality = quality_spline
    algorithm == :houdayer && (quality = quality_houdayer)

    # do a scan of the parameter space to find good starting points
    verbose && @info "Scanning parameter space..."
    p_combos = _parameter_combinations(p_space)
    S_p_space = zeros(size(p_combos))
    for (i, ps) in enumerate(p_combos)
        S_p_space[i] = quality(
            data, sf.f, ps; kwargs...
        )
    end
    loc_starting_ps = local_minima(S_p_space, p_combos)
    verbose && @info "Found $(length(loc_starting_ps)) local minima."
    verbose && @info "Optimizing each starting point..."

    # optimize each starting point
    loc_minima = zeros(length(loc_starting_ps))
    loc_optimal_ps = Vector{Vector{Float64}}(undef, length(loc_starting_ps))
    for (i, l_ps) in enumerate(loc_starting_ps)
        result = optimize(
            ps -> quality(
                data, sf.f, ps;
                check_bounds=true,
                p_space=p_space,
                kwargs...),
            l_ps,
            NelderMead(; initial_simplex=Optim.AffineSimplexer(; b=0.1))
        )
        loc_minima[i] = result.minimum
        loc_optimal_ps[i] = result.minimizer
    end

    # find global minimum
    indexmin = argmin(loc_minima)
    minimum = loc_minima[indexmin]
    optimal_ps = loc_optimal_ps[indexmin]

    if algorithm == :hybrid
        # optimize from the global minimum of the quality_spline function
        # using the quality_houdayer function

        verbose && @info "Optimizing global minimum with Houdayer method..."
        result = optimize(
            ps -> quality_houdayer(
                data, sf.f, ps;
                check_bounds=true,
                p_space=p_space,
                kwargs...),
            optimal_ps,
            NelderMead(; initial_simplex=Optim.AffineSimplexer(; b=0.1))
        )
        minimum = result.minimum
        optimal_ps = result.minimizer
    end

    return optimal_ps, minimum
end

function optimize_parameters(
    data::Vector{Data},
    sf::ScalingFunction,
    quality,
    p_space,
    starting_ps,
    verbose;
    kwargs...
)
    verbose && @info "Starting optimzation of provided startign parameters..."
    # optimize starting parameters
    result = optimize(
        ps -> quality(
            data, sf.f, ps;
            check_bounds=true,
            p_space=p_space,
            kwargs...),
        starting_ps,
        NelderMead(; initial_simplex=Optim.AffineSimplexer(; b=0.1))
    )

    minimum = result.minimum
    optimal_ps = result.minimizer

    verbose && @info "Found optimal parameters."
    return optimal_ps, minimum
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
