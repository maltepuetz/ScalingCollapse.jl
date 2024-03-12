function optimize_parameters!(sp)

    # scan parameter space
    if !sp.skip_scan
        parameter_scan!(sp)
    end

    # optimize parameters
    optimize_starting_ps!(sp)

    return nothing
end


function parameter_scan!(sp)

    sp.verbose && @info "Scanning parameter space..."
    p_combos = _parameter_combinations(sp.p_space)
    S_p_space = zeros(size(p_combos))
    for (i, ps) in enumerate(p_combos)
        S_p_space[i] = sp.quality_scan(sp, ps; check_bounds=false)
    end
    loc_starting_ps = local_minima(S_p_space, p_combos)
    sp.verbose && @info "Found $(length(loc_starting_ps)) local minima."
    sp.verbose && @info "Optimizing each starting point..."

    # optimize each starting point
    loc_minima = zeros(length(loc_starting_ps))
    loc_optimal_ps = Vector{Vector{Float64}}(undef, length(loc_starting_ps))
    for (i, l_ps) in enumerate(loc_starting_ps)
        result = optimize(
            ps -> sp.quality_scan(sp, ps; check_bounds=true),
            l_ps,
            NelderMead(; initial_simplex=Optim.AffineSimplexer(; b=0.1))
        )
        loc_minima[i] = result.minimum
        loc_optimal_ps[i] = result.minimizer
    end

    # find global minimum
    indexmin = argmin(loc_minima)
    sp.starting_ps = loc_optimal_ps[indexmin]

    return nothing
end

function optimize_starting_ps!(sp)

    sp.verbose && @info "Starting optimzation of provided startign parameters..."

    # optimize starting parameters
    result = optimize(
        ps -> sp.quality(sp, ps; check_bounds=true),
        sp.starting_ps,
        NelderMead(; initial_simplex=Optim.AffineSimplexer(; b=0.1))
    )

    sp.minimum = result.minimum
    sp.optimal_ps = result.minimizer

    sp.verbose && @info "Found optimal parameters."
    return nothing
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
    elseif length(X) == 4
        return [
            [X[1][i], X[2][j], X[3][k], X[4][l]]
            for i in eachindex(X[1]), j in eachindex(X[2]), k in eachindex(X[3]), l in eachindex(X[4])
        ]
    elseif length(X) == 5
        return [
            [X[1][i], X[2][j], X[3][k], X[4][l], X[5][m]]
            for i in eachindex(X[1]), j in eachindex(X[2]), k in eachindex(X[3]), l in eachindex(X[4]), m in eachindex(X[5])
        ]
    else
    error(
        "Not implemented for more than 5 parameters - " *
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
