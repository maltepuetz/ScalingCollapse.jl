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
### small helper function to span the parameter space
# iterative function to generate all combinations of parameters
function _parameter_combinations!(p_space, current_combination=Float64[], all_combinations=Vector{Float64}[])
    # Base case: if p_space is empty, add the current combination to all_combinations
    if isempty(p_space)
        push!(all_combinations, current_combination)
        return
    end

    # Recursive case: iterate over the first list and combine its elements with the rest
    first_list = p_space[end]
    for value in first_list
        new_combination = [current_combination..., value] # Create a new combination
        _parameter_combinations!(p_space[1:end-1], new_combination, all_combinations) # Recursive call with the rest of p_space
    end

    return all_combinations
end

# Wrapper function to call the iterative function
function _parameter_combinations(p_space)
    all_combinations = _parameter_combinations!(p_space)
    return reshape(all_combinations .|> reverse, [length(p_space[i]) for i in eachindex(p_space)]...)
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
