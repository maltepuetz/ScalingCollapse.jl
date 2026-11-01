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
        S_p_space[i] = _nan_to_inf(sp.quality_scan(sp, ps; check_bounds=false))
    end
    all(isinf, S_p_space) && error(
        "The quality function is infinite for every parameter combination in p_space. " *
        "This usually means that the rescaled data of the different system sizes do not " *
        "overlap within the optimization interval dx (or contain too few points there). " *
        "Check p_space and dx, or provide starting_ps."
    )
    loc_starting_ps = local_minima(S_p_space, p_combos)
    if isempty(loc_starting_ps)
        # without a strict local minimum, the global minimum is tied with a neighbor,
        # i.e. the quality function is flat along some parameters
        i_min = argmin(S_p_space)
        flat_dims = [
            d for d in 1:ndims(S_p_space)
            if any(S_p_space[j] == S_p_space[i_min] && j[d] != i_min[d]
                   for j in get_neighbors(i_min, size(S_p_space)))
        ]
        p_names = hasproperty(sp, :sf) ? scaled_p_names(sp.sf) :
                  ["p$d" for d in 1:ndims(S_p_space)]
        error(
            "The parameter scan found no isolated minimum of the quality function: its " *
            "smallest value $(S_p_space[i_min]) at $(p_combos[i_min]) is shared by " *
            "neighboring grid points along $(join(p_names[flat_dims], ", ")). The data " *
            "does not seem to constrain these parameters. Consider fixing them in the " *
            "ScalingFunction, or adjust p_space or dx."
        )
    end
    sp.verbose && @info "Found $(length(loc_starting_ps)) local minima."
    sp.verbose && @info "Optimizing each starting point..."

    # optimize each starting point
    loc_minima = zeros(length(loc_starting_ps))
    loc_optimal_ps = Vector{Vector{Float64}}(undef, length(loc_starting_ps))
    f = ps -> _nan_to_inf(sp.quality_scan(sp, ps; check_bounds=true))
    for (i, l_ps) in enumerate(loc_starting_ps)
        result = optimize(f, l_ps, _bounded_optimizer(sp, f), sp.optimizer_options)
        loc_minima[i] = result.minimum
        loc_optimal_ps[i] = result.minimizer
    end

    # find global minimum
    indexmin = argmin(loc_minima)
    sp.starting_ps = loc_optimal_ps[indexmin]

    return nothing
end

function optimize_starting_ps!(sp)

    sp.verbose && @info "Starting optimization of provided starting parameters..."

    # optimize starting parameters
    f = ps -> _nan_to_inf(sp.quality(sp, ps; check_bounds=true))
    result = optimize(f, sp.starting_ps, _bounded_optimizer(sp, f), sp.optimizer_options)

    isfinite(result.minimum) || error(
        "The optimization did not find any parameters with a finite quality. Starting " *
        "from $(sp.starting_ps), the rescaled data of the different system sizes do not " *
        "overlap within the optimization interval dx (or contain too few points there). " *
        "Check starting_ps, p_space and dx."
    )
    if !Optim.converged(result) && sp.optimizer_options.show_warnings
        hint = if Optim.iteration_limit_reached(result)
            "Consider increasing the iteration limit via " *
            "optimizer_options=Optim.Options(iterations=...)."
        else
            "This can happen if the initial simplex around the starting point " *
            "$(sp.starting_ps) has vertices where the quality function is not finite (e.g. " *
            "no overlap of the rescaled data within dx). Consider widening p_space or " *
            "using a smaller initial simplex, e.g. " *
            "optimizer=Optim.NelderMead(initial_simplex=Optim.AffineSimplexer(b=0.01)) " *
            "or one step per parameter with StepSimplexer(steps)."
        end
        @warn "The Nelder-Mead optimization stopped after " *
              "$(Optim.iterations(result)) iterations without converging. The optimal " *
              "parameters might be inaccurate. " * hint
    end

    sp.minimum = result.minimum
    sp.optimal_ps = result.minimizer

    sp.verbose && @info "Found optimal parameters."
    return nothing
end



"""
    StepSimplexer(steps)

Initial Nelder-Mead simplex with one step per free parameter: vertex `i + 1` is the
starting point with `steps[i]` added to parameter `i`. Use it as
`optimizer=Optim.NelderMead(initial_simplex=StepSimplexer(steps))` when the free parameters
live on very different scales. `Optim.AffineSimplexer(a, b)` moves every parameter by
`b * p + a`, which follows the value of the parameter rather than its scale in the problem;
e.g. the default (`a = 0.025`, `b = 0.1`) moves a critical point `p_c ≈ 0.5` by 0.075, even
if the data only span a window of width 1e-3 around it.

`steps` must have one finite, nonzero entry per free parameter (in the order of
`scaled_p_names(sf)`; fixed parameters have no entry). A step of the order of the expected
error of the parameter is a good choice. Steps for which the quality is not finite are
halved, as for any initial simplex. In the profile error analysis (`error=:profile`), the
inner fits use the steps of the parameters they re-fit. This matters for narrow quality
valleys: if the default steps are much larger than the valley of a re-fitted parameter,
the inner fits can miss it, and the profile error comes out close to the conditional one.

# Examples
```julia
optimizer = ScalingCollapse.Optim.NelderMead(initial_simplex=StepSimplexer([1e-3, 0.05]))
ScalingProblem(xs, ys, Ls; sf=ScalingFunction(:x), p_space=[0.48:0.002:0.52, 0.8:0.1:2.0],
    optimizer=optimizer)
```
"""
struct StepSimplexer <: Optim.Simplexer
    steps::Vector{Float64}
    function StepSimplexer(steps::AbstractVector{<:Real})
        all(s -> isfinite(s) && !iszero(s), steps) || throw(ArgumentError(
            "The steps of a StepSimplexer must be finite and nonzero, got $(steps)."
        ))
        return new(Float64.(steps))
    end
end

function Optim.simplexer(S::StepSimplexer, initial_x::Tx) where {Tx}
    length(S.steps) == length(initial_x) || throw(DimensionMismatch(
        "The StepSimplexer has $(length(S.steps)) steps, but the starting point has " *
        "$(length(initial_x)) parameters."
    ))
    simplex = Tx[copy(initial_x) for _ in 1:length(initial_x)+1]
    for (k, j) in enumerate(eachindex(initial_x))
        simplex[k+1][j] += S.steps[k]
    end
    return simplex
end

# the initial simplex for a fit of the free parameters `indices` only
_restrict_simplexer(S::Optim.Simplexer, indices) = S
_restrict_simplexer(S::StepSimplexer, indices) = StepSimplexer(S.steps[indices])



##### Little Helper Functions #####
### keep the initial Nelder-Mead simplex at finite quality values
# The quality functions return Inf outside p_space (check_bounds=true), but also inside p_space
# when the rescaled data do not overlap within dx. If any vertex of the initial simplex is Inf,
# Optim's Nelder-Mead (>= 1.12) reports a non-finite convergence measure and stops at iteration
# 0 without moving. This wrapper moves such vertices back into p_space and then towards the
# starting point until the quality function f is finite there.
const _MAX_SIMPLEX_HALVINGS = 20

struct BoundedSimplexer{S<:Optim.Simplexer,F} <: Optim.Simplexer
    simplexer::S
    lower::Vector{Float64}
    upper::Vector{Float64}
    f::F
end

function Optim.simplexer(S::BoundedSimplexer, initial_x)
    simplex = Optim.simplexer(S.simplexer, initial_x)
    for vertex in simplex, i in eachindex(vertex)
        lo, hi = S.lower[i], S.upper[i]
        lo <= vertex[i] <= hi && continue
        # reflect the step through the starting point, clamp if it still leaves p_space
        vertex[i] = clamp(2 * initial_x[i] - vertex[i], lo, hi)
    end
    for vertex in simplex
        vertex == initial_x && continue
        # halve the step from the starting point until the quality is finite
        n = 0
        while !isfinite(S.f(vertex)) && n < _MAX_SIMPLEX_HALVINGS
            vertex .= initial_x .+ (vertex .- initial_x) ./ 2
            n += 1
        end
    end
    return simplex
end

function _bounded_optimizer(sp, f)
    lower = Float64[first(p) for p in sp.p_space]
    upper = Float64[last(p) for p in sp.p_space]
    return Optim.NelderMead(
        BoundedSimplexer(sp.optimizer.initial_simplex, lower, upper, f),
        sp.optimizer.parameters,
    )
end

# treat NaN values of the quality function as invalid parameters
_nan_to_inf(S) = isnan(S) ? oftype(S, Inf) : S

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
