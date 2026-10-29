"""
    residuals(sp; kwargs...)

Calculate the residuals of the `ScalingProblem` around the optimal parameters or for a
given parameter space.

# Arguments
- `sp::ScalingProblem`: The `ScalingProblem` to calculate the residuals for.

# Keyword Arguments
- `p_space::Vector{Vector{Float64}}`: The parameter space to calculate the residuals for.
    If not provided, the parameter space is calculated from the optimal parameters,
    their errors and the number of steps.
- `N_steps::Int=50`: The number of steps to use for the default parameter space.
- `dims::Vector{Int}=[1, ...]`: The dimensions to calculate the residuals for. If not
    specified, all dimensions are used.

# Returns
- `p_space::Vector{StepRangeLen}`: The parameter space used for the calculation.
- `residuals::Array{Float64}`: The residuals for the given parameter space.

# Example
```julia
# lets say our xs, ys and Ls are data for the susceptibility of the Ising model
using ScalingCollapse
sp = ScalingProblem(xs, ys, Ls;
    sf=ScalingFunction(:xny; p_names=["T_c", "nu", "gamma"]),
    dx=[-1.0, 2.0],
)

# now we can calculate the residual landscape around the optimal parameters as follows:
p_space, residuals = residuals(sp)
# in the above case length(p_space) == 3 and size(residuals) == (50, 50, 50)

# we can specify the dimensions we are interested in as follows:
p_space, residuals = residuals(sp; dims=[1, 2])
# now length(p_space) == 2 and size(residuals) == (50, 50)
# we fixed the third parameter (in this case gamma) to its optimal value and calculated
# a 2D cut through the 3D residual landscape
```
"""
function residuals(sp; kwargs...)  # p_space, N_steps=50, dims

    # get active dimensions
    dims = get(kwargs, :dims, [d for d in 1:n_parameters(sp.sf)]) |> sort
    N = length(dims)
    N == 0 && throw(ArgumentError("No active dimensions provided."))
    maximum(dims) > n_parameters(sp.sf) && throw(ArgumentError("Dimension out of bounds."))
    length(unique(dims)) != N && throw(ArgumentError("Duplicate dimensions provided."))

    # set parameter space
    N_steps = get(kwargs, :N_steps, 51)
    p_space = get(kwargs, :p_space) do
        # default: optimal parameters ± 3 errors, which requires valid errors
        p_names = scaled_p_names(sp.sf)
        for i in dims
            e = sp.optimal_ps_error[i]
            (isfinite(e) && e > 0.0) || throw(ArgumentError(
                "Cannot construct the default p_space around the optimal parameters: " *
                "the error of $(p_names[i]) is $e. Pass p_space explicitly, or solve " *
                "the ScalingProblem with error=true."
            ))
        end
        [
            i in dims ? range(
                sp.optimal_ps[i] - 3 * sp.optimal_ps_error[i],
                sp.optimal_ps[i] + 3 * sp.optimal_ps_error[i],
                length=N_steps
            ) : _single_point_range(sp.optimal_ps[i])
            for i in 1:n_parameters(sp.sf)
        ]
    end

    # get adapted p_space
    adapted_p_space = [
        i ∈ dims ? p_space[i] : _single_point_range(sp.optimal_ps[i])
        for i in 1:n_parameters(sp.sf)
    ]

    # generate parameter combinations
    p_combos = _parameter_combinations(adapted_p_space)

    # set up arrays to store residual landscapes
    residuals = zeros(size(p_combos))
    for (i, ps) in enumerate(p_combos)
        residuals[i] = sp.quality(sp, ps)
    end

    return p_space[dims], reshape(residuals, [size(residuals)[i] for i in sort(dims)]...)
end

# range containing only x; unlike x:x, this does not search for a rational approximation
# of x, which can overflow (DivideError) on 32-bit systems
_single_point_range(x) = StepRangeLen(Float64(x), 0.0, 1)
