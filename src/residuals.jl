"""
    residuals(sp::ScalingProblem; p_space, N_steps=51, dims) -> (p_space, S)

Evaluate the quality function `sp.quality` on a grid of free parameters, e.g. to plot the
quality landscape around a fit or to check that a minimum is not a numerical artifact.

Returns the ranges of the active dimensions and the array `S` of quality values, with one
axis per active dimension: `S[i, j, ...]` belongs to `p_space[1][i]`, `p_space[2][j]`, ...
The parameters of the other dimensions stay at their values in `sp.optimal_ps`. The
parameters are not checked against the `p_space` of `sp`.

# Keyword arguments
- `dims::Vector{Int}=1:n`: the indices of the active free parameters (in the order of
    `ScalingCollapse.scaled_p_names(sp.sf)`), e.g. `dims=[1, 2]` for a 2D cut through a
    problem with three free parameters.
- `p_space`: one range per free parameter; the ranges of inactive dimensions are ignored.
    By default, every active parameter spans `optimal_ps ± 3 optimal_ps_error` with
    `N_steps` points, which needs finite, nonzero errors.
- `N_steps::Int=51`: the number of points per dimension of the default `p_space`.

# Examples
```julia
sp = ScalingProblem(Ts, chi, Ls; sf=ScalingFunction(:xny; p_names=["T_c", "nu", "gamma"]))

# 3D landscape: length(ps) == 3, size(S) == (51, 51, 51)
ps, S = residuals(sp)

# 2D cut in T_c and nu, with gamma at its optimum
ps, S = residuals(sp; dims=[1, 2], N_steps=101)

# 1D cut along nu on a given range
ps, S = residuals(sp; dims=[2], p_space=[0:0, 0.8:0.005:1.2, 0:0])
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
