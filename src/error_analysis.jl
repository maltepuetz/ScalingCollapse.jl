function error_analysis!(sp; max_restarts::Int=3)
    if sp.error === :conditional
        _conditional_errors!(sp)
    elseif sp.error === :profile
        _profile_errors!(sp; max_restarts)
    else
        throw(ArgumentError("unknown error analysis method $(repr(sp.error))."))
    end
    return nothing
end

# Profile errors (`error = :profile`): for parameter i, the profile quality at displacement x
# is the quality minimized over the other ("nuisance") parameters with p[i] = optimal_ps[i] + x.
# Starting at the conditional roots, each side is bracketed outward and bisected to the first
# valid finite point where profile / minimum reaches error_threshold[i].
const _PROFILE_GROWTH = 1.5             # outward bracket growth factor
const _PROFILE_PRECISION = 1e-3         # relative bisection precision
const _PROFILE_MAX_STEPS = 60           # outward bracket steps per side
const _PROFILE_MAX_BISECTIONS = 100
const _PROFILE_LOWER_TOLERANCE = 1e-6   # relative; a lower quality means a better fit exists
const _PROFILE_MIN_START = 1e-6         # relative floor of the first outward bracket step

# bookkeeping of the profile evaluations of one parameter
mutable struct _ProfileRecord
    n_failed::Int                   # inner fits that did not converge
    lowest_q::Float64               # lowest quality seen inside p_space ...
    lowest_ps::Vector{Float64}      # ... and its full parameter vector
    lower_outside::Bool             # a lower quality was seen outside p_space
    failures::Vector{Tuple{Symbol,String}}  # (reason, side) of each side without a width
end
_ProfileRecord(sp) = _ProfileRecord(0, sp.minimum, copy(sp.optimal_ps), false, [])

# whether q is lower than the minimum by more than the tolerance: relative, but never below
# the absolute accuracy at which Nelder-Mead stops (g_abstol), so optimizer noise is ignored
_is_lower(sp, q) = q < sp.minimum -
    max(_PROFILE_LOWER_TOLERANCE * sp.minimum, sp.optimizer_options.g_abstol)

# whether a profile evaluation found a better fit inside p_space
_found_lower(sp, record::_ProfileRecord) = _is_lower(sp, record.lowest_q)

_inside_p_space(sp, ps) =
    all(j -> first(sp.p_space[j]) <= ps[j] <= last(sp.p_space[j]), eachindex(ps))

function _warn_profile_failure(p_name, threshold, reason, sides)
    reach = "error_threshold = $(threshold) times its minimum"
    problem = if reason === :optimum
        "the quality at the optimum is not below $(reach)"
    elseif reason === :steps
        "the profile quality does not reach $(reach) $(sides) the optimum within " *
        "$(_PROFILE_MAX_STEPS) steps"
    else # :undefined
        "the profile quality becomes undefined (no overlap of the rescaled data within dx) " *
        "$(sides) the optimum before reaching $(reach)"
    end
    @warn "Error analysis failed for parameter $(p_name): $(problem). The error is set to NaN."
    return nothing
end

function _profile_errors!(sp; max_restarts::Int=3)

    # with one free parameter, the profile is the plain quality
    length(sp.optimal_ps) == 1 && return _conditional_errors!(sp)

    sp.verbose && @info "Starting error analysis..."
    _has_valid_minimum(sp) || return nothing

    # stop at the first parameter whose profile finds a better fit inside p_space
    lower = nothing
    for i in eachindex(sp.optimal_ps)
        sp.optimal_ps_error[i], record = _profile_error(sp, i)
        if _found_lower(sp, record)
            lower = record
            break
        end
    end

    if !isnothing(lower)
        old_minimum = sp.minimum
        if max_restarts > 0
            # re-optimize from the better point and redo the analysis from scratch
            _reoptimize!(sp, lower.lowest_ps, lower.lowest_q)
            @warn "The error analysis found a lower quality than the minimum: the fit is " *
                  "re-optimized, which moves the minimum from $(old_minimum) to " *
                  "$(sp.minimum) (optimal parameters $(sp.optimal_ps)), and the error " *
                  "analysis is restarted."
            return _profile_errors!(sp; max_restarts=max_restarts - 1)
        end
        sp.optimal_ps = lower.lowest_ps
        sp.minimum = lower.lowest_q
        sp.optimal_ps_error .= NaN
        @warn "The minimum could not be confirmed: the error analysis found a lower quality " *
              "than the minimum $(old_minimum) and no restarts are left. The optimal " *
              "parameters are set to the best fit found, $(sp.optimal_ps) with quality " *
              "$(sp.minimum), and all parameter errors are set to NaN."
        return nothing
    end
    sp.verbose && @info "Found errors."

    return nothing
end

# final optimization started from ps (with quality q); sp.starting_ps keeps its value. The
# result is never worse than ps itself.
function _reoptimize!(sp, ps, q)
    starting_ps = sp.starting_ps
    sp.starting_ps = ps
    try
        optimize_starting_ps!(sp)
    finally
        sp.starting_ps = starting_ps
    end
    if sp.minimum > q
        sp.optimal_ps, sp.minimum = copy(ps), q
    end
    return nothing
end

# profile error of parameter i and the record of its profile evaluations
function _profile_error(sp, i::Int)
    record = _ProfileRecord(sp)

    # the conditional roots are the starting points of the profile brackets
    deltas, f0 = _conditional_roots(sp, i)
    any(isnothing, deltas) && return _no_crossing_error(sp, i), record

    widths = Union{Nothing,Float64}[]
    for d in deltas
        push!(widths, _profile_side(sp, i, d, f0; record))
        # the analysis is redone from a better fit (see _profile_errors!)
        _found_lower(sp, record) && return NaN, record
    end

    p_name = scaled_p_names(sp.sf)[i]
    if record.lower_outside
        @warn "Error analysis failed for parameter $(p_name): the profile quality is lower " *
              "than the minimum for $(p_name) outside p_space. A better fit may exist " *
              "outside p_space, consider widening p_space. The error is set to NaN."
        return NaN, record
    end

    # one warning per parameter and reason
    for reason in unique(first.(record.failures))
        sides = join([side for (r, side) in record.failures if r === reason], " and ")
        _warn_profile_failure(p_name, sp.error_threshold[i], reason, sides)
    end
    if record.n_failed > 0
        @warn "The error of parameter $(p_name) may be underestimated: " *
              "$(record.n_failed) fits of the other parameters (profile quality) did not " *
              "converge."
    end
    any(isnothing, widths) && return NaN, record
    return maximum(abs, widths), record
end

# Displacement (signed) of the first valid finite point outward from the fit, on the side of
# the conditional root d, where the profile quality reaches error_threshold[i] times the
# minimum. f0 = quality / minimum - error_threshold[i] at the fit. Returns nothing if no such
# point is found, and adds the reason to record.failures.
function _profile_side(sp, i::Int, d::Float64, f0::Float64;
    record::_ProfileRecord=_ProfileRecord(sp))

    side = d < 0 ? "below" : "above"
    threshold = sp.error_threshold[i]
    nuisance = _nuisance_indices(sp, i)
    if !(f0 < 0.0)
        push!(record.failures, (:optimum, side))
        return nothing
    end

    # evaluated displacements and their fitted nuisance values (warm starts)
    visited = [(0.0, sp.optimal_ps[nuisance])]

    # profile / minimum - threshold; NaN if the profile quality is invalid
    function f(x)
        warm = last(visited[argmin([abs(v[1] - x) for v in visited])])
        q, nu, converged = _profile_value(sp, i, x, warm)
        converged || (record.n_failed += 1)
        isfinite(q) || return NaN
        push!(visited, (x, nu))
        if q < sp.minimum
            ps = copy(sp.optimal_ps)
            ps[i] += x
            ps[nuisance] .= nu
            if !_inside_p_space(sp, ps)
                record.lower_outside |= _is_lower(sp, q)
            elseif q < record.lowest_q
                record.lowest_q, record.lowest_ps = q, ps
            end
        end
        return q / sp.minimum - threshold
    end

    # bracket outward from the conditional root (at least a small step relative to the
    # parameter, e.g. when the minimum sits on a jump of the quality); lo stays a valid point
    # below threshold, hi is a valid point at or above threshold (hi_valid) or an invalid point
    lo = 0.0
    x = copysign(max(abs(d), _PROFILE_MIN_START * max(abs(sp.optimal_ps[i]), 1.0)), d)
    hi, hi_valid = NaN, false
    for _ in 1:_PROFILE_MAX_STEPS
        fx = f(x)
        if isnan(fx) || fx >= 0.0
            hi, hi_valid = x, !isnan(fx)
            break
        end
        lo = x
        x *= _PROFILE_GROWTH
    end
    if isnan(hi)
        push!(record.failures, (:steps, side))
        return nothing
    end

    # bisect; invalid midpoints become the outer end, so the result is a valid crossing
    # before the quality becomes undefined
    for _ in 1:_PROFILE_MAX_BISECTIONS
        m = (lo + hi) / 2.0
        if abs((hi - lo) / m) <= _PROFILE_PRECISION
            hi_valid && return m
            break
        end
        fm = f(m)
        if isnan(fm) || fm >= 0.0
            hi, hi_valid = m, !isnan(fm)
        else
            lo = m
        end
    end
    push!(record.failures, (:undefined, side))
    return nothing
end

# Profile quality of parameter i at displacement x: the quality minimized over the nuisance
# parameters (kept inside their p_space bounds), starting from the better of the conditional
# candidate (nuisance at the optimum) and the warm start. Returns the quality (Inf if invalid),
# the nuisance values and whether the inner fit converged. The quality is never above either
# candidate's quality, so the profile width is never below the conditional width.
function _profile_value(sp, i::Int, x::Float64, warm::Vector{Float64})
    nuisance = _nuisance_indices(sp, i)
    lower = Float64[first(sp.p_space[j]) for j in nuisance]
    upper = Float64[last(sp.p_space[j]) for j in nuisance]
    ps = copy(sp.optimal_ps)
    ps[i] += x
    function f(nu)
        all(k -> lower[k] <= nu[k] <= upper[k], eachindex(nu)) || return Inf
        ps[nuisance] .= nu
        return _nan_to_inf(sp.quality(sp, ps; check_bounds=false))
    end

    conditional = sp.optimal_ps[nuisance]
    q_conditional, q_warm = f(conditional), f(warm)
    start_conditional = q_conditional <= q_warm
    q, nu = start_conditional ? (q_conditional, conditional) : (q_warm, copy(warm))
    # no finite starting point: the point is invalid
    isfinite(q) || return Inf, nu, true

    optimizer = Optim.NelderMead(
        BoundedSimplexer(sp.optimizer.initial_simplex, lower, upper, f),
        sp.optimizer.parameters,
    )
    result = optimize(f, copy(nu), optimizer, sp.optimizer_options)
    converged = Optim.converged(result)
    if !converged
        # retry once from the conditional candidate (from where the first fit stopped, if it
        # already started there)
        retry_start = start_conditional ? copy(Optim.minimizer(result)) : copy(conditional)
        if Optim.minimum(result) < q
            q, nu = Optim.minimum(result), copy(Optim.minimizer(result))
        end
        result = optimize(f, retry_start, optimizer, sp.optimizer_options)
        converged = Optim.converged(result)
    end
    if Optim.minimum(result) < q
        q, nu = Optim.minimum(result), copy(Optim.minimizer(result))
    end

    return q, nu, converged
end

_nuisance_indices(sp, i) = [j for j in eachindex(sp.optimal_ps) if j != i]

function _conditional_errors!(sp)

    sp.verbose && @info "Starting error analysis..."
    _has_valid_minimum(sp) || return nothing

    for i in eachindex(sp.optimal_ps)
        deltas, _ = _conditional_roots(sp, i)
        sp.optimal_ps_error[i] = any(isnothing, deltas) ? _no_crossing_error(sp, i) :
                                 maximum(abs, deltas)
    end
    sp.verbose && @info "Found errors."

    return nothing
end

# the error analysis compares the quality relative to its minimum
function _has_valid_minimum(sp)
    isfinite(sp.minimum) && sp.minimum > 0.0 && return true
    @warn "Skipping error analysis: the minimum of the quality function is " *
          "$(sp.minimum), but a finite, positive minimum is required. The " *
          "parameter errors are set to NaN."
    sp.optimal_ps_error .= NaN
    return false
end

# Roots (d₋, d₊) of the quality with all other parameters at their optimum (nothing for a side
# without a root), and f_root(0.0) = quality / minimum - error_threshold[i] at the fit.
function _conditional_roots(sp, i)

    # define the function to find roots of
    function f_root(x)
        y = sp.quality(
            sp,
            sp.optimal_ps .+ x .* _one_vector(i, length(sp.optimal_ps))
        ) / sp.minimum - sp.error_threshold[i]
        return y
    end

    # f_root(0) < 0 since error_threshold > 1
    f0 = f_root(0.0)

    # initial step relative to the parameter (absolute if the parameter is zero)
    step = iszero(sp.optimal_ps[i]) ? 0.1 : 0.1 * abs(sp.optimal_ps[i])

    # find starting intervals and roots on both sides
    deltas = map((-step, step)) do s
        bracket = _bracket(s, f_root)
        isnothing(bracket) && return nothing
        x, fx = bracket
        return _root(x, fx, 0.0, f0, f_root)
    end

    return deltas, f0
end

function _no_crossing_error(sp, i)
    @warn "Error analysis failed for parameter $(scaled_p_names(sp.sf)[i]): the quality " *
          "function does not reach error_threshold = $(sp.error_threshold[i]) " *
          "times its minimum in both directions (or becomes NaN). The error is " *
          "set to NaN."
    return NaN
end

function _one_vector(i, L)
    res = zeros(L)
    res[i] = 1.0
    return res
end

# find x (in the direction of `step`) with f(x) >= 0; returns (x, f(x)) or nothing
function _bracket(step, f::Function; max_shrink=100, max_grow=100)

    # decrease starting step while function returns NaN
    fx = f(step)
    n_shrink = 0
    while isnan(fx)
        n_shrink += 1
        n_shrink > max_shrink && return nothing
        step *= 0.1^0.1
        fx = f(step)
    end

    # increase step until the threshold is crossed
    n_grow = 0
    while fx < 0.0
        n_grow += 1
        n_grow > max_grow && return nothing
        step *= 2.0
        fx = f(step)
        isnan(fx) && return nothing
    end

    return step, fx
end


function _root(l, fl, r, fr, f::Function; precision=1e-10, max_iterations=200)

    # iterate until relative precision is reached
    m = (l + r) / 2.0
    n = 0
    while abs((r - l) / m) > precision && n < max_iterations
        # replace left or right border by middle point
        fm = f(m)
        if fm * fl > 0.0
            l = m
            fl = fm
        else
            r = m
            fr = fm
        end
        m = (l + r) / 2.0
        n += 1
    end

    return m
end
