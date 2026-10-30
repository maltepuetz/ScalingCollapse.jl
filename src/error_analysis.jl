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

function _profile_errors!(sp; max_restarts::Int=3)
    # TODO Task 3
    return _conditional_errors!(sp)
end

function _conditional_errors!(sp)

    sp.verbose && @info "Starting error analysis..."

    p_names = scaled_p_names(sp.sf)

    # the error analysis compares the quality relative to its minimum
    if !(isfinite(sp.minimum) && sp.minimum > 0.0)
        @warn "Skipping error analysis: the minimum of the quality function is " *
              "$(sp.minimum), but a finite, positive minimum is required. The " *
              "parameter errors are set to NaN."
        sp.optimal_ps_error .= NaN
        return nothing
    end

    for i in eachindex(sp.optimal_ps)

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

        if any(isnothing, deltas)
            @warn "Error analysis failed for parameter $(p_names[i]): the quality " *
                  "function does not reach error_threshold = $(sp.error_threshold[i]) " *
                  "times its minimum in both directions (or becomes NaN). The error is " *
                  "set to NaN."
            sp.optimal_ps_error[i] = NaN
            continue
        end

        # set error to be the largest of l and r error
        sp.optimal_ps_error[i] = maximum(abs, deltas)
    end
    sp.verbose && @info "Found errors."

    return nothing
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
