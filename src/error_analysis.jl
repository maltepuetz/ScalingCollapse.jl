function error_analysis!(sp)

    sp.verbose && @info "Starting error analysis..."

    for i in eachindex(sp.optimal_ps)

        # define the function to find roots of
        function f_root(x)
            y = sp.quality(
                sp,
                sp.optimal_ps .+ x .* _one_vector(i, length(sp.optimal_ps))
            ) / sp.minimum - sp.error_threshold[i]
            return y
        end

        # find starting intervals
        l, r = _start_deltas(sp.optimal_ps[i], f_root)

        # find roots
        delta_l = _root(l, f_root(l), 0.0, -1.0, f_root)
        delta_r = _root(0.0, -1.0, r, f_root(r), f_root)

        # set error to be the largest of l and r error
        sp.optimal_ps_error[i] = max(abs(delta_l), abs(delta_r))
    end
    sp.verbose && @info "Found errors."

    return nothing
end

function _one_vector(i, L)
    res = zeros(L)
    res[i] = 1.0
    return res
end

function _start_deltas(optimal_pi, f::Function)

    l = -0.1 * optimal_pi
    while f(l) < 0.0
        l *= 2.0
    end

    r = 0.1 * optimal_pi
    while f(r) < 0.0
        r *= 2.0
    end

    return l, r
end

function _root(l, fl, r, fr, f::Function; precision=1e-10)

    # iterate until relative precision is reached
    m = (l + r) / 2.0
    while abs((r - l) / m) > precision
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
    end

    return m
end
