function residual_landscapes(sp::ScalingProblem; kwargs...)

    if sp.optimal_ps_error == zeros(size(sp.optimal_ps)) && !haskey(kwargs, p_space)
        error("Either run error analysis or provide p_space to plot residual landscapes.")
    end

    # set parameter space
    N_steps = get(kwargs, :N_steps, 50)
    p_space = get(
        kwargs,
        :p_space,
        [
            range(
                sp.optimal_ps[i] - 3 * sp.optimal_ps_error[i],
                sp.optimal_ps[i] + 3 * sp.optimal_ps_error[i],
                length=N_steps
            ) for i in 1:n_parameters(sp.sf)
        ]
    )

    # set up arrays to store residual landscapes
    residuals = [zeros(length(p_space[i])) for i in 1:n_parameters(sp.sf)]
    for i in 1:n_parameters(sp.sf)
        for (j, p) in enumerate(p_space[i])
            ps = copy(sp.optimal_ps)
            ps[i] = p
            residuals[i][j] = sp.quality(sp, ps)
        end
    end

    return p_space, residuals
end
