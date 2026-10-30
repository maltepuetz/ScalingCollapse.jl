"""
    ScalingProblem(args...; kwargs...)

Create a scaling problem which solves on initialization.

# Arguments
- `xs`: x values of the data
- `ys`: y values of the data
- `es`: y error values of the data (optional)
- `Ls`: system sizes of the data
For more information on the arguments, see methods(ScalingCollapse.unzip_data).

# Keyword Arguments
- `sf::ScalingFunction=ScalingFunction(preset; kwargs...)`: scaling function
- `p_space::Vector{StepRangeLen}=[0.1:0.1:3.0 for _ in sf.p_names]`: parameter search space
- `dx::Vector{Float64}=[-Inf, Inf]`: optimization interval
- `quality::QualityFunction=Linear()`: quality function
- `quality_scan::QualityFunction=MultipleSplines(scan_mode=true)`: quality function for the
    parameter space scan (it is highly recommended to use the default function here)
- `verbose::Bool=false`: print information during optimization
- `starting_ps::Vector{Float64}`: If `starting_ps` are given, there will be no initial
    parameter space scan. Instead, the optimization will start at the given points. This is
    much faster, but might not find the global minimum.
- `error::Union{Bool,Symbol}=true`: error analysis method. `true`/`:profile` (default),
    `:conditional`, `:bootstrap`, or `false` (no error analysis); anything else throws an
    `ArgumentError`.
    With `:profile`, the error of a free parameter is a local profile width: the parameter
    is displaced from its optimum, all other free parameters are re-fitted (inside
    `p_space`) at every displacement, and the error is the larger of the two displacements
    (one on each side) at which the re-fitted quality first reaches `error_threshold` times
    its minimum when moving outward from the fit. The profiled parameter itself is not
    bounded by `p_space`. With `:conditional`, the other parameters are held at their
    optimal values instead, which ignores correlations between parameters and gives smaller
    or equal errors (equal with one free parameter or uncorrelated parameters). With one free parameter both methods agree. An error is `NaN` if the threshold
    is never reached or the quality becomes undefined (no overlap within `dx`) first. If the
    profile finds a quality lower than the minimum inside `p_space`, the fit is re-optimized
    from there and the analysis restarts (up to 3 times, with a warning), so `optimal_ps`
    and `minimum` can change during the error analysis. If the restarts run out, all errors
    are `NaN`. A lower quality outside `p_space` gives `NaN` for that parameter and a
    warning that `p_space` may need widening. With `error=:profile`, `optimizer` and
    `optimizer_options` (including `time_limit`, `show_trace` and `callback`) also apply to
    every inner fit of the error analysis. The errors are widths defined by a quality
    threshold, not standard deviations with a known coverage. Profiling accounts for
    correlations between parameters, but the threshold is not calibrated and correlated
    data points (e.g. curves built from the same Monte Carlo samples) are not accounted for.
    With `:bootstrap`, the error of a free parameter is a standard deviation: `n_bootstrap`
    synthetic datasets are drawn by replacing every y value with a Gaussian random number
    centred on it with its error bar as width, each one is re-fitted (Nelder-Mead from
    `optimal_ps`, inside `p_space`, with `optimizer` and `optimizer_options`), and the error
    is the standard deviation of the re-fitted parameter. This requires error bars (`es` or
    `Measurement` y values) and assumes they are Gaussian and independent; correlated data
    points are not accounted for. `optimal_ps` and `minimum` remain those of the original
    data, and `error_threshold` is not used. Re-fits without a finite quality are dropped
    (with a warning); if more than 10% are dropped, all errors are `NaN`.
- `error_threshold::Union{Float64,Vector{Float64}}=1.3`: The error of a parameter is the
    largest deviation from its optimal value, with the other free parameters re-fitted
    (`:profile`) or fixed (`:conditional`), at which the quality function reaches
    `error_threshold` times its minimum. Must be larger than 1.
- `n_bootstrap::Int=100`: number of synthetic datasets for `error=:bootstrap` (at least 2).
    The bootstrap errors scatter between random seeds by a relative amount of about
    `1/sqrt(2(n_bootstrap - 1))` (about 7% for 100), which is enough to quote errors to one
    significant digit; use about 500 or more for a reliable second digit.
- `rng::AbstractRNG=Random.default_rng()`: random number generator for `error=:bootstrap`;
    pass a seeded one (e.g. `Random.Xoshiro(42)`) for reproducible errors.
- `optimizer::Optim.NelderMead`: Nelder-Mead method configuration, including the
    initial simplex and adaptive or fixed algorithm parameters.
- `optimizer_options::Optim.Options`: Optimization runtime and stopping configuration,
    including tolerances, iteration and function-call limits, callbacks, tracing, and
    warnings.

# Fields
- `data::Vector{Data}`: data to be scaled
- `sf::ScalingFunction`: scaling function
- `p_space::Vector{StepRangeLen}`: search parameter space
- `dx::Vector{Float64}`: optimization interval
- `optimal_ps::Vector{Float64}`: optimal parameters
- `optimal_ps_error::Vector{Float64}`: uncertainties of the optimal parameters
- `minimum::Float64`: minimum of the objective function
- `solved::Bool`: `true` if optimization was called
- `verbose::Bool`: print information during optimization
- `error::Symbol`: error analysis method (`:profile`, `:conditional`, `:bootstrap` or `:none`)
- `quality_scan::QualityFunction`: quality function for the parameter space scan
- `quality::QualityFunction`: quality function
- `skip_scan::Bool`: `true` if `starting_ps` are given by the user
- `starting_ps::Vector{Float64}`: starting points for the optimization
- `errors_defined::Bool`: `true` if errors are given by the user
- `optimizer::Optim.NelderMead`: Nelder-Mead method configuration
- `optimizer_options::Optim.Options`: optimization runtime and stopping configuration
- `n_bootstrap::Int`: number of synthetic datasets for `error=:bootstrap`
- `rng::AbstractRNG`: random number generator for `error=:bootstrap`

# Examples
```julia
using ScalingCollapse
```

### rescale the x axis only
```julia
# (e.g. for ys == binder cumulant of Ising model)
ScalingProblem(xs, ys, Ls;
    sf=ScalingFunction(:x; p_names=["T_c", "nu"]),
    p_space=[0.1:0.1:3.0, 0.1:0.1:3.0],
    dx=[-1.0, 1.0],
)

# or the same but with errors
ScalingProblem(xs, ys, es, Ls;
    sf=ScalingFunction(:x; p_names=["T_c", "nu"]),
    p_space=[0.1:0.1:3.0, 0.1:0.1:3.0],
    dx=[-1.0, 1.0],
)
```

### rescale the x and y axis
```julia
# (e.g. for ys == magnetization of Ising model)
ScalingProblem(xs, ys, Ls;
    sf=ScalingFunction(:xy; p_names=["T_c", "nu", "beta"]),
    dx=[-1.0, 1.0],
)
```

### rescale the x and y axis, but y axis with negative exponent
```julia
# (e.g. for ys == susceptibility of Ising model)
ScalingProblem(xs, ys, Ls;
    sf=ScalingFunction(:xny; p_names=["T_c", "nu", "gamma"]),
    dx=[-1.0, 1.0],
)
```

### configure the Nelder-Mead optimizer
```julia
optimizer = ScalingCollapse.Optim.NelderMead(
    initial_simplex=ScalingCollapse.Optim.AffineSimplexer(a=0.01, b=0.1),
)
optimizer_options = ScalingCollapse.Optim.Options(
    iterations=5_000,
    x_abstol=1e-10,
    f_reltol=1e-12,
    time_limit=60.0,
)

ScalingProblem(xs, ys, Ls;
    sf=ScalingFunction(:x; p_names=["T_c", "nu"]),
    starting_ps=[2.27, 1.0],
    optimizer=optimizer,
    optimizer_options=optimizer_options,
)
```

"""
mutable struct ScalingProblem

    # data to be scaled
    data::Vector{Data}

    # scaling function, parameter space, and optimization interval
    sf::ScalingFunction
    p_space::Vector{StepRangeLen}
    dx::Vector{Float64}

    # optimization results
    optimal_ps::Vector{Float64}
    optimal_ps_error::Vector{Float64}
    minimum::Float64

    # other options
    solved::Bool
    verbose::Bool
    error::Symbol
    quality_scan::QualityFunction
    quality::QualityFunction
    skip_scan::Bool
    starting_ps::Vector{Float64}
    errors_defined::Bool
    error_threshold::Vector{Float64}
    optimizer::Optim.NelderMead
    optimizer_options::Optim.Options
    n_bootstrap::Int
    rng::AbstractRNG



    function ScalingProblem(
        data::Vector{Data},
        sf,
        p_space,
        dx,
        optimal_ps,
        optimal_ps_error,
        minimum,
        solved,
        verbose,
        error,
        quality_scan,
        quality,
        skip_scan,
        starting_ps,
        errors_defined,
        error_threshold,
        optimizer,
        optimizer_options,
        n_bootstrap,
        rng,
    )
        return new(
            data,
            sf,
            p_space,
            dx,
            optimal_ps,
            optimal_ps_error,
            minimum,
            solved,
            verbose,
            error,
            quality_scan,
            quality,
            skip_scan,
            starting_ps,
            errors_defined,
            error_threshold,
            optimizer,
            optimizer_options,
            n_bootstrap,
            rng,
        )
    end

    function ScalingProblem(args...; kwargs...)

        ##### unzip kwargs... #####
        if haskey(kwargs, :sf)
            _check_kwargs(
                "ScalingProblem",
                keys(kwargs),
                SCALINGPROBLEM_KWARGS;
                hint="Keyword arguments of the ScalingFunction (e.g. p_names or fixed " *
                     "parameters) must be passed to ScalingFunction(...) when `sf` is given.",
            )
            sf = kwargs[:sf]
        else
            # remaining kwargs are passed to (and checked by) the ScalingFunction
            sf_kwargs = (; (k => v for (k, v) in kwargs if k ∉ SCALINGPROBLEM_KWARGS)...)
            sf = try
                ScalingFunction(; sf_kwargs...)
            catch e
                e isa ArgumentError || rethrow()
                throw(ArgumentError(
                    "The keyword arguments $(join(keys(sf_kwargs), ", ")) are not " *
                    "ScalingProblem keywords, so they were passed to the default " *
                    "ScalingFunction, which rejected them: $(e.msg) " *
                    "Valid ScalingProblem keyword arguments are: " *
                    "$(join(SCALINGPROBLEM_KWARGS, ", "))."
                ))
            end
        end
        sf isa ScalingFunction || throw(ArgumentError(
            "sf must be a ScalingFunction, got $(typeof(sf))."
        ))
        n_ps = n_parameters(sf)
        p_space = get(kwargs, :p_space, [0.1:0.1:3.0 for _ in 1:n_ps])
        dx = get(kwargs, :dx, [-Inf, Inf])
        verbose = get(kwargs, :verbose, false)
        error = _error_method(get(kwargs, :error, true))
        skip_scan = false
        starting_ps = get(kwargs, :starting_ps, zeros(n_ps))
        if haskey(kwargs, :starting_ps)
            length(starting_ps) == n_ps || throw(ArgumentError(
                "starting_ps has $(length(starting_ps)) entries, but the scaling function " *
                "has $(n_ps) free parameters $(scaled_p_names(sf))."
            ))
            # if starting_ps are given, there will be no initial parameter space scan
            p_space = [-1_000_000:1_000_000 for _ in 1:n_ps]
            skip_scan = true
        end
        quality_scan = get(kwargs, :quality_scan, MultipleSplines(scan_mode=true))
        quality = get(kwargs, :quality, Linear())
        error_threshold = get(kwargs, :error_threshold, 1.3 .* ones(n_ps))
        if !isa(error_threshold, Vector)
            error_threshold = error_threshold .* ones(n_ps)
        end
        length(error_threshold) == n_ps || throw(ArgumentError(
            "error_threshold has $(length(error_threshold)) entries, but the scaling " *
            "function has $(n_ps) free parameters $(scaled_p_names(sf))."
        ))
        all(>(1), error_threshold) || throw(ArgumentError(
            "error_threshold must be larger than 1 (it is the factor by which the " *
            "quality function may increase), got $(error_threshold)."
        ))
        (length(dx) == 2 && dx[1] < dx[2]) || throw(ArgumentError(
            "dx must be an interval [x_min, x_max] with x_min < x_max, got $(dx)."
        ))
        optimizer = get(
            kwargs,
            :optimizer,
            NelderMead(; initial_simplex=Optim.AffineSimplexer(; b=0.1)),
        )
        optimizer isa Optim.NelderMead || throw(
            ArgumentError("optimizer must be an Optim.NelderMead instance"),
        )
        optimizer_options = get(kwargs, :optimizer_options, Optim.Options())
        optimizer_options isa Optim.Options || throw(
            ArgumentError("optimizer_options must be an Optim.Options instance"),
        )
        n_bootstrap = get(kwargs, :n_bootstrap, 100)
        (n_bootstrap isa Integer && n_bootstrap >= 2) || throw(ArgumentError(
            "n_bootstrap must be an integer of at least 2, got $(repr(n_bootstrap))."
        ))
        rng = get(kwargs, :rng, Random.default_rng())
        rng isa AbstractRNG || throw(ArgumentError(
            "rng must be an AbstractRNG (e.g. Random.Xoshiro(seed)), got $(typeof(rng))."
        ))

        length(p_space) == n_ps || throw(ArgumentError(
            "p_space has $(length(p_space)) entries, but the scaling function has " *
            "$(n_ps) free parameters $(scaled_p_names(sf))" *
            (isempty(sf.fixed_idxs) ? "." : " (fixed: $(fixed_p_names(sf)))."
            )
        ))
        all(!isempty, p_space) || throw(ArgumentError("p_space contains an empty range."))
        for (p, r) in zip(scaled_p_names(sf), p_space)
            (issorted(r) && allunique(r)) || throw(ArgumentError(
                "p_space for $p must be strictly increasing, got $r."
            ))
        end

        # unzip args...
        data = unzip_data(args...)
        _check_data(data, quality, quality_scan, skip_scan)
        _check_scaling_function(sf, first(data), skip_scan ? starting_ps : first.(p_space))

        # check if errors are defined
        errors_defined = false
        for d in data
            if any(d.es .!= 0.0)
                errors_defined = true
                break
            end
        end
        (error === :bootstrap && !errors_defined) && throw(ArgumentError(
            "error=:bootstrap resamples the y values within their error bars, but the data " *
            "have no error bars. Pass y errors `es` or Measurement y values."
        ))

        # construct unsolved ScalingProblem
        optimal_ps = zeros(n_parameters(sf))
        optimal_ps_error = zeros(n_parameters(sf))
        minimum = 0.0
        solved = false
        sp = ScalingProblem(
            data,
            sf,
            p_space,
            dx,
            optimal_ps,
            optimal_ps_error,
            minimum,
            solved,
            verbose,
            error,
            quality_scan,
            quality,
            skip_scan,
            starting_ps,
            errors_defined,
            error_threshold,
            optimizer,
            optimizer_options,
            n_bootstrap,
            rng,
        )

        # solve ScalingProblem
        solve!(sp)
        return sp
    end

end


const SCALINGPROBLEM_KWARGS = (
    :sf,
    :p_space,
    :dx,
    :verbose,
    :error,
    :starting_ps,
    :quality_scan,
    :quality,
    :error_threshold,
    :optimizer,
    :optimizer_options,
    :n_bootstrap,
    :rng,
)

function _check_data(data, quality, quality_scan, skip_scan)
    length(data) >= 2 || throw(ArgumentError(
        "A scaling collapse needs data for at least two system sizes, got $(length(data))."
    ))
    for d in data
        for (name, vals) in (("x", d.xs), ("y", d.ys), ("error", d.es))
            all(isfinite, vals) || throw(ArgumentError(
                "The $name values for L = $(d.L) contain non-finite values (NaN or Inf)."
            ))
        end
        # xs are sorted, so repeated values are adjacent
        i = findfirst(i -> d.xs[i] == d.xs[i+1], 1:length(d.xs)-1)
        isnothing(i) || throw(ArgumentError(
            "The x value $(d.xs[i]) appears more than once for L = $(d.L). Every system " *
            "size must have distinct x values; average or remove the duplicates."
        ))
    end

    # splines of degree 3 need at least 4 data points per system size
    uses_splines = quality isa MultipleSplines || (!skip_scan && quality_scan isa MultipleSplines)
    if uses_splines
        for d in data
            length(d.xs) >= 4 || throw(ArgumentError(
                "The MultipleSplines quality function needs at least 4 data points per " *
                "system size, but L = $(d.L) has $(length(d.xs)). Use a different " *
                "quality / quality_scan function or provide starting_ps to skip the scan."
            ))
        end
    end
    return nothing
end

# call the scaling function once to catch signature and return type mistakes early
function _check_scaling_function(sf, d, ps)
    signature = "f(d::ScalingCollapse.Data, $(join(sf.p_names, ", "))) -> ScalingCollapse.Data"
    scaled = try
        sf.f(d, Float64.(ps)...)
    catch e
        # other exceptions may be specific to these parameters and surface later
        e isa MethodError || return nothing
        throw(ArgumentError(
            "Calling the scaling function failed with a MethodError. It must have the " *
            "signature $signature. Original error: $(sprint(showerror, e))"
        ))
    end
    scaled isa Data || throw(ArgumentError(
        "The scaling function returned a $(typeof(scaled)), but it must return a " *
        "ScalingCollapse.Data object: $signature."
    ))
    return nothing
end


function solve!(sp::ScalingProblem)

    # find global minimum
    optimize_parameters!(sp)

    # perform error analysis
    if sp.error !== :none
        error_analysis!(sp)
    end
end

# map the user-facing `error` keyword to the internal method symbol
function _error_method(x)
    (x === true || x === :profile) && return :profile
    x === :conditional && return :conditional
    x === :bootstrap && return :bootstrap
    x === false && return :none
    throw(ArgumentError(
        "error must be one of true, false, :profile, :conditional, :bootstrap, got " *
        "$(repr(x))."
    ))
end


function Base.show(io::IO, sp::ScalingProblem)
    println(io, "ScalingProblem")
    println(io, "    Scaling function:")
    println(io, "        x -> $(sp.sf.x_scale)")
    println(io, "        y -> $(sp.sf.y_scale)")
    if n_parameters(sp.sf) < length(sp.sf.p_names)
        println(io, "    Fixed parameters:")
        for i in eachindex(sp.sf.fixed_idxs)
            idx = sp.sf.fixed_idxs[i]
            p = sp.sf.p_names[idx]
            println(io, "        $p = $(sp.sf.fixed_ps[i])")
        end
    end
    println(io, "    Optimal parameters:")
    for (i, p) in enumerate(scaled_p_names(sp.sf))
        println(
            io,
            "        $p = $(sp.optimal_ps[i])" *
            (
                if sp.optimal_ps_error[i] == 0.0
                    ""
                else
                    " ± $(sp.optimal_ps_error[i])"
                end
            )
        )
    end
    println(io, "    Optimization interval dx:")
    println(io, "        $(sp.dx[1]) < x < $(sp.dx[2])")

    if sp.p_space != [-1_000_000:1_000_000 for _ in 1:n_parameters(sp.sf)]
        println(io, "    Searched parameter space:")
        for (i, p) in enumerate(scaled_p_names(sp.sf))
            println(io, "        $p ∈ $(sp.p_space[i])")
        end
    else
        println(io, "    There was no scan of the parameter space.")
        println(io, "    The optimization started at the given starting points.")
        println(io, "    This is faster, but might not find the global minimum!")
    end
end

"""
    scaled_data(sp::ScalingProblem)
    scaled_data(sp::ScalingProblem, ps)

Return the scaled data of the `ScalingProblem` `sp`.

# Arguments
- `sp::ScalingProblem`: scaling problem
- `ps::Vector{Float64}=sp.optimal_ps`: parameters to scale the data with (manually)

# Returns
- `xs::Vector{Vector{Float64}}`: scaled x values
- `ys::Vector{Vector{Float64}}`: scaled y values
- `es::Vector{Vector{Float64}}`: scaled y error values
- `Ls::Vector{Float64}`: scaled system sizes
"""
function scaled_data(sp::ScalingProblem; kwargs...)
    return scaled_data(sp, sp.optimal_ps; kwargs...)
end

function scaled_data(sp::ScalingProblem, ps; kwargs...)
    scaled_data = [sp.sf.f(sp.data[i], ps...) for i in eachindex(sp.data)]
    xs = [scaled_data[i].xs for i in eachindex(scaled_data)]
    ys = [scaled_data[i].ys for i in eachindex(scaled_data)]
    es = [scaled_data[i].es for i in eachindex(scaled_data)]
    Ls = [scaled_data[i].L for i in eachindex(scaled_data)]

    if get(kwargs, :splines, false)
        if typeof(sp.quality) != MultipleSplines
            throw(ArgumentError("Splines are available only for the MultipleSplines quality function."))
        end

        y_splines = [
            Spline1D(xs[i], ys[i], k=3) for i in eachindex(scaled_data)
        ]
        e_splines = [
            Spline1D(xs[i], es[i], k=1) for i in eachindex(scaled_data)
        ]

        # set interval
        interval = Vector{Float64}(undef, 2)
        interval[1] = max(
            maximum(scaled_data[i].xs[1] for i in eachindex(scaled_data)),
            sp.dx[1]
        )
        interval[2] = min(
            minimum(scaled_data[i].xs[end] for i in eachindex(scaled_data)),
            sp.dx[2]
        )


        xvals = [
            range(interval[1], interval[2], length=sp.quality.N_steps)
            for _ in eachindex(scaled_data)
        ]
        yvals = [zeros(sp.quality.N_steps) for _ in eachindex(scaled_data)]
        evals = [zeros(sp.quality.N_steps) for _ in eachindex(scaled_data)]

        for l in eachindex(scaled_data)
            for (i, x) in enumerate(xvals[l])
                yvals[l][i] = y_splines[l](x)
                evals[l][i] = e_splines[l](x)
            end
        end
        return xs, ys, es, Ls, xvals, yvals, evals
    end

    return xs, ys, es, Ls
end
