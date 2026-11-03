"""
    ScalingProblem(xs, ys, Ls; kwargs...)
    ScalingProblem(xs, ys, es, Ls; kwargs...)

Find the parameters of a scaling function that collapse the data `ys(xs)` of all system
sizes `Ls` onto one master curve.

The constructor solves the problem right away: it scans `p_space` for minima of the quality
function, refines the best one with Nelder–Mead and estimates the error of every free
parameter. The results are stored in the returned object:

- `sp.optimal_ps::Vector{Float64}`: the optimal values of the free parameters, in the order
    of `ScalingCollapse.scaled_p_names(sp.sf)` (fixed parameters are not included);
- `sp.optimal_ps_error::Vector{Float64}`: their errors (see `error` below; zeros with
    `error=false`);
- `sp.minimum::Float64`: the value of the quality function at `optimal_ps` (lower is
    better).

The input is kept in `sp.data`, one [`ScalingCollapse.Data`](@ref) per system size, and the
scaling function in `sp.sf`.

Use [`scaled_data`](@ref) to get the rescaled data, e.g. for a plot, and
[`parameter_errors`](@ref) to compute errors with another method after the fit.

# Data
`Ls` is a vector of integer system sizes. The x values `xs`, the y values `ys` and the
optional y errors `es` can each be given as
- a vector: the same values for every system size (typical for `xs`, e.g. a common
    temperature grid);
- a vector of vectors: one vector per system size, so different sizes can have different
    x values and different numbers of points;
- a matrix with one column or one row per system size. The orientation is inferred from
    `length(Ls)`; a square matrix is read column by column (with a warning).

`ys` can also contain `Measurement`s (Measurements.jl), whose uncertainties are then used as
`es`. If any error is nonzero, the quality functions use their error-weighted form.

The data of every system size are sorted by x. Within a system size the x values must be
distinct, and all values must be finite. At least two system sizes are needed.

# Keyword arguments
**Scaling form**
- `sf::ScalingFunction=ScalingFunction(:xy)`: the scaling function that rescales the data
    of every system size, see [`ScalingFunction`](@ref). Without `sf`, all keywords that
    are not listed here are passed on to `ScalingFunction`, so `p_names=["T_c", "nu"],
    nu=1.0` is short for `sf=ScalingFunction(:x; p_names=["T_c", "nu"], nu=1.0)` (two
    names select the preset `:x`, three names the preset `:xy`).
- `dx::Vector{Float64}=[-Inf, Inf]`: the window `[x_min, x_max]` of *rescaled* x values in
    which the collapse is measured. The quality function only compares the curves inside
    this window, intersected with the x range that all sizes cover. Choose it around the
    critical point, where the scaling form holds.

**Search**
- `p_space::Vector{<:AbstractRange}=[0.1:0.1:3.0, ...]`: one range per free parameter.
    The quality function is evaluated on every point of this grid, a Nelder–Mead fit
    starts from every grid point that is lower than its neighbours, and the best result
    starts the final fit. All fits are bounded by `p_space`, so the optimum must lie
    inside it. The scan costs `prod(length, p_space)` evaluations.
- `starting_ps`: skip the scan and start the final fit at this point (one value per free
    parameter), or at each of several points (a vector of such vectors, e.g.
    `[[2.25, 1.0], [2.3, 0.8]]`; the fit with the lowest finite quality is kept). The
    fits are then not bounded by `p_space`. This is faster than a scan, but only finds
    the minimum close to the start.
- `quality::QualityFunction=MultipleSplines()`: the quality function of the final fit and
    of the error analysis. See [`MultipleSplines`](@ref), [`Linear`](@ref),
    [`SingleSpline`](@ref) and [`Houdayer`](@ref).
- `quality_scan::QualityFunction=MultipleSplines(scan_mode=true)`: the quality function of
    the scan. The default ignores error bars, which makes the scan more robust; there is
    rarely a reason to change it.
- `optimizer::Optim.NelderMead`: the Nelder–Mead configuration of every fit, including the
    initial simplex and the algorithm parameters. The default initial simplex,
    `Optim.AffineSimplexer(a=0.025, b=0.1)`, moves every parameter `p` by `0.1p + 0.025`.
    For parameters on very different scales, use [`StepSimplexer`](@ref).
- `optimizer_options::Optim.Options=Optim.Options()`: tolerances, iteration and time
    limits, tracing and callbacks of every fit, including the fits of the error analysis.
- `verbose::Bool=false`: log the progress of the fit.
- `threads::Bool=true`: if Julia runs with several threads (e.g. `julia --threads=auto`),
    spread the scan, the fits from several starting points and the error analysis over
    them. The results are the same as with `threads=false`. A custom quality or scaling
    function is then called from several threads at once, so it must not modify shared
    state.

**Error analysis**
- `error::Union{Bool,Symbol}=true`: the method of the error analysis.
    - `true` or `:profile` (default): profile width. The parameter is displaced from its
        optimum, the other free parameters are re-fitted at every displacement, and the
        error is the larger of the two displacements (one on each side) at which the
        quality first reaches `error_threshold` times its minimum.
    - `:conditional`: as `:profile`, but the other parameters stay at their optimum. This
        ignores correlations between the parameters and gives smaller or equal errors. With
        one free parameter both methods agree.
    - `:bootstrap`: standard deviation of the parameters fitted to `n_bootstrap` synthetic
        data sets, in which every y value is redrawn from a Gaussian with its error bar as
        width. Needs error bars.
    - `false`: no error analysis.
- `error_threshold::Union{Real,Vector{<:Real}}=1.3`: the factor by which the quality may
    exceed its minimum within the error of a parameter (`:profile` and `:conditional`); one
    number, or one per free parameter. Must be larger than 1.
- `n_bootstrap::Integer=100`: the number of synthetic data sets of `error=:bootstrap` (at
    least 2).
- `rng::AbstractRNG=Random.default_rng()`: the random number generator of
    `error=:bootstrap`. Pass a seeded one, e.g. `Random.Xoshiro(42)`, for reproducible
    errors.

# Error analysis in detail
The profile and conditional errors are widths defined by a threshold on the quality, not
standard deviations with a known coverage: the threshold is not calibrated, and correlations
between data points (e.g. curves computed from the same Monte Carlo samples) are not taken
into account. The profile accounts for correlations between the parameters.

*Profile.* The re-fitted parameters stay inside `p_space`; the profiled parameter does not.
The crossing closest to the fit is used, so the error is a local width. Every inner fit
uses `optimizer` and `optimizer_options`. With an affine initial simplex (the default), it
is followed by a second fit whose initial steps are the conditional errors of the re-fitted
parameters, which resolves quality valleys much narrower than the affine steps; any other
initial simplex is used alone. If the profile finds a lower quality than `minimum` inside
`p_space`, the fit is re-optimized from there and the analysis restarts (with a warning, up
to 3 times), so `optimal_ps` and `minimum` can change. When the restarts run out, all errors
are `NaN`. A lower quality outside `p_space` sets the error of that parameter to `NaN` and
warns that `p_space` may be too narrow.

*Profile and conditional.* An error is `NaN` (with a warning) if the quality never reaches
the threshold, or becomes undefined (no overlap of the rescaled data within `dx`) before.

*Bootstrap.* Every synthetic data set is fitted with Nelder–Mead from `optimal_ps`, inside
`p_space`. `optimal_ps` and `minimum` stay those of the original data. The method assumes
Gaussian, independent errors. Fits without a finite quality are dropped (with a warning);
if more than 10% are dropped, all errors are `NaN`. The bootstrap costs about `n_bootstrap`
fits, and its errors scatter between random seeds by a relative amount of about
`1/sqrt(2(n_bootstrap - 1))`: about 7% for 100 samples, enough for one significant digit.
Use 500 or more for a reliable second digit.

# Examples
```julia
# Binder cumulant: rescale x only, U(T, L) = F((T - T_c) L^(1/ν))
sp = ScalingProblem(Ts, binder, Ls;
    sf=ScalingFunction(:x; p_names=["T_c", "nu"]),
    p_space=[2.0:0.02:2.5, 0.5:0.1:1.5],
    dx=[-2.0, 2.0],
)
sp.optimal_ps, sp.optimal_ps_error

# magnetization with error bars: M(T, L) = L^(-β/ν) F((T - T_c) L^(1/ν))
ScalingProblem(Ts, mag, mag_err, Ls; sf=ScalingFunction(:xy; p_names=["T_c", "nu", "beta"]))

# susceptibility, with ν fixed: χ(T, L) = L^(γ/ν) F((T - T_c) L^(1/ν))
ScalingProblem(Ts, chi, Ls; sf=ScalingFunction(:xny; p_names=["T_c", "nu", "gamma"], nu=1))

# no scan: start at a guess, with bootstrap errors
ScalingProblem(Ts, binder, binder_err, Ls;
    sf=ScalingFunction(:x; p_names=["T_c", "nu"]),
    starting_ps=[2.27, 1.0],
    error=:bootstrap, n_bootstrap=500, rng=Random.Xoshiro(1),
)

# configure Nelder–Mead
using ScalingCollapse: Optim
ScalingProblem(Ts, binder, Ls;
    sf=ScalingFunction(:x; p_names=["T_c", "nu"]),
    optimizer=Optim.NelderMead(initial_simplex=StepSimplexer([0.002, 0.05])),
    optimizer_options=Optim.Options(iterations=5_000, time_limit=60.0),
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
    threads::Bool



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
        threads,
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
            threads,
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
        starts = [zeros(n_ps)]
        if haskey(kwargs, :starting_ps)
            starts = _starting_points(kwargs[:starting_ps], sf)
            # if starting_ps are given, there will be no initial parameter space scan
            p_space = [-1_000_000:1_000_000 for _ in 1:n_ps]
            skip_scan = true
        end
        starting_ps = first(starts)
        quality_scan = get(kwargs, :quality_scan, MultipleSplines(scan_mode=true))
        quality = get(kwargs, :quality, MultipleSplines())
        error_threshold = _check_error_threshold(get(kwargs, :error_threshold, 1.3), sf)
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
        simplex = optimizer.initial_simplex
        (!(simplex isa StepSimplexer) || length(simplex.steps) == n_ps) || throw(
            ArgumentError(
                "The StepSimplexer of the optimizer has $(length(simplex.steps)) steps, " *
                "but the scaling function has $(n_ps) free parameters $(scaled_p_names(sf))."
            ),
        )
        optimizer_options = get(kwargs, :optimizer_options, Optim.Options())
        optimizer_options isa Optim.Options || throw(
            ArgumentError("optimizer_options must be an Optim.Options instance"),
        )
        n_bootstrap = _check_n_bootstrap(get(kwargs, :n_bootstrap, 100))
        rng = _check_rng(get(kwargs, :rng, Random.default_rng()))
        threads = get(kwargs, :threads, true)
        threads isa Bool || throw(ArgumentError(
            "threads must be true or false, got $(repr(threads))."
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
        _check_scaling_function(sf, first(data), skip_scan ? first(starts) : first.(p_space))

        # check if errors are defined
        errors_defined = false
        for d in data
            if any(d.es .!= 0.0)
                errors_defined = true
                break
            end
        end
        _check_bootstrap_data(error, errors_defined)

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
            threads,
        )

        # solve ScalingProblem
        solve!(sp, starts)
        return sp
    end

end


# normalize `starting_ps` (one vector or a vector of vectors) to a list of starting points
function _starting_points(starting_ps, sf)
    is_point(s) = s isa AbstractVector && all(v -> v isa Real, s)
    starts = if is_point(starting_ps) && !isempty(starting_ps)
        [starting_ps]
    elseif starting_ps isa AbstractVector && !isempty(starting_ps) &&
           all(is_point, starting_ps)
        starting_ps
    else
        throw(ArgumentError(
            "starting_ps must be a vector of parameter values or a non-empty vector of " *
            "such vectors (several starting points), got $(repr(starting_ps))."
        ))
    end
    n_ps = n_parameters(sf)
    for (i, s) in enumerate(starts)
        name = length(starts) == 1 ? "starting_ps" : "starting_ps[$i]"
        length(s) == n_ps || throw(ArgumentError(
            "$name has $(length(s)) entries, but the scaling function has $(n_ps) free " *
            "parameters $(scaled_p_names(sf))."
        ))
    end
    return [Vector{Float64}(s) for s in starts]
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
    :threads,
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

    # the interpolations need a minimum number of data points per system size (cubic
    # splines: 4)
    used = [q for q in (quality, skip_scan ? nothing : quality_scan) if q isa MultipleSplines]
    if !isempty(used)
        n_min = maximum(min_points, used)
        for d in data
            length(d.xs) >= n_min || throw(ArgumentError(
                "The MultipleSplines quality function needs at least $(n_min) data points " *
                "per system size, but L = $(d.L) has $(length(d.xs)). Use a different " *
                "quality / quality_scan function or provide starting_ps to skip the scan."
            ))
        end
    end
    return nothing
end

_workspace(sp::ScalingProblem) = Workspace(sp)

# whether to split the work into tasks (see _tmap)
_threaded(sp::ScalingProblem) = sp.threads && Threads.nthreads() > 1

# shallow copy of sp, e.g. to give it other data in a task
_copy(sp::ScalingProblem) = ScalingProblem((getfield(sp, k) for k in fieldnames(ScalingProblem))...)

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


function solve!(sp::ScalingProblem, starts=[sp.starting_ps])

    # find global minimum
    optimize_parameters!(sp, starts)

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

# error_threshold as one Float64 per free parameter (a number applies to all of them)
function _check_error_threshold(error_threshold, sf)
    n_ps = n_parameters(sf)
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
    return Float64.(error_threshold)
end

function _check_n_bootstrap(n_bootstrap)
    (n_bootstrap isa Integer && n_bootstrap >= 2) || throw(ArgumentError(
        "n_bootstrap must be an integer of at least 2, got $(repr(n_bootstrap))."
    ))
    return n_bootstrap
end

function _check_rng(rng)
    rng isa AbstractRNG || throw(ArgumentError(
        "rng must be an AbstractRNG (e.g. Random.Xoshiro(seed)), got $(typeof(rng))."
    ))
    return rng
end

function _check_bootstrap_data(error, errors_defined)
    (error === :bootstrap && !errors_defined) && throw(ArgumentError(
        "error=:bootstrap resamples the y values within their error bars, but the data " *
        "have no error bars. Pass y errors `es` or Measurement y values."
    ))
    return nothing
end


"""
    parameter_errors(sp::ScalingProblem; method, error_threshold, n_bootstrap, rng)

Compute the errors of the optimal parameters of the solved problem `sp`, without changing
`sp`. Use it to compare error methods, or to add errors to a fit made with `error=false`.
Returns one error per free parameter, in the order of `ScalingCollapse.scaled_p_names(sp.sf)`.

The keywords are those of the [`ScalingProblem`](@ref) constructor and default to the
settings stored in `sp`. `method` takes the values of the `error` keyword except `false`;
it defaults to `sp.error`, or to `:profile` if `sp` was fitted with `error=false`. The fits
of the analysis use the optimizer settings of `sp`.

The bootstrap (`method=:bootstrap`) advances `rng` (`sp.rng` by default). Pass a seeded one,
e.g. `rng=Random.Xoshiro(1)`, for reproducible errors; for the same `rng` state, the errors
equal those of a problem constructed with the same keywords.

If the profile (`method=:profile`) finds a lower quality than `sp.minimum`, `sp` is not at a
minimum, and `parameter_errors` throws an `ErrorException`. [`parameter_errors!`](@ref)
re-optimizes the fit in that case.

# Examples
```julia
using Random
sp = ScalingProblem(xs, ys, es, Ls; error=false)
σ_profile = parameter_errors(sp)
σ_conditional = parameter_errors(sp; method=:conditional)
σ_bootstrap = parameter_errors(sp; method=:bootstrap, n_bootstrap=500, rng=Xoshiro(1))
```
"""
function parameter_errors(sp::ScalingProblem; kwargs...)
    work = _error_problem(sp; kwargs...)
    error_analysis!(work; refit=false)
    return work.optimal_ps_error
end

"""
    parameter_errors!(sp::ScalingProblem; method, error_threshold, n_bootstrap, rng)

Like [`parameter_errors`](@ref), but store the errors in `sp.optimal_ps_error` and the
keywords in `sp` (`method` in `sp.error`), as if `sp` had been constructed with them.
Returns a copy of the errors.

If the profile finds a lower quality than `sp.minimum`, the fit is re-optimized from there
(with a warning) and the analysis restarts, as in the constructor, so `sp.optimal_ps` and
`sp.minimum` can change.
"""
function parameter_errors!(sp::ScalingProblem; kwargs...)
    work = _error_problem(sp; kwargs...)
    sp.error, sp.error_threshold = work.error, work.error_threshold
    sp.n_bootstrap, sp.rng = work.n_bootstrap, work.rng
    error_analysis!(sp)
    return copy(sp.optimal_ps_error)
end

# copy of sp with the given error settings (validated as in the constructor), whose fit
# results can change without changing sp
function _error_problem(sp::ScalingProblem;
    method=(sp.error === :none ? :profile : sp.error),
    error_threshold=sp.error_threshold,
    n_bootstrap=sp.n_bootstrap,
    rng=sp.rng,
)
    error = _error_method(method)
    error === :none && throw(ArgumentError(
        "method must be one of true, :profile, :conditional, :bootstrap, got $(repr(method))."
    ))
    _check_bootstrap_data(error, sp.errors_defined)
    return ScalingProblem(
        sp.data,
        sp.sf,
        sp.p_space,
        sp.dx,
        copy(sp.optimal_ps),
        copy(sp.optimal_ps_error),
        sp.minimum,
        sp.solved,
        sp.verbose,
        error,
        sp.quality_scan,
        sp.quality,
        sp.skip_scan,
        copy(sp.starting_ps),
        sp.errors_defined,
        _check_error_threshold(error_threshold, sp.sf),
        sp.optimizer,
        sp.optimizer_options,
        _check_n_bootstrap(n_bootstrap),
        _check_rng(rng),
        sp.threads,
    )
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
    scaled_data(sp::ScalingProblem[, ps]; splines=false) -> (xs, ys, es, Ls)

Rescale the data of `sp` with the scaling function `sp.sf` and the free parameters `ps`
(default: `sp.optimal_ps`). Returns one vector per system size, in the order of the sizes
given to `ScalingProblem`:

- `xs::Vector{Vector{Float64}}`: rescaled x values (sorted);
- `ys::Vector{Vector{Float64}}`: rescaled y values;
- `es::Vector{Vector{Float64}}`: rescaled y errors (zeros without error bars);
- `Ls::Vector{Int}`: the system sizes.

The result contains all data points, also those outside the window `sp.dx`.

With `splines=true` (only for the [`MultipleSplines`](@ref) quality function), three more
vectors are returned: `xvals`, `yvals` and `evals` hold, for every system size, the
`sp.quality.N_steps` evaluation points in the window that the quality function compares,
the interpolated y values there (as used by the quality function, before normalization),
and the linearly interpolated errors.

# Examples
```julia
using CairoMakie
xs, ys, es, Ls = scaled_data(sp)
fig = Figure()
ax = Axis(fig[1, 1])
for (x, y, L) in zip(xs, ys, Ls)
    scatterlines!(ax, x, y; label="L = \$L")
end
axislegend(ax)
fig
```
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
        # the same interpolation as the quality function
        yvals = [
            _interpolate(sp.quality.interp, xs[l], ys[l], collect(xvals[l]))
            for l in eachindex(scaled_data)
        ]
        evals = [zeros(sp.quality.N_steps) for _ in eachindex(scaled_data)]

        for l in eachindex(scaled_data)
            for (i, x) in enumerate(xvals[l])
                evals[l][i] = e_splines[l](x)
            end
        end
        return xs, ys, es, Ls, xvals, yvals, evals
    end

    return xs, ys, es, Ls
end
