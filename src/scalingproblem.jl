"""
    ScalingProblem(args...; kwargs...)

Create a scaling problem which solves on initialization.

# Arguments
- `xs`: x values of the data
- `ys`: y values of the data
- `es`: y error values of the data (optional)
- `Ls`: system sizes of the data
For more information on the arguments, see methods(Scaling.unzip_data).

# Keyword Arguments
- `sf::ScalingFunction=ScalingFunction(preset; kwargs...)`: scaling function
- `p_space::Vector{StepRangeLen}=[0.1:0.1:3.0 for _ in sf.p_names]`: parameter search space
- `dx::Vector{Float64}=[-Inf, Inf]`: optimization interval
- `quality::QualityFunction=Spline()`: quality function
- `quality_scan::QualityFunction=Spline(scan_mode=true)`: quality function for the
    parameter space scan (it is highly recommended to use the default function here)
- `verbose::Bool=false`: print information during optimization
- `starting_ps::Vector{Float64}`: If `starting_ps` are given, there will be no initial
    parameter space scan. Instead, the optimization will start at the given points. This is
    much faster, but might not find the global minimum.
- `error::Bool=true`: If `error=true`, the error analysis will be performed to give
    estimates of the uncertainties of the optimal parameters.

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
- `error::Bool`: perform error analysis
- `quality_scan::QualityFunction`: quality function for the parameter space scan
- `quality::QualityFunction`: quality function
- `skip_scan::Bool`: `true` if `starting_ps` are given by the user
- `starting_ps::Vector{Float64}`: starting points for the optimization
- `errors_defined::Bool`: `true` if errors are given by the user

# Examples
```julia
using Scaling
```

### rescale the x axis only
```julia
# (e.g. for ys == binder cumulant of Ising model)
ScalingProblem(xs, ys, Ls;
    sc=ScalingFunction(:x; p_names=["T_c", "nu"]),
    p_space=[0.1:0.1:3.0, 0.1:0.1:3.0],
    dx=[-1.0, 1.0],
)

# or the same but with errors
ScalingProblem(xs, ys, es, Ls;
    sc=ScalingFunction(:x; p_names=["T_c", "nu"]),
    p_space=[0.1:0.1:3.0, 0.1:0.1:3.0],
    dx=[-1.0, 1.0],
)
```

### rescale the x and y axis
```julia
# (e.g. for ys == magnetization of Ising model)
ScalingProblem(xs, ys, Ls;
    sc=ScalingFunction(:xy; p_names=["T_c", "nu", "beta"]),
    dx=[-1.0, 1.0],
)
```

### rescale the x and y axis, but y axis with negative exponent
```julia
# (e.g. for ys == susceptibility of Ising model)
ScalingProblem(xs, ys, Ls;
    sc=ScalingFunction(:xny; p_names=["T_c", "nu", "gamma"]),
    dx=[-1.0, 1.0],
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
    error::Bool
    quality_scan::QualityFunction
    quality::QualityFunction
    skip_scan::Bool
    starting_ps::Vector{Float64}
    errors_defined::Bool



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
        )
    end

    function ScalingProblem(args...; kwargs...)

        ##### unzip kwargs... #####
        sf = get(kwargs, :sf, ScalingFunction(; kwargs...))
        p_space = get(kwargs, :p_space, [0.1:0.1:3.0 for _ in 1:n_parameters(sf)])
        dx = get(kwargs, :dx, [-Inf, Inf])
        verbose = get(kwargs, :verbose, false)
        error = get(kwargs, :error, true)
        skip_scan = false
        starting_ps = get(kwargs, :starting_ps, zeros(n_parameters(sf)))
        if haskey(kwargs, :starting_ps)
            # if starting_ps are given, there will be no initial parameter space scan
            p_space = [-1_000_000:1_000_000 for _ in 1:n_parameters(sf)]
            skip_scan = true
        end
        quality_scan = get(kwargs, :quality_scan, Spline(scan_mode=true))
        quality = get(kwargs, :quality, Spline())

        @assert length(p_space) == n_parameters(sf)

        # unzip args...
        data = unzip_data(args...)

        # check if errors are defined
        errors_defined = false
        for d in data
            if any(d.es .!= 0.0)
                errors_defined = true
                break
            end
        end

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
        )

        # solve ScalingProblem
        solve!(sp)
        return sp
    end

end


function solve!(sp::ScalingProblem)

    # find global minimum
    optimize_parameters!(sp)

    # perform error analysis
    if sp.error
        error_analysis!(sp)
    end
end


function Base.show(io::IO, sp::ScalingProblem)
    println(io, "ScalingProblem")
    println(io, "    Scaling function:")
    println(io, "        x -> $(sp.sf.x_scale)")
    println(io, "        y -> $(sp.sf.y_scale)")
    if n_parameters(sp.sf) < length(sp.sf.p_names)
        println(io, "    Fixed parameters:")
        for (i, p) in enumerate(sp.sf.p_names)
            if sp.sf.fixed_ps[i] != Inf
                println(io, "        $p = $(sp.sf.fixed_ps[i])")
            end
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
        if typeof(sp.quality) != Spline
            throw(ArgumentError("Splines are available only for Spline quality function."))
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

        xvals = range(interval[1], interval[2], length=sp.quality.N_steps)
        yvals = zeros(sp.quality.N_steps, length(scaled_data))
        evals = zeros(sp.quality.N_steps, length(scaled_data))
        for l in eachindex(scaled_data)
            for (i, x) in enumerate(xvals)
                yvals[i, l] = y_splines[l](x)
                evals[i, l] = e_splines[l](x)
            end
        end
        return xs, ys, es, Ls, xvals, yvals, evals
    end

    return xs, ys, es, Ls
end
