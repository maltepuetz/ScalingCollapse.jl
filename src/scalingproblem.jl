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
- `sf::ScalingFunction=ScalingFunction(; kwargs...)`: scaling function
- `p_space::Vector{StepRangeLen}=[0.1:0.1:3.0 for _ in sf.p_names]`: search parameter space
- `dx::Vector{Float64}=[-Inf, Inf]`: optimization interval
- `verbose::Bool=false`: print information during optimization
- `starting_ps::Vector{Float64}`: If `starting_ps` are given, there will be no initial
    parameter space scan. Instead, the optimization will start at the given points. This is
    much faster, but might not find the global minimum.
- `error::Bool=false`: If `error=true`, the error analysis will be performed to give
    estimates of the uncertainties of the optimal parameters.

# Fields
- `data::Vector{Data}`: data to be scaled
- `sf::ScalingFunction`: scaling function
- `p_space::Vector{StepRangeLen}`: search parameter space
- `dx::Vector{Float64}`: optimization interval
- `optimal_ps::Vector{Float64}`: optimal parameters
- `optimal_ps_error::Vector{Float64}`: uncertainties of the optimal parameters
- `minimum::Float64`: minimum of the objective function

# Examples
```julia
using Scaling
```

### rescale the x axis only
    (e.g. for ys == binder cumulant of Ising model)
```julia
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
    (e.g. for ys == magnetization of Ising model)
```julia
ScalingProblem(xs, ys, Ls;
    sc=ScalingFunction(:xy; p_names=["T_c", "nu", "beta"]),
    dx=[-1.0, 1.0],
)
```

### rescale the x and y axis, but y axis with negative exponent
    (e.g. for ys == susceptibility of Ising model)
```julia
ScalingProblem(xs, ys, Ls;
    sc=ScalingFunction(:xny; p_names=["T_c", "nu", "gamma"]),
    dx=[-1.0, 1.0],
)
```

"""
struct ScalingProblem

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

    function ScalingProblem(args...; kwargs...)

        # unzip kwargs...
        sf = get(kwargs, :sf, ScalingFunction(; kwargs...))
        p_space = get(kwargs, :p_space, [0.1:0.1:3.0 for _ in 1:n_parameters(sf)])
        dx = get(kwargs, :dx, [-Inf, Inf])
        verbose = get(kwargs, :verbose, false)
        error = get(kwargs, :error, false)
        if haskey(kwargs, :starting_ps)
            # if starting_ps are given, there will be no initial parameter space scan
            p_space = [-1_000_000:1_000_000 for _ in 1:n_parameters(sf)]
        end

        @assert length(p_space) == n_parameters(sf)

        # unzip args...
        data = unzip_data(args...)

        # find global minimum
        optimal_ps, minimum = optimize_parameters(data, sf, p_space, verbose; kwargs...)

        # perform error analysis
        optimal_ps_error = zeros(size(optimal_ps))
        if error
            optimal_ps_error = error_analysis(data, sf, optimal_ps, minimum; kwargs...)
        end

        new(data, sf, p_space, dx, optimal_ps, optimal_ps_error, minimum)
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
    end
end


function scaled_data(sp::ScalingProblem)
    scaled_data = [sp.sf.f(sp.data[i], sp.optimal_ps...) for i in eachindex(sp.data)]
    xs = [scaled_data[i].xs for i in eachindex(scaled_data)]
    ys = [scaled_data[i].ys for i in eachindex(scaled_data)]
    es = [scaled_data[i].es for i in eachindex(scaled_data)]
    Ls = [scaled_data[i].L for i in eachindex(scaled_data)]

    return xs, ys, es, Ls
end