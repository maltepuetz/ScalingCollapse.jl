"""
    ScalingFunction(preset::Symbol; kwargs...)
    ScalingFunction(f::Function; kwargs...)

A `ScalingFunction` is used in a `ScalingProblem` to rescale the data. It is
defined by a function `f` and a set of parameters `p_names`. The function `f` is
called with the data `d` and the parameters `p1`, `p2`, ... and returns a new
`Data` object.

# Arguments
- `preset::Symbol`: a preset for the `ScalingFunction`. See `Presets` below.
- `f::Function`: the function `f` to be used in the `ScalingFunction`. Pass this
    argument to create a custom `ScalingFunction`. f should take a `Data` object
    and a set of parameters `p1`, `p2`, ... and return a new `Data` object:
    f(d::Data, p1, p2, ...) -> Data.

    Note that you might want to set the function `f` as a keyword argument to
    use the "fixing parameters" feature (see below).

# Presets
The `preset` argument can be used to create a `ScalingFunction` with a preset
function `f` and a preset set of parameter names `p_names`. The following
presets are available:
- `:x` rescale the x values of the data according to x -> (x - p1)/p1 * L^(1/p2)
- `:xy` rescale the x and y values of the data according to x -> (x - p1)/p1 * L^(1/p2)
    and y -> y * L^(p3/p2)
- `:xny` rescale the x and y values of the data according to x -> (x - p1)/p1 * L^(1/p2)
    and y -> y * L^(-p3/p2)

# Keyword arguments
- `p_names::Vector{String}`: the parameter names `p_names` to be used in the
    `ScalingFunction`. This kwarg can be used to overwrite the preset parameter
    names `p_names` (p1, p2, ...).
- `f::Function`: the function `f` to be used in the `ScalingFunction`. This kwarg can
    be used to overwrite the preset function `f`. If you pass f as a kwarg, instead of
    an argument, you can use the "fixing parameters" feature (see below). In that case
    the number of paramters of your function f should match the preset, i.e. if your
    function f takes 2 parameters, you should use the preset `:x` and if your function
    f takes 3 parameters, you should use the preset `:xy`.
- `x_scale::String`: the scaling function for the x values. This is just used to
    display the scaling function in plots and julias show function.
- `y_scale::String`: the scaling function for the y values. This is just used to
    display the scaling function in plots and julias show function.
- `p::Float64`: fix a parameter to a value. If you know the analytical value of a scaling
    parameter you might want to fix it to this value. This can speed up the optimization
    and give better results. Note, that `p` is either "p1", "p2", ... or the name you
    specified in `p_names`.

    Note that this feature is only available if you use one of the presets. You can still
    modify the scaling function by setting `f` as a kwarg.

# Examples
### Presets
```julia
# rescale the x axis only
ScalingFunction(:x; p_names=["T_c", "nu"])

# lets fix nu to 1.0
ScalingFunction(:x; p_names=["T_c", "nu"], nu=1.0)

# rescale the x and y axis
ScalingFunction(:xy; p_names=["T_c", "nu", "beta"])

# rescale the x and y axis with negative exponent on y
ScalingFunction(:xny; p_names=["T_c", "nu", "gamma"])
```

### Custom scaling function
```julia
# define the function that scales the data
function myfunction(d::ScalingCollapse.Data, p1, p2)

    #  initialize arrays for scaled data
    xs = zeros(length(d.xs))
    ys = zeros(length(d.ys))
    es = zeros(length(d.es))

    # scale data according to p1 and p2
    for (i, x, y, e) in zip(eachindex(d.xs), d.xs, d.ys, d.es)
        xs[i] = (x - p1) / p1 * log(d.L)
        ys[i] = y * d.L^(p2)
        es[i] = e * d.L^(p2)
    end

    # create new Data object with scaled data
    return ScalingCollapse.Data(d.L, xs, ys, es)
end

ScalingFunction(myfunction; p_names=["myp1", "myp2"])

# lets say we want to fix the parameter "myp1" to 1.0
# we use the preset :x (because it also has 2 parameters) and overwrite the function f
ScalingFunction(:x; f=myfunction, p_names=["myp1", "myp2"], myp1=1.0)
```
"""
struct ScalingFunction

    f::Function
    p_names::Vector{String}
    x_scale::String
    y_scale::String
    fixed_ps::Vector{Float64}

    function ScalingFunction(preset::Symbol; kwargs...)

        #=
        Presets:
            :x    -> x scaling with positive exponent
            :xy   -> x and y scaling with positive exponents
            :xny  -> x and y scaling with (pos x), (neg y) exponent
        =#

        if preset == :x
            f = get(kwargs, :f, _power_scaling)
            p_names = get(kwargs, :p_names, ["p1", "p2"])
            x_scale = get(
                kwargs,
                :x_scale,
                "(x - $(p_names[1]))/$(p_names[1]) * L^(1/$(p_names[2]))"
            )
            y_scale = get(kwargs, :y_scale, "y")

            # check whether a parameter was provided as kwarg
            fixed_ps = [Inf for _ in 1:length(p_names)]
            for (i, p) in enumerate(p_names)
                if haskey(kwargs, Symbol(p))
                    fixed_ps[i] = kwargs[Symbol(p)]
                end
            end

            return new(
                (d, args...) -> f(
                    d,
                    (_i = 0; fixed_ps[1] == Inf) ? (_i += 1; args[_i]) : fixed_ps[1],
                    fixed_ps[2] == Inf ? (_i += 1; args[_i]) : fixed_ps[2]
                ),
                p_names,
                x_scale,
                y_scale,
                fixed_ps
            )

        elseif preset == :xy
            f = get(kwargs, :f, _power_scaling)
            p_names = get(kwargs, :p_names, ["p1", "p2", "p3"])
            x_scale = get(
                kwargs,
                :x_scale,
                "(x - $(p_names[1]))/$(p_names[1]) * L^(1/$(p_names[2]))")
            y_scale = get(
                kwargs,
                :y_scale,
                "y * L^($(p_names[3])/$(p_names[2]))"
            )

            # check whether a parameter was provided as kwarg
            fixed_ps = [Inf for _ in 1:length(p_names)]
            for (i, p) in enumerate(p_names)
                if haskey(kwargs, Symbol(p))
                    fixed_ps[i] = kwargs[Symbol(p)]
                end
            end

            return new(
                (d, args...) -> f(
                    d,
                    (_i = 0; fixed_ps[1] == Inf) ? (_i += 1; args[_i]) : fixed_ps[1],
                    fixed_ps[2] == Inf ? (_i += 1; args[_i]) : fixed_ps[2],
                    fixed_ps[3] == Inf ? (_i += 1; args[_i]) : fixed_ps[3]
                ),
                p_names,
                x_scale,
                y_scale,
                fixed_ps
            )

            #return new(f, p_names, x_scale, y_scale, zeros(3)) # TODO fixed_ps

        elseif preset == :xny
            # put negative exponent on y
            f = get(
                kwargs,
                :f,
                (d::Data, p1, p2, p3) -> _power_scaling(d::Data, p1, p2, -p3)
            )
            p_names = get(kwargs, :p_names, ["p1", "p2", "p3"])
            x_scale = get(
                kwargs,
                :x_scale,
                "(x - $(p_names[1]))/$(p_names[1]) * L^(1/$(p_names[2]))"
            )
            y_scale = get(
                kwargs,
                :y_scale,
                "y * L^(-$(p_names[3])/$(p_names[2]))"
            )

            # check whether a parameter was provided as kwarg
            fixed_ps = [Inf for _ in 1:length(p_names)]
            for (i, p) in enumerate(p_names)
                if haskey(kwargs, Symbol(p))
                    fixed_ps[i] = kwargs[Symbol(p)]
                end
            end

            return new(
                (d, args...) -> f(
                    d,
                    (_i = 0; fixed_ps[1] == Inf) ? (_i += 1; args[_i]) : fixed_ps[1],
                    fixed_ps[2] == Inf ? (_i += 1; args[_i]) : fixed_ps[2],
                    fixed_ps[3] == Inf ? (_i += 1; args[_i]) : fixed_ps[3]
                ),
                p_names,
                x_scale,
                y_scale,
                fixed_ps
            )
        else
            throw(ArgumentError("Unknown preset: $preset"))
        end
    end

    function ScalingFunction(f::Function; kwargs...)
        p_names = if haskey(kwargs, :p_names)
            kwargs[:p_names]
        else
            if haskey(kwargs, :N_parameters)
                ["p$i" for i in 1:kwargs[:N_parameters]]
            else
                throw(ArgumentError("Please specify either :p_names or :N_parameters"))
            end
        end
        x_scale = get(kwargs, :x_scale, "tmp")
        y_scale = get(kwargs, :y_scale, "tmp")
        return new(f, p_names, x_scale, y_scale, [Inf for _ in p_names]) # TODO fixed_ps
    end

    function ScalingFunction(p_names::Vector{String}; kwargs...)
        length(p_names) == 2 && return ScalingFunction(:x, p_names=p_names; kwargs...)
        length(p_names) == 3 && return ScalingFunction(:xy, p_names=p_names; kwargs...)
    end

    function ScalingFunction(; kwargs...)
        haskey(kwargs, :p_names) && return ScalingFunction(kwargs[:p_names]; kwargs...)
        return ScalingFunction(:xy; kwargs...)
    end
end

function n_parameters(sf::ScalingFunction)
    n = length(sf.fixed_ps)
    for p in sf.fixed_ps
        if p != Inf
            n -= 1
        end
    end
    return n
end

function scaled_p_names(sf::ScalingFunction)
    res = sf.p_names |> copy
    for i in length(sf.fixed_ps):-1:1
        if sf.fixed_ps[i] != Inf
            deleteat!(res, i)
        end
    end
    return res
end

function fixed_p_names(sf::ScalingFunction)
    res = sf.p_names |> copy
    for i in length(sf.fixed_ps):-1:1
        if sf.fixed_ps[i] == Inf
            deleteat!(res, i)
        end
    end
    return res
end

# standard scaling function for 2 parameters
function _power_scaling(d::Data, p1, p2)

    xs = zeros(length(d.xs))
    ys = zeros(length(d.ys))
    es = zeros(length(d.es))

    for (i, x) in enumerate(d.xs)
        xs[i] = (x - p1) / p1 * d.L^(1 / p2)
    end
    for (i, y) in enumerate(d.ys)
        ys[i] = y
    end
    for (i, e) in enumerate(d.es)
        es[i] = e
    end

    return Data(d.L, xs, ys, es)
end

# standard scaling function for 3 parameters
function _power_scaling(d::Data, p1, p2, p3)

    xs = zeros(length(d.xs))
    ys = zeros(length(d.ys))
    es = zeros(length(d.es))

    for (i, x) in enumerate(d.xs)
        xs[i] = (x - p1) / p1 * d.L^(1 / p2)
    end
    for (i, y) in enumerate(d.ys)
        ys[i] = y * d.L^(p3 / p2)
    end
    for (i, e) in enumerate(d.es)
        es[i] = e * d.L^(p3 / p2)
    end

    return Data(d.L, xs, ys, es)
end
