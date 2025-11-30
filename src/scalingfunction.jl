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
    fixed_idxs::Vector{Int}
    free_idxs::Vector{Int}

    function ScalingFunction(f::Function, p_names::Vector{String}, x_scale::String, y_scale::String; kwargs...)
        fixed_idxs = Tuple(filter(i -> haskey(kwargs, Symbol(p_names[i])), eachindex(p_names)))
        free_idxs = Tuple(setdiff(eachindex(p_names), fixed_idxs))
        fixed_vals = Tuple(Float64(kwargs[Symbol(p_names[i])]) for i in fixed_idxs)
        return new(
            specialize_wrapper(f, fixed_idxs, free_idxs, fixed_vals),
            p_names,
            x_scale,
            y_scale,
            collect(fixed_vals),
            collect(fixed_idxs),
            collect(free_idxs)
        )
    end
end

function specialize_wrapper(f::Function, fixed_idxs::NTuple{N1,Int}, free_idxs::NTuple{N2,Int}, fixed_vals::NTuple{N1,Float64}) where {N1,N2}
    return (input, args::Vararg{Float64,N2}) -> begin
        parameter_tuple = ntuple(Val(N1 + N2)) do index
            if index in fixed_idxs
                idx = findfirst(==(index), fixed_idxs)
                return fixed_vals[idx]
            else
                idx = findfirst(==(index), free_idxs)
                return args[idx]
            end
        end
        f(input, parameter_tuple...)
    end
end

"""
    handle_preset(::Val{preset}; kwargs...)
Helper function to handle presets for ScalingFunction's. You can add your own by
adding a new method for `handle_preset`. The method should return the function `f`, the
parameter names `p_names`, the x scaling function `x_scale` and the y scaling function
`y_scale` as a tuple. The function f should take a `ScalingCollapse.Data` object and the
parameters of your scaling function and return a new `ScalingCollapse.Data` object which
contains the scaled data. Like so:
```julia
using ScalingCollapse
function _my_custom_scaling(d::ScalingCollapse.Data, p1, p2)

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

    return ScalingCollapse.Data(d.L, xs, ys, es)
end
```
`p_names` should be a vector of strings containing the names of the parameters used in the
scaling function. `x_scale` and `y_scale` are strings containing the scaling functions
for the x and y values respectively. These are just used for display purposes.

A new preset can then be added like so:
```julia
import ScalingCollapse: handle_preset
function handle_preset(::Val{:my_custom_preset}; kwargs...)
    f = _my_custom_scaling
    p_names = ["p1", "p2"]
    x_scale = "(x - p1)/p1 * L^(1/p2)"
    y_scale = "y"
    return f, p_names, x_scale, y_scale
end
```

We can now use the new preset like so:
```julia
sf = ScalingFunction(:my_custom_preset)
```
"""
function handle_preset end

function handle_preset(::Val{:x}; kwargs...)
    f = get(kwargs, :f, _power_scaling)
    p_names = get(kwargs, :p_names, ["p1", "p2"])
    x_scale = get(kwargs, :x_scale) do
        "(x - $(p_names[1]))/$(p_names[1]) * L^(1/$(p_names[2]))"
    end
    y_scale = get(kwargs, :y_scale, "y")
    return f, p_names, x_scale, y_scale
end

# x and y scaling with (pos x), (pos y) exponent
function handle_preset(::Val{:xy}; kwargs...)
    f = get(kwargs, :f, _power_scaling)
    p_names = get(kwargs, :p_names, ["p1", "p2", "p3"])
    x_scale = get(kwargs, :x_scale) do
        "(x - $(p_names[1]))/$(p_names[1]) * L^(1/$(p_names[2]))"
    end
    y_scale = get(kwargs, :y_scale) do
        "y * L^($(p_names[3])/$(p_names[2]))"
    end
    return f, p_names, x_scale, y_scale
end

# x and y scaling with (pos x), (neg y) exponent
function handle_preset(::Val{:xny}; kwargs...)
    f = get(kwargs, :f) do
        (d, p1, p2, p3) -> _power_scaling(d, p1, p2, -p3)
    end
    p_names = get(kwargs, :p_names, ["p1", "p2", "p3"])
    x_scale = get(kwargs, :x_scale) do
        "(x - $(p_names[1]))/$(p_names[1]) * L^(1/$(p_names[2]))"
    end
    y_scale = get(kwargs, :y_scale) do
        "y * L^(-$(p_names[3])/$(p_names[2]))"
    end
    return f, p_names, x_scale, y_scale
end

function handle_preset(::Val{T}; kwargs...) where {T}
    throw(ArgumentError("Unknown preset: $(T)"))
end

function ScalingFunction(preset::Symbol; kwargs...)
    f, p_names, x_scale, y_scale = handle_preset(Val(preset); kwargs...)
    ScalingFunction(f, p_names, x_scale, y_scale; kwargs...)
end

function ScalingFunction(p_names::Vector{String}; kwargs...)
    if length(p_names) == 2
        preset = :x
    elseif length(p_names) == 3
        preset = :xy
    else
        error("no preset supports %$(length(p_names)) parameters")
    end
    ScalingFunction(preset, p_names=p_names; kwargs...)
end

function ScalingFunction(; kwargs...)
    haskey(kwargs, :p_names) && return ScalingFunction(kwargs[:p_names]; kwargs...)
    return ScalingFunction(:xy; kwargs...)
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
    return ScalingFunction(f, p_names, x_scale, y_scale; kwargs...)
end

function n_parameters(sf::ScalingFunction)
    return length(sf.free_idxs)
end

function scaled_p_names(sf::ScalingFunction)
    return sf.p_names[sf.free_idxs]
end

function fixed_p_names(sf::ScalingFunction)
    return sf.p_names[sf.fixed_idxs]
end

# standard scaling function for 2 parameters
function _power_scaling(d::Data, p1, p2)
    xs = similar(d.xs)
    prefactor = d.L^(1 / p2) / p1
    shift = -d.L^(1 / p2)
    map!(x -> fma(x, prefactor, shift), xs, d.xs)
    return Data(d.L, xs, copy(d.ys), copy(d.es))
end

# standard scaling function for 3 parameters
function _power_scaling(d::Data, p1, p2, p3)
    xs = similar(d.xs)
    ys = similar(d.ys)
    es = similar(d.es)

    prefactor_x = d.L^(1 / p2) / p1
    prefactor_y = d.L^(p3 / p2)
    shift_x = -d.L^(1 / p2)
    map!(x -> fma(x, prefactor_x, shift_x), xs, d.xs)

    @assert length(ys) == length(es)
    for i in eachindex(ys)
        ys[i] = d.ys[i] * prefactor_y
        es[i] = d.es[i] * prefactor_y
    end

    return Data(d.L, xs, ys, es)
end
