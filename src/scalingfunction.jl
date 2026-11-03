"""
    ScalingFunction(preset::Symbol=:xy; p_names, relative=false, <name>=value, ...)
    ScalingFunction(f::Function; p_names, <name>=value, ..., x_scale, y_scale)
    ScalingFunction(p_names::Vector{String}; kwargs...)

The scaling form of a [`ScalingProblem`](@ref): how the data of one system size are
rescaled for given values of the scaling parameters.

# Presets
The presets rescale x around a critical point `p1` with the correlation length exponent
`p2`, and optionally y with a power of `L`:

| preset | x is rescaled to      | y is rescaled to | typical observable                     |
|:-------|:----------------------|:-----------------|:---------------------------------------|
| `:x`   | `(x - p1) * L^(1/p2)` | `y`              | dimensionless: Binder cumulant, `ξ/L`  |
| `:xy`  | `(x - p1) * L^(1/p2)` | `y * L^(p3/p2)`  | `y ∼ L^(-p3/p2)`: order parameter (β)  |
| `:xny` | `(x - p1) * L^(1/p2)` | `y * L^(-p3/p2)` | `y ∼ L^(p3/p2)`: susceptibility (γ)    |

Errors are rescaled like y. `:xy` is the default. Two `p_names` without a preset select
`:x`, three select `:xy`.

# Custom scaling functions
Any other scaling form is a function `f(d::ScalingCollapse.Data, p1, p2, ...)` that receives
the data `d` of one system size and the values of *all* parameters (also the fixed ones) and
returns the rescaled data as a new `ScalingCollapse.Data` (see [`ScalingCollapse.Data`](@ref)).
It must not modify `d`. Give the parameter names with `p_names` (or only their number with
`N_parameters`), and describe the scaling with `x_scale` and `y_scale` for display. To
make your own scaling form available as a preset symbol, add a method to
[`ScalingCollapse.handle_preset`](@ref).

# Keyword arguments
- `p_names::Vector{String}`: the names of the parameters (default: `"p1"`, `"p2"`, ...).
    They name the parameters in the output and select parameters to fix.
- `<name>=value`: fix the parameter `<name>` (one of `p_names`) to `value`. A fixed
    parameter is not fitted, has no entry in `p_space`, `starting_ps`, `optimal_ps` or
    `optimal_ps_error`, and is passed to the scaling function with its fixed value.
- `relative::Bool=false`: with `true`, the presets rescale x relative to the critical
    point, `(x - p1)/p1 * L^(1/p2)` (e.g. with the reduced temperature `(T - T_c)/T_c`);
    this was the default up to version 0.2. A window `dx` of the relative form equals
    `dx .* p1` in the absolute form, which also works for `p1 ≤ 0`. Only for presets.
- `x_scale::String`, `y_scale::String`: the rescaling of x and y as text, used to display
    the problem. The presets set them.
- `N_parameters::Int`: the number of parameters of a custom `f`, if `p_names` is not given.
- `f::Function`: replace the function of a preset, keeping its parameter names and
    `x_scale`/`y_scale` (the function must take as many parameters as the preset).

# Examples
```julia
ScalingFunction(:x; p_names=["T_c", "nu"])            # Binder cumulant
ScalingFunction(:x; p_names=["T_c", "nu"], nu=1.0)    # ν fixed to 1
ScalingFunction(:xy; p_names=["T_c", "nu", "beta"])   # magnetization
ScalingFunction(:xny; p_names=["T_c", "nu", "gamma"]) # susceptibility
ScalingFunction(:x; p_names=["T_c", "nu"], relative=true)

# custom: logarithmic shift of y, as for the specific heat of the 2D Ising model
function log_shift(d::ScalingCollapse.Data, T_c, nu, a)
    xs = @. (d.xs - T_c) * d.L^(1 / nu)
    ys = @. d.ys - a * log(d.L)
    return ScalingCollapse.Data(d.L, xs, ys, d.es)
end
ScalingFunction(log_shift; p_names=["T_c", "nu", "a"], nu=1.0,
    x_scale="(T - T_c) * L^(1/nu)", y_scale="C - a * log(L)")
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
    # for the presets, (L, free parameters...) -> (a, b, c) of _linear_scaling, so that the
    # quality functions can rescale in place; nothing for any other f
    coefficients::Union{Nothing,Function}

    function ScalingFunction(f::Function, p_names::Vector{String}, x_scale::String, y_scale::String; kwargs...)
        _check_kwargs(
            "ScalingFunction",
            keys(kwargs),
            (SCALINGFUNCTION_KWARGS..., Symbol.(p_names)...);
            hint="Parameters can be fixed by passing their name, one of $(p_names).",
        )
        fixed_idxs = Tuple(filter(i -> haskey(kwargs, Symbol(p_names[i])), eachindex(p_names)))
        free_idxs = Tuple(setdiff(eachindex(p_names), fixed_idxs))
        fixed_vals = Tuple(Float64(kwargs[Symbol(p_names[i])]) for i in fixed_idxs)
        coefficients = _linear_coefficients(f)
        return new(
            specialize_wrapper(f, fixed_idxs, free_idxs, fixed_vals),
            p_names,
            x_scale,
            y_scale,
            collect(fixed_vals),
            collect(fixed_idxs),
            collect(free_idxs),
            isnothing(coefficients) ? nothing :
            specialize_wrapper(coefficients, fixed_idxs, free_idxs, fixed_vals),
        )
    end
end

const SCALINGFUNCTION_KWARGS = (:f, :p_names, :x_scale, :y_scale, :N_parameters, :relative)

# throw an ArgumentError listing the valid keywords if any keyword is unknown
function _check_kwargs(caller, given, valid; hint="")
    unknown = [k for k in given if k ∉ valid]
    isempty(unknown) && return nothing
    throw(ArgumentError(
        "Unknown keyword argument$(length(unknown) > 1 ? "s" : "") " *
        "$(join(unknown, ", ")) for $caller. " *
        "Valid keyword arguments are: $(join(valid, ", ")). " * hint
    ))
end

function specialize_wrapper(f::Function, fixed_idxs::NTuple{N1,Int}, free_idxs::NTuple{N2,Int}, fixed_vals::NTuple{N1,Float64}) where {N1,N2}
    # position of every parameter in (fixed_vals..., args...), found here once rather than in
    # every call (Julia 1.10 does not fold the search into a constant and allocates)
    positions = ntuple(Val(N1 + N2)) do index
        index in fixed_idxs && return findfirst(==(index), fixed_idxs)
        return N1 + findfirst(==(index), free_idxs)
    end
    return (input, args::Vararg{Float64,N2}) -> begin
        values = (fixed_vals..., args...)
        parameter_tuple = ntuple(i -> values[positions[i]], Val(N1 + N2))
        f(input, parameter_tuple...)
    end
end

"""
    ScalingCollapse.handle_preset(::Val{preset}; kwargs...) -> (f, p_names, x_scale, y_scale)

Return the definition of the [`ScalingFunction`](@ref) preset `preset`. Add a method to
define a new preset, which `ScalingFunction(preset; kwargs...)` then accepts like the
built-in ones, including `p_names` and fixed parameters.

The method receives the keywords given to `ScalingFunction` and returns a tuple of
- `f`: the scaling function `f(d::ScalingCollapse.Data, p1, p2, ...) -> ScalingCollapse.Data`;
- `p_names::Vector{String}`: the default parameter names, usually overridden by
    `kwargs[:p_names]`;
- `x_scale::String`, `y_scale::String`: the rescaling of x and y as text, for display.

# Examples
```julia
using ScalingCollapse

# x rescaled with L^(1/nu), y shifted by a log(L)
function log_shift(d::ScalingCollapse.Data, x_c, nu, a)
    xs = @. (d.xs - x_c) * d.L^(1 / nu)
    ys = @. d.ys - a * log(d.L)
    return ScalingCollapse.Data(d.L, xs, ys, d.es)
end

function ScalingCollapse.handle_preset(::Val{:log_shift}; kwargs...)
    p_names = get(kwargs, :p_names, ["x_c", "nu", "a"])
    x_scale = "(x - \$(p_names[1])) * L^(1/\$(p_names[2]))"
    y_scale = "y - \$(p_names[3]) * log(L)"
    return log_shift, p_names, x_scale, y_scale
end

sf = ScalingFunction(:log_shift; p_names=["T_c", "nu", "a"], nu=1.0)
```
"""
function handle_preset end

function handle_preset(::Val{:x}; kwargs...)
    relative = _relative(kwargs)
    f = get(kwargs, :f, relative ? _power_scaling : _absolute_power_scaling)
    p_names = get(kwargs, :p_names, ["p1", "p2"])
    x_scale = get(kwargs, :x_scale) do
        _x_scale(p_names, relative)
    end
    y_scale = get(kwargs, :y_scale, "y")
    return f, p_names, x_scale, y_scale
end

# x and y scaling with (pos x), (pos y) exponent
function handle_preset(::Val{:xy}; kwargs...)
    relative = _relative(kwargs)
    f = get(kwargs, :f, relative ? _power_scaling : _absolute_power_scaling)
    p_names = get(kwargs, :p_names, ["p1", "p2", "p3"])
    x_scale = get(kwargs, :x_scale) do
        _x_scale(p_names, relative)
    end
    y_scale = get(kwargs, :y_scale) do
        "y * L^($(p_names[3])/$(p_names[2]))"
    end
    return f, p_names, x_scale, y_scale
end

# x and y scaling with (pos x), (neg y) exponent
function handle_preset(::Val{:xny}; kwargs...)
    relative = _relative(kwargs)
    f = get(kwargs, :f) do
        relative ? _negative_power_scaling : _absolute_negative_power_scaling
    end
    p_names = get(kwargs, :p_names, ["p1", "p2", "p3"])
    x_scale = get(kwargs, :x_scale) do
        _x_scale(p_names, relative)
    end
    y_scale = get(kwargs, :y_scale) do
        "y * L^(-$(p_names[3])/$(p_names[2]))"
    end
    return f, p_names, x_scale, y_scale
end

# the `relative` kwarg of the presets: x absolute (default) or relative to p1
function _relative(kwargs)
    haskey(kwargs, :relative) || return false
    relative = kwargs[:relative]
    relative isa Bool || throw(ArgumentError(
        "relative must be true or false, got $(repr(relative))."
    ))
    haskey(kwargs, :f) && throw(ArgumentError(
        "relative only selects the preset scaling function and cannot be combined with " *
        "a custom f."
    ))
    return relative
end

function _x_scale(p_names, relative)
    relative && return "(x - $(p_names[1]))/$(p_names[1]) * L^(1/$(p_names[2]))"
    return "(x - $(p_names[1])) * L^(1/$(p_names[2]))"
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
        throw(ArgumentError(
            "No preset supports $(length(p_names)) parameters (p_names = $(p_names)). " *
            "Pass a custom scaling function, e.g. ScalingFunction(f; p_names=...)."
        ))
    end
    ScalingFunction(preset, p_names=p_names; kwargs...)
end

function ScalingFunction(; kwargs...)
    haskey(kwargs, :p_names) && return ScalingFunction(kwargs[:p_names]; kwargs...)
    return ScalingFunction(:xy; kwargs...)
end

function ScalingFunction(f::Function; kwargs...)
    haskey(kwargs, :relative) && throw(ArgumentError(
        "relative only selects the preset scaling function and cannot be combined with " *
        "a custom f."
    ))
    p_names = if haskey(kwargs, :p_names)
        kwargs[:p_names]
    else
        if haskey(kwargs, :N_parameters)
            ["p$i" for i in 1:kwargs[:N_parameters]]
        else
            throw(ArgumentError("Please specify either :p_names or :N_parameters"))
        end
    end
    x_scale = get(kwargs, :x_scale, "custom scaling (set x_scale)")
    y_scale = get(kwargs, :y_scale, "custom scaling (set y_scale)")
    return ScalingFunction(f, p_names, x_scale, y_scale; kwargs...)
end

"""
    ScalingCollapse.n_parameters(sf::ScalingFunction)

The number of free (not fixed) parameters of `sf`.
"""
function n_parameters(sf::ScalingFunction)
    return length(sf.free_idxs)
end

"""
    ScalingCollapse.scaled_p_names(sf::ScalingFunction)

The names of the free (not fixed) parameters of `sf`, in the order of `p_space`,
`starting_ps`, `optimal_ps` and `optimal_ps_error`.
"""
function scaled_p_names(sf::ScalingFunction)
    return sf.p_names[sf.free_idxs]
end

"""
    ScalingCollapse.fixed_p_names(sf::ScalingFunction)

The names of the fixed parameters of `sf`.
"""
function fixed_p_names(sf::ScalingFunction)
    return sf.p_names[sf.fixed_idxs]
end

# coefficients (a, b) of x -> fma(x, a, b) = (x - p1)/p1 * L^(1/p2)
_relative_x_coefficients(L, p1, p2) = (L^(1 / p2) / p1, -L^(1 / p2))

# coefficients (a, b) of x -> fma(x, a, b) = (x - p1) * L^(1/p2)
function _absolute_x_coefficients(L, p1, p2)
    s = L^(1 / p2)
    return s, -p1 * s
end

# coefficients (a, b, c) of _linear_scaling of the presets for 2 and 3 parameters, x relative
# to p1 (relative=true) or absolute (the default); the 2-parameter forms leave y unchanged
_power_coefficients(L, p1, p2) = (_relative_x_coefficients(L, p1, p2)..., nothing)
_power_coefficients(L, p1, p2, p3) = (_relative_x_coefficients(L, p1, p2)..., L^(p3 / p2))
_absolute_power_coefficients(L, p1, p2) = (_absolute_x_coefficients(L, p1, p2)..., nothing)
function _absolute_power_coefficients(L, p1, p2, p3)
    return (_absolute_x_coefficients(L, p1, p2)..., L^(p3 / p2))
end
# negative y exponent (:xny)
_negative_power_coefficients(L, p1, p2, p3) = _power_coefficients(L, p1, p2, -p3)
function _absolute_negative_power_coefficients(L, p1, p2, p3)
    return _absolute_power_coefficients(L, p1, p2, -p3)
end

# the scaling functions of the presets (Vararg{Any,N} makes Julia specialize on N)
function _power_scaling(d::Data, ps::Vararg{Any,N}) where {N}
    return _linear_scaling(d, _power_coefficients(d.L, ps...)...)
end
function _absolute_power_scaling(d::Data, ps::Vararg{Any,N}) where {N}
    return _linear_scaling(d, _absolute_power_coefficients(d.L, ps...)...)
end
function _negative_power_scaling(d::Data, ps::Vararg{Any,N}) where {N}
    return _linear_scaling(d, _negative_power_coefficients(d.L, ps...)...)
end
function _absolute_negative_power_scaling(d::Data, ps::Vararg{Any,N}) where {N}
    return _linear_scaling(d, _absolute_negative_power_coefficients(d.L, ps...)...)
end

# the coefficient function of a preset scaling function f, nothing for any other f
_linear_coefficients(f) = nothing
_linear_coefficients(::typeof(_power_scaling)) = _power_coefficients
_linear_coefficients(::typeof(_absolute_power_scaling)) = _absolute_power_coefficients
_linear_coefficients(::typeof(_negative_power_scaling)) = _negative_power_coefficients
function _linear_coefficients(::typeof(_absolute_negative_power_scaling))
    return _absolute_negative_power_coefficients
end

# x -> fma(x, a, b), y -> y * c, e -> e * c (c > 0; y and e unchanged for c = nothing). For
# a < 0 the points are stored in reverse order, so that x stays sorted.
function _linear_scaling(d::Data, a, b, c=nothing)
    scaled = Data(d.L, similar(d.xs), similar(d.ys), similar(d.es), Val(:sorted))
    return _linear_scaling!(scaled, d, a, b, c)
end

# _linear_scaling writing into the arrays of scaled, which have the lengths of those of d
function _linear_scaling!(scaled::Data, d::Data, a, b, c)
    n = length(d.xs)
    for i in eachindex(d.xs)
        j = a < 0 ? n + 1 - i : i
        scaled.xs[j] = fma(d.xs[i], a, b)
        scaled.ys[j] = isnothing(c) ? d.ys[i] : d.ys[i] * c
        scaled.es[j] = isnothing(c) ? d.es[i] : d.es[i] * c
    end
    return scaled
end
