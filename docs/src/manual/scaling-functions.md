```@meta
CurrentModule = ScalingCollapse
```

# Scaling functions

A [`ScalingFunction`](@ref) defines how the data of one system size are rescaled for given
values of the parameters. The fit searches the parameters that make the rescaled curves of
all sizes collapse.

## Presets

Three presets cover the usual scaling forms (see [Finite-size scaling](@ref) for which
observable needs which):

| preset | x is rescaled to | y is rescaled to |
|:-------|:-----------------|:-----------------|
| `:x`   | ``(x - p_1)\, L^{1/p_2}`` | ``y`` |
| `:xy` (default) | ``(x - p_1)\, L^{1/p_2}`` | ``y\, L^{p_3/p_2}`` |
| `:xny` | ``(x - p_1)\, L^{1/p_2}`` | ``y\, L^{-p_3/p_2}`` |

Error bars are rescaled like the y values. `p_names` gives the parameters meaningful names,
which appear in the output of the problem:

```@example sf
using ScalingCollapse

sf = ScalingFunction(:xny; p_names=["T_c", "nu", "gamma"])
sf.x_scale, sf.y_scale
```

As a shortcut, the keywords of the scaling function can be passed directly to
[`ScalingProblem`](@ref) when no `sf` is given. Two names select the preset `:x`, three
names `:xy`:

```julia
ScalingProblem(xs, ys, Ls; p_names=["T_c", "nu"])  # same as sf=ScalingFunction(:x; p_names=["T_c", "nu"])
```

## Fixing parameters

Pass a parameter's name and value to fix it:

```@example sf
sf = ScalingFunction(:xy; p_names=["T_c", "nu", "beta"], nu=1.0)
ScalingCollapse.scaled_p_names(sf)
```

A fixed parameter is not fitted. `p_space`, `starting_ps`, `optimal_ps` and
`optimal_ps_error` only contain the free parameters, in the order given by
`ScalingCollapse.scaled_p_names(sf)`. Fixing parameters works for presets and custom scaling
functions alike, e.g. `nu=1.0` for the 2D Ising model, see the
[Ising example](@ref "Fixing ν").

## The relative form

With `relative=true`, the presets rescale x relative to the critical point, as
``(x - p_1)/p_1 \cdot L^{1/p_2}``, e.g. with the reduced temperature ``(T - T_c)/T_c``:

```@example sf
ScalingFunction(:x; p_names=["T_c", "nu"], relative=true).x_scale
```

This was the default up to version 0.2. Both forms give the same parameters, but a window
`dx` means different things: the window `dx` of the relative form is the window `dx .* p1` of
the absolute form. To reproduce an old fit, pass `relative=true` or multiply `dx` by the
critical point. The absolute form also works for ``p_1 \le 0``.

## Custom scaling functions

Any other scaling form is a function `f(d, p1, p2, ...)` that takes the data
`d::ScalingCollapse.Data` of one system size and the values of **all** parameters (fixed ones
included), and returns the rescaled data as a new [`ScalingCollapse.Data`](@ref):

- `d.L` is the system size, `d.xs`, `d.ys`, `d.es` are the x values, y values and errors;
- return `ScalingCollapse.Data(d.L, new_xs, new_ys, new_es)` and do not modify `d`;
- rescale the errors like the y values.

As an example, consider strong corrections to scaling: the order parameter of the synthetic
data below follows ``Y = L^{-\beta/\nu} (1 + b L^{-\omega}) \mathcal{F}\big((x - x_c)
L^{1/\nu}\big)``, with ``b = 2`` and ``\omega = 1``, i.e. a 25% correction for ``L = 8``:

```@example sf
using Random

x_c, ν, β, b, ω = 0.5, 1.2, 0.3, 2.0, 1.0
F(z) = 1 / (1 + exp(2z)) + 0.2
Ls = [8, 16, 32, 64]
xs = collect(range(0.2, 0.8; length=41))
rng = Xoshiro(7)
ys = [L^(-β / ν) * (1 + b * L^(-ω)) .* F.((xs .- x_c) .* L^(1 / ν)) .*
      (1 .+ 0.002 .* randn(rng, length(xs))) for L in Ls]
nothing # hide
```

The preset `:xy` ignores the correction, and ``\beta`` comes out several standard deviations
too large:

```@example sf
sp_preset = ScalingProblem(xs, ys, Ls;
    sf=ScalingFunction(:xy; p_names=["x_c", "nu", "beta"]),
    p_space=[0.4:0.01:0.6, 0.5:0.1:2.0, 0.0:0.05:0.6],
    dx=[-3.0, 3.0],
)
```

A custom scaling function divides the correction out. We fit its amplitude `b` and fix the
correction exponent `omega`; `x_scale` and `y_scale` describe the scaling for the output:

```@example sf
function corrected(d::ScalingCollapse.Data, x_c, nu, beta, b, omega)
    xs = @. (d.xs - x_c) * d.L^(1 / nu)
    c = d.L^(beta / nu) / (1 + b * d.L^(-omega))
    return ScalingCollapse.Data(d.L, xs, d.ys .* c, d.es .* c)
end

sf = ScalingFunction(corrected;
    p_names=["x_c", "nu", "beta", "b", "omega"],
    omega=1.0,
    x_scale="(x - x_c) * L^(1/nu)",
    y_scale="y * L^(beta/nu) / (1 + b * L^(-omega))",
)
sp_custom = ScalingProblem(xs, ys, Ls;
    sf,
    p_space=[0.4:0.02:0.6, 0.5:0.1:2.0, 0.0:0.05:0.6, 0.0:0.5:4.0],
    dx=[-3.0, 3.0],
)
```

All four free parameters agree with the true values. Every additional parameter makes the scan
more expensive (here ``11 \times 16 \times 13 \times 9`` grid points); with a good guess,
`starting_ps` skips it, see [Fitting](@ref).

!!! tip "Type stability"
    The scaling function is called many times during a fit. Broadcasting as above, or loops
    over preallocated vectors, keeps it fast. The parameters are passed as `Float64`.

## Custom presets

To reuse a custom scaling form under a name, add a method to
[`ScalingCollapse.handle_preset`](@ref). It returns the function, the default parameter
names and the display strings:

```@example sf
function ScalingCollapse.handle_preset(::Val{:corrected}; kwargs...)
    p_names = get(kwargs, :p_names, ["x_c", "nu", "beta", "b", "omega"])
    x_scale = "(x - $(p_names[1])) * L^(1/$(p_names[2]))"
    y_scale = "y * L^($(p_names[3])/$(p_names[2])) / (1 + $(p_names[4]) * L^(-$(p_names[5])))"
    return corrected, p_names, x_scale, y_scale
end

sf = ScalingFunction(:corrected; omega=1.0)
sf.y_scale
```

The new preset accepts `p_names` and fixed parameters like the built-in ones.
