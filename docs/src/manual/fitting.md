```@meta
CurrentModule = ScalingCollapse
```

# Fitting

The constructor of [`ScalingProblem`](@ref) runs the whole fit. This page explains its steps
and the keywords that control them.

## How the fit works

1. **Scan.** The quality function `quality_scan` is evaluated on every point of the grid
   `p_space`. Every grid point that is lower than all its neighbours is a candidate.
2. **Local fits.** A Nelder–Mead fit starts from every candidate, bounded by `p_space`. The
   best result is the starting point of the final fit.
3. **Final fit.** A Nelder–Mead fit of the quality function `quality` from that starting
   point gives `optimal_ps` and `minimum`.
4. **Error analysis.** The errors of the parameters, see [Error analysis](@ref).

Pass `verbose=true` to log the progress of these steps. Throughout this page, we use
synthetic data of a dimensionless observable with ``x_c = 1`` and ``\nu = 0.8``:

```@example fit
using ScalingCollapse
using Random

F(z) = 1 / (1 + exp(z))
Ls = [10, 20, 40, 80]
xs = collect(range(0.8, 1.2; length=31))
rng = Xoshiro(1)
ys = [F.((xs .- 1.0) .* L^(1 / 0.8)) .+ 0.005 .* randn(rng, length(xs)) for L in Ls]

sp = ScalingProblem(xs, ys, Ls;
    sf=ScalingFunction(:x; p_names=["x_c", "nu"]),
    p_space=[0.9:0.02:1.1, 0.4:0.1:1.6],
    dx=[-2.0, 2.0],
)
```

## The parameter space

`p_space` holds one range per free parameter, in the order of
`ScalingCollapse.scaled_p_names(sp.sf)`. It serves two purposes:

- **The scan grid.** The grid should be fine enough to resolve the minimum you are looking
  for, but the scan costs `prod(length, p_space)` evaluations of the quality function. For
  three parameters, 20 to 30 points per parameter are usually plenty.
- **The bounds.** All fits stay inside `p_space`, so the optimum must lie inside it. An
  optimum at an edge of `p_space` is a sign that the range is too narrow.

The default, `0.1:0.1:3.0` for every parameter, is a generic first guess. A critical point
outside ``[0.1, 3]``, or a negative exponent, needs its own range.

## The window `dx`

`dx = [x_min, x_max]` limits the comparison of the curves to rescaled x values (after the
scaling function) in this window; the default `[-Inf, Inf]` compares them wherever all sizes
overlap. Restricting the window to the critical region reduces the bias from corrections to
scaling, at the price of fewer points and larger errors. Checking that the parameters do not
change by more than their errors when the window changes is a good test of a fit:

```@example fit
for dx in ([-1.0, 1.0], [-2.0, 2.0], [-4.0, 4.0])
    s = ScalingProblem(xs, ys, Ls;
        sf=ScalingFunction(:x; p_names=["x_c", "nu"]),
        p_space=[0.9:0.02:1.1, 0.4:0.1:1.6], dx=dx)
    println("dx = $dx: x_c = $(round(s.optimal_ps[1]; digits=4)) ± ",
        "$(round(s.optimal_ps_error[1]; digits=4)), nu = $(round(s.optimal_ps[2]; digits=3)) ± ",
        "$(round(s.optimal_ps_error[2]; digits=3))")
end
```

Remember that `dx` refers to the rescaled x axis, so its meaning depends on the scaling
function (see [the relative form](@ref "The relative form")).

## Starting points

If you have a good guess, `starting_ps` skips the scan and starts the final fit there:

```@example fit
sp = ScalingProblem(xs, ys, Ls;
    sf=ScalingFunction(:x; p_names=["x_c", "nu"]),
    starting_ps=[1.0, 1.0],
    dx=[-2.0, 2.0],
)
```

This is much faster, especially with many parameters, but Nelder–Mead only finds the minimum
close to the start, and the fit is no longer bounded by `p_space`. To be safer, pass several
starting points; the fit with the lowest quality is kept, and `sp.starting_ps` tells which
start it came from:

```@example fit
sp = ScalingProblem(xs, ys, Ls;
    sf=ScalingFunction(:x; p_names=["x_c", "nu"]),
    starting_ps=[[0.95, 0.6], [1.0, 1.0], [1.05, 1.4]],
    dx=[-2.0, 2.0],
)
sp.starting_ps
```

## The quality function

`quality` selects the measure of the collapse used by the final fit and the error analysis.
The default [`MultipleSplines`](@ref) is the right choice in most cases;
`MultipleSplines(interp=:cubic)` is more accurate for noise-free data with peaks. See
[Quality functions](@ref) for the alternatives. `quality_scan` is the measure used by the
scan; its default ignores error bars, which makes the scan more robust, and there is rarely a
reason to change it.

## Configuring Nelder–Mead

All fits use [Optim.jl](https://github.com/JuliaNLSolvers/Optim.jl)'s Nelder–Mead method.
`optimizer` sets the method (the initial simplex and the algorithm parameters),
`optimizer_options` the tolerances, iteration and time limits, tracing and callbacks:

```julia
using ScalingCollapse: Optim

ScalingProblem(xs, ys, Ls;
    sf=ScalingFunction(:x; p_names=["x_c", "nu"]),
    optimizer=Optim.NelderMead(initial_simplex=Optim.AffineSimplexer(a=0.01, b=0.05)),
    optimizer_options=Optim.Options(iterations=5_000, x_abstol=1e-10, time_limit=60.0),
)
```

The initial simplex sets the first steps of the fit. The default,
`Optim.AffineSimplexer(a=0.025, b=0.1)`, moves every parameter ``p`` by ``0.1 p + 0.025``.
This follows the value of a parameter, not its scale in the problem, which is a poor fit when
e.g. the critical point is ``0.5`` but the data only span ``10^{-3}`` around it. Then
[`StepSimplexer`](@ref) sets one step per free parameter, of the order of its expected error:

```julia
ScalingProblem(xs, ys, Ls;
    sf=ScalingFunction(:x; p_names=["p_c", "nu"]),
    p_space=[0.4995:0.00005:0.5005, 0.5:0.1:2.0],
    optimizer=Optim.NelderMead(initial_simplex=StepSimplexer([1e-4, 0.1])),
)
```

Vertices of the initial simplex outside `p_space`, or where the quality is not finite, are
moved towards the starting point automatically.

## Troubleshooting

**"The quality function is infinite for every parameter combination in p_space."**
For no grid point do the rescaled curves of all sizes overlap within `dx`. Check that the
x values of the sizes overlap near the critical point, widen `dx` or move `p_space`.

**"The parameter scan found no isolated minimum"**. The quality is flat along some
parameters, i.e. the data do not constrain them. Fix these parameters, or check `dx` and
`p_space`.

**"The Nelder-Mead optimization stopped after ... iterations without converging"**. Raise
the iteration limit with `optimizer_options=Optim.Options(iterations=...)`, or use a smaller
initial simplex (see above).

**The result lies at an edge of `p_space`.** The minimum is probably outside. Widen the range
of that parameter.

**The parameters depend strongly on `dx` or on the smallest sizes.** Corrections to scaling
are important. Use larger sizes, a narrower window, or a scaling function with a correction
term (see [Custom scaling functions](@ref)).

**The errors are `NaN`.** The error analysis warns about the reason. See
[Error analysis](@ref).
