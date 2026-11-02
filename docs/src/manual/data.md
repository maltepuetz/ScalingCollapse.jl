```@meta
CurrentModule = ScalingCollapse
```

# Input data

A [`ScalingProblem`](@ref) takes the data as positional arguments,

```julia
ScalingProblem(xs, ys, Ls; kwargs...)      # without error bars
ScalingProblem(xs, ys, es, Ls; kwargs...)  # with error bars es of the y values
```

where `Ls` is a `Vector` of integer system sizes. `xs`, `ys` and `es` can each be given in
one of three forms, which can be mixed freely.

| form | meaning | typical use |
|:-----|:--------|:------------|
| `Vector{<:Real}` | the same values for every system size | a common grid of `xs` |
| `Vector{<:Vector}` | one vector per system size | sizes measured at different x values |
| `Matrix` | one column (or one row) per system size | data loaded from a file |

To illustrate, we build the same small data set in each form:

```@example data
using ScalingCollapse

Ls = [8, 16, 32]
xs = collect(range(0.0, 1.0; length=21))          # a common x grid
ys = [tanh.((xs .- 0.5) .* L^(1 / 1.5)) for L in Ls] # one vector per size (x_c = 0.5, ν = 1.5)
nothing # hide
```

A shared vector of x values with one vector of y values per size:

```@example data
sp1 = ScalingProblem(xs, ys, Ls; sf=ScalingFunction(:x), error=false)
sp1.optimal_ps
```

The same, with a matrix of y values with one column per system size:

```@example data
Y = reduce(hcat, ys)
size(Y)
```

```@example data
sp2 = ScalingProblem(xs, Y, Ls; sf=ScalingFunction(:x), error=false)
sp2.optimal_ps ≈ sp1.optimal_ps
```

The orientation of a matrix is inferred from `length(Ls)`: here `Y` has 3 columns, one per
size. A matrix with one row per size (`permutedims(Y)`) works as well. A square matrix is
ambiguous and read column by column, with a warning.

Different sizes may be measured at different x values, and may have different numbers of
points; then `xs` is a vector of vectors as well:

```@example data
xs3 = [collect(range(0.0, 1.0; length=n)) for n in (21, 31, 41)]
ys3 = [tanh.((x .- 0.5) .* L^(1 / 1.5)) for (x, L) in zip(xs3, Ls)]
sp3 = ScalingProblem(xs3, ys3, Ls; sf=ScalingFunction(:x), error=false)
sp3.optimal_ps
```

!!! warning "Ranges and system sizes"
    `xs` must be a `Vector` or a `Matrix`; `collect` ranges such as `0:0.1:1` first. `Ls`
    must be a `Vector` of integers; convert floating-point sizes with `Int.(Ls)`.

## Error bars

Pass the errors of the y values as a fourth positional argument `es`, in any of the forms
above. Alternatively, `ys` can contain `Measurement`s from
[Measurements.jl](https://github.com/JuliaPhysics/Measurements.jl), whose uncertainties are
then used as error bars:

```@example data
using Measurements

ys_m = [y .± 0.01 for y in ys]
sp4 = ScalingProblem(xs, ys_m, Ls; sf=ScalingFunction(:x), error=false)
sp4.data[1].es[1:3] # the stored errors of the smallest size
```

If any error is nonzero, the quality functions weight the deviations between the curves by
the error bars, and `error=:bootstrap` becomes available (see [Error analysis](@ref)).

## Requirements and preprocessing

- At least two system sizes.
- All values finite (no `NaN` or `Inf`).
- Within a system size, every x value at most once. Average repeated measurements first.
- Enough points per size for the quality function: 2 for the default
  `MultipleSplines()`, 4 for `MultipleSplines(interp=:cubic)`.

The data of every size are sorted by x, so the input order does not matter. The data are
stored as one [`ScalingCollapse.Data`](@ref) per system size in `sp.data`; the raw data of the
first size are `sp.data[1].xs`, `sp.data[1].ys` and `sp.data[1].es`.
