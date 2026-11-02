```@meta
CurrentModule = ScalingCollapse
```

# Uniform grids and quality landscapes

A fit is only as good as the quality function it minimizes. This example uses
[`residuals`](@ref) to look at the quality landscape of all quality functions, and shows why
the default [`MultipleSplines`](@ref) is smooth where the others are not.

## The data

We use the spanning probability of bond percolation on the square lattice (``p_c = 1/2``,
``\nu = 4/3``), computed with the Newman–Ziff algorithm for ``L = 16`` to ``256``, with error
bars. The data are in `test/testdata/percolation` of the repository. The probabilities ``p``
were chosen on a grid that is *uniform in the scaling variable*
``(p - 1/2)/(1/2) \cdot L^{0.7}``, a common choice to place the same number of points in the
critical region of every size.

```@example land
using ScalingCollapse
using HDF5

file = joinpath(pkgdir(ScalingCollapse), "test", "testdata", "percolation", "bond_square.h5")
Ls, xs, ys, es = h5open(file) do h
    Ls = read(h, "Ls")
    curves(v) = [read(h, "spanning/L$L/$v") for L in Ls]  # one vector per size
    Ls, curves("x"), curves("y"), curves("e")
end
nothing # hide
```

## Fits with every quality function

```@example land
sf = ScalingFunction(:x; p_names=["p_c", "nu"])
setup = (; sf, p_space=[0.48:0.002:0.52, 0.8:0.1:2.0], dx=[-0.5, 0.5])
qualities = ["MultipleSplines" => MultipleSplines(), "Linear" => Linear(),
    "SingleSpline" => SingleSpline(), "Houdayer" => Houdayer()]
fits = [name => ScalingProblem(xs, ys, es, Ls; setup..., quality) for (name, quality) in qualities]
for (name, sp) in fits
    println(rpad(name, 16), "p_c = ", round(sp.optimal_ps[1]; digits=5), " ± ",
        round(sp.optimal_ps_error[1]; digits=5), ",  nu = ", round(sp.optimal_ps[2]; digits=3),
        " ± ", round(sp.optimal_ps_error[2]; digits=3))
end
```

All four fits agree with each other, and with the exact values to within about two errors.
Their landscapes do not.

## The landscapes

[`residuals`](@ref) evaluates the quality function of a problem on a parameter grid. We plot
``\log_{10}(S/S_\text{min})`` on the same grid for every quality function:

```@example land
using CairoMakie

grid = [0.497:0.0001:0.503, 1.0:0.005:1.7]
fig = Figure(size=(780, 620))
for (k, (name, sp)) in enumerate(fits)
    ps, S = residuals(sp; p_space=grid)
    ax = Axis(fig[fldmod1(k, 2)...]; title=name, xlabel=L"p_c", ylabel=L"\nu")
    heatmap!(ax, ps[1], ps[2], log10.(S ./ minimum(S));
        colormap=Reverse(:Blues), colorrange=(0, 1.5))
    scatter!(ax, [0.5], [4 / 3]; color=:white, strokecolor=:black, strokewidth=1.5,
        markersize=10)
end
Colorbar(fig[1:2, 3]; colormap=Reverse(:Blues), colorrange=(0, 1.5),
    label=L"\log_{10}(S / S_\mathrm{min})")
fig
```

The white dot marks the exact values. `MultipleSplines` has a single smooth valley close to
them.
`SingleSpline` and `Houdayer` are rippled, and `Linear` has a sharp cross-shaped structure at
``p_c = 1/2``, ``\nu \approx 1.43``.

A cut along ``\nu`` at the exact ``p_c = 1/2`` shows what happens there:

```@example land
fig = Figure(size=(780, 360))
ax = Axis(fig[1, 1]; xlabel=L"\nu", ylabel=L"S / S_\mathrm{min}", yscale=log10,
    yticks=[1, 2, 5, 10, 20, 50], title=L"quality at $p_c = 1/2$")
for (name, sp) in fits
    ps, S = residuals(sp; dims=[1, 2], p_space=[0.5:0.5, 1.0:0.002:1.8])
    lines!(ax, ps[2], S[1, :] ./ minimum(S); linewidth=2, label=name)
end
vlines!(ax, [4 / 3, 1 / 0.7]; color=:gray, linestyle=:dash)
Legend(fig[1, 2], ax)
fig
```

The dashed lines mark ``\nu = 4/3`` and ``\nu = 1/0.7``. At ``p_c = 1/2`` and
``\nu = 1/0.7``, the rescaled x values of all sizes coincide, because the data grid is uniform
in ``(p - 1/2) L^{0.7}``. `Linear` then compares points of different sizes at the same x,
without any interpolation error, and its quality drops to a false minimum, lower than the one
near the true ``\nu``. The quality functions that compare individual points (`Linear`,
`SingleSpline`, `Houdayer`) also jump whenever points of different sizes pass each other,
which produces the ripples.

`MultipleSplines` compares smooth interpolations of every size on a fixed set of points, so it
does not care how the grids of the sizes line up. Its landscape is smooth, which makes both
the fit and the [error analysis](@ref "Error analysis") reliable.

!!! tip
    Whatever the quality function, plotting the landscape around a fit with
    [`residuals`](@ref) is a quick check that the minimum is real and not an artifact of the
    data grid.
