```@meta
CurrentModule = ScalingCollapse
```

# 2D Ising model

The ferromagnetic Ising model on the square lattice is the textbook example of a continuous
phase transition, and its critical parameters are known exactly:

```math
T_c = \frac{2}{\ln(1 + \sqrt{2})} \approx 2.2692, \qquad \nu = 1, \qquad \beta = \frac{1}{8},
\qquad \gamma = \frac{7}{4}.
```

We determine them from Monte Carlo data for ``L = 8, 16, 32, 64`` at 81 temperatures between
1.5 and 3.5. The data ship with the repository of ScalingCollapse.jl (in `test/testdata`),
so the code below runs as it is.

## The data

```@example ising
using ScalingCollapse
using HDF5

file = joinpath(pkgdir(ScalingCollapse), "test", "testdata", "IsingMonteCarlo.h5")
Ts, Ls, binder, M, chi = h5open(file) do h
    read(h, "Ts"), read(h, "Ls"), read(h, "binder"), read(h, "M_abs_mean"), read(h, "susceptibility")
end
size(binder)
```

Every observable is a matrix with one column per system size and one row per temperature;
[`ScalingProblem`](@ref) recognizes the orientation from the number of sizes. To compare
raw and rescaled data throughout this page, we use a small plotting function. It shows the
rescaled data only within `window`, and marks the window `dx` of the fit with dashed lines.

```@example ising
using CairoMakie

function plot_collapse(sp; xlabel, ylabel, xlabel_scaled, ylabel_scaled, window=(-Inf, Inf))
    colors = cgrad(:Blues)[range(0.35, 1.0; length=length(sp.data))]
    fig = Figure(size=(780, 320))
    ax1 = Axis(fig[1, 1]; xlabel, ylabel, title="data")
    ax2 = Axis(fig[1, 2]; xlabel=xlabel_scaled, ylabel=ylabel_scaled, title="collapse")
    xs, ys, _, Ls = scaled_data(sp)
    for (i, d) in enumerate(sp.data)
        scatterlines!(ax1, d.xs, d.ys; color=colors[i], markersize=5, label="L = $(Ls[i])")
        inside = window[1] .<= xs[i] .<= window[2]
        scatter!(ax2, xs[i][inside], ys[i][inside]; color=colors[i], markersize=7)
    end
    vlines!(ax2, filter(isfinite, sp.dx); color=:gray, linestyle=:dash)
    Legend(fig[1, 3], ax1)
    return fig
end
nothing # hide
```

## Binder cumulant

The Binder cumulant ``U = 1 - \langle m^4 \rangle / (3 \langle m^2 \rangle^2)`` is
dimensionless: at ``T_c`` it does not depend on ``L``, and close to ``T_c`` it follows
``U(T, L) = \mathcal{F}\big((T - T_c) L^{1/\nu}\big)``. Only the x axis needs rescaling, so
we use the preset `:x`:

```@example ising
sp_binder = ScalingProblem(Ts, binder, Ls;
    sf=ScalingFunction(:x; p_names=["T_c", "nu"]),
    p_space=[2.0:0.02:2.6, 0.5:0.05:1.5],
    dx=[-2.0, 2.0],
)
```

Both parameters agree with the exact values within their errors.

```@example ising
plot_collapse(sp_binder; xlabel=L"T", ylabel=L"U",
    xlabel_scaled=L"(T - T_c)\,L^{1/\nu}", ylabel_scaled=L"U", window=(-6, 6))
```

### Fixing ν

Since ``\nu = 1`` is known, we can fix it and fit only ``T_c``, which reduces its error:

```@example ising
sp_binder_nu = ScalingProblem(Ts, binder, Ls;
    sf=ScalingFunction(:x; p_names=["T_c", "nu"], nu=1),
    p_space=[2.0:0.01:2.6],
    dx=[-2.0, 2.0],
)
```

## Magnetization

The order parameter vanishes at ``T_c`` as ``\langle |m| \rangle \sim L^{-\beta/\nu}``, so
``\langle |m| \rangle L^{\beta/\nu}`` is a function of ``(T - T_c) L^{1/\nu}`` only. This is
the preset `:xy`, which rescales y as ``y L^{p_3/p_2}``:

```@example ising
sp_mag = ScalingProblem(Ts, M, Ls;
    sf=ScalingFunction(:xy; p_names=["T_c", "nu", "beta"]),
    p_space=[2.1:0.02:2.4, 0.6:0.1:1.4, 0.05:0.025:0.3],
    dx=[-2.0, 2.0],
)
```

```@example ising
plot_collapse(sp_mag; xlabel=L"T", ylabel=L"\langle |m| \rangle",
    xlabel_scaled=L"(T - T_c)\,L^{1/\nu}", ylabel_scaled=L"\langle |m| \rangle\, L^{\beta/\nu}",
    window=(-6, 6))
```

## Susceptibility

The susceptibility diverges as ``\chi \sim L^{\gamma/\nu}``. The preset `:xny` rescales y as
``y L^{-p_3/p_2}``, so that ``\gamma`` comes out positive:

```@example ising
sp_chi = ScalingProblem(Ts, chi, Ls;
    sf=ScalingFunction(:xny; p_names=["T_c", "nu", "gamma"]),
    p_space=[2.1:0.02:2.4, 0.6:0.1:1.4, 1.2:0.1:2.2],
    dx=[-2.0, 2.0],
)
```

```@example ising
plot_collapse(sp_chi; xlabel=L"T", ylabel=L"\chi",
    xlabel_scaled=L"(T - T_c)\,L^{1/\nu}", ylabel_scaled=L"\chi\, L^{-\gamma/\nu}",
    window=(-6, 6))
```

## Summary

```@example ising
exact = Dict("T_c" => 2 / log(1 + sqrt(2)), "nu" => 1.0, "beta" => 1 / 8, "gamma" => 7 / 4)
for (name, sp) in ("Binder" => sp_binder, "Binder, ν = 1" => sp_binder_nu,
    "magnetization" => sp_mag, "susceptibility" => sp_chi)
    println(name)
    for (p, v, σ) in zip(ScalingCollapse.scaled_p_names(sp.sf), sp.optimal_ps, sp.optimal_ps_error)
        println("    ", rpad(p, 6), round(v; digits=4), " ± ", round(σ; digits=4),
            "   (exact ", round(exact[p]; digits=4), ", ", round(abs(v - exact[p]) / σ; digits=1), "σ)")
    end
end
```

All parameters lie within about two errors of the exact values. The largest deviation, ``\nu``
from the magnetization, comes from corrections to scaling, which are largest for the smallest
size, ``L = 8``. The data have no error bars, so the errors are the default profile widths
(see [Error analysis](@ref)).
