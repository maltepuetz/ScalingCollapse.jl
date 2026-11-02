```@meta
CurrentModule = ScalingCollapse
```

# Specific heat: a custom scaling function

The specific heat of the 2D Ising model does not diverge with a power law: its critical
exponent is ``\alpha = 0``, and the divergence is logarithmic. Near ``T_c``, the specific heat
per spin of a system of size ``L`` follows

```math
C(T, L) = A \ln L + \mathcal{F}\big((T - T_c)\, L^{1/\nu}\big),
```

with the exact amplitude ``A = \frac{2}{\pi} \ln^2(1 + \sqrt{2}) \approx 0.4945``
(Ferdinand and Fisher, 1969). None of the presets shifts y by a logarithm, so we write a
custom scaling function, using the Monte Carlo data of the [2D Ising example](@ref "2D Ising model").

```@example cv
using ScalingCollapse
using HDF5

file = joinpath(pkgdir(ScalingCollapse), "test", "testdata", "IsingMonteCarlo.h5")
Ts, Ls, C = h5open(file) do h
    read(h, "Ts"), read(h, "Ls"), read(h, "specific_heat")
end
nothing # hide
```

## The scaling function

A scaling function receives the data `d` of one system size and the values of all parameters,
and returns the rescaled data as a new `ScalingCollapse.Data`:

```@example cv
function log_shift(d::ScalingCollapse.Data, T_c, nu, A)
    xs = @. (d.xs - T_c) * d.L^(1 / nu)
    ys = @. d.ys - A * log(d.L)
    return ScalingCollapse.Data(d.L, xs, ys, d.es)
end

sf = ScalingFunction(log_shift;
    p_names=["T_c", "nu", "A"],
    x_scale="(T - T_c) * L^(1/nu)",
    y_scale="C - A * log(L)",
)
nothing # hide
```

The errors are not changed by a shift, so `d.es` is passed on unchanged. `x_scale` and
`y_scale` only describe the scaling in the output.

## The fit

```@example cv
sp = ScalingProblem(Ts, C, Ls;
    sf,
    p_space=[2.1:0.02:2.4, 0.6:0.1:1.4, 0.0:0.1:1.0],
    dx=[-2.0, 2.0],
)
```

All three parameters agree with the exact values ``T_c = 2.2692``, ``\nu = 1`` and
``A = 0.4945``.

```@example cv
using CairoMakie

colors = cgrad(:Blues)[range(0.35, 1.0; length=length(Ls))]
xs_s, ys_s, _, _ = scaled_data(sp)
fig = Figure(size=(780, 320))
ax1 = Axis(fig[1, 1]; xlabel=L"T", ylabel=L"C", title="data")
ax2 = Axis(fig[1, 2]; xlabel=L"(T - T_c)\,L^{1/\nu}", ylabel=L"C - A \ln L", title="collapse")
for i in eachindex(Ls)
    scatterlines!(ax1, Ts, C[:, i]; color=colors[i], markersize=5, label="L = $(Ls[i])")
    inside = -6 .<= xs_s[i] .<= 6
    scatter!(ax2, xs_s[i][inside], ys_s[i][inside]; color=colors[i], markersize=7)
end
vlines!(ax2, sp.dx; color=:gray, linestyle=:dash)
Legend(fig[1, 3], ax1)
fig
```

## Fixing parameters of a custom function

Parameters of custom scaling functions are fixed in the same way as those of the presets.
With ``T_c`` and ``\nu`` fixed to their exact values, only the amplitude is fitted:

```@example cv
sf_A = ScalingFunction(log_shift;
    p_names=["T_c", "nu", "A"], T_c=2 / log(1 + sqrt(2)), nu=1.0,
    x_scale="(T - T_c) * L^(1/nu)", y_scale="C - A * log(L)",
)
sp_A = ScalingProblem(Ts, C, Ls; sf=sf_A, p_space=[0.0:0.05:1.0], dx=[-2.0, 2.0])
```

To use the scaling form like a preset, e.g. `ScalingFunction(:log_shift)`, see
[Custom presets](@ref).
