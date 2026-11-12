```@raw html
---
# https://vitepress.dev/reference/default-theme-home-page
layout: home

hero:
  name: "ScalingCollapse.jl"
  text: "Automatic finite-size scaling"
  tagline: Find the critical point and exponents that collapse data from many system sizes onto one master curve, in pure Julia.
  image:
    light: /logo-light.svg
    dark: /logo-dark.svg
    alt: ScalingCollapse.jl
  actions:
    - theme: brand
      text: Getting Started
      link: /getting-started
    - theme: alt
      text: Examples
      link: /examples/ising
    - theme: alt
      text: View on GitHub
      link: https://github.com/maltepuetz/ScalingCollapse.jl

features:
  - icon: "📉"
    title: One call, full fit
    details: Pass your data and system sizes; a parameter scan and a Nelder–Mead refinement find the best collapse.
    link: /getting-started
  - icon: "🧩"
    title: Flexible scaling forms
    details: Use the presets for the usual observables, fix any parameter by name, or write your own scaling function.
    link: /manual/scaling-functions
  - icon: "📏"
    title: Parameter errors
    details: Profile, conditional or bootstrap errors for every fitted parameter, and the quality landscape around the fit.
    link: /manual/errors
---
```

```@meta
CurrentModule = ScalingCollapse
```

```@raw html
<div style="height: 3rem" aria-hidden="true"></div>
```

# What is ScalingCollapse.jl?

Near a continuous phase transition, an observable ``Y`` measured in systems of linear size
``L`` follows a scaling form such as

```math
Y(T, L) = L^{-\beta/\nu} \, \mathcal{F}\big((T - T_c)\, L^{1/\nu}\big).
```

Plotted against ``(T - T_c) L^{1/\nu}`` and multiplied by ``L^{\beta/\nu}``, the curves of
all system sizes then fall onto one master curve ``\mathcal{F}``, but only for the right
values of ``T_c``, ``\nu`` and ``\beta``. ScalingCollapse.jl finds these values for you: it
measures how well the rescaled curves collapse, minimizes this measure over the parameters,
and estimates their errors.

```@example home
using ScalingCollapse
using Random

# synthetic data with T_c = 2.0, ν = 1.0 and β = 0.125
F(z) = 1 / (1 + exp(2z))
Ts = collect(range(1.8, 2.2; length=41))
Ls = [16, 32, 64, 128]
ys = [L^(-0.125) .* F.((Ts .- 2.0) .* L) .+ 0.002 .* randn(Xoshiro(L), 41) for L in Ls]

sp = ScalingProblem(Ts, ys, Ls;
    sf=ScalingFunction(:xy; p_names=["T_c", "nu", "beta"]),
    p_space=[1.9:0.02:2.1, 0.5:0.1:1.5, 0.0:0.05:0.5],
)
```

```@example home
using CairoMakie

colors = cgrad(:Blues)[range(0.35, 1.0; length=length(Ls))]
fig = Figure(size=(760, 320))
ax1 = Axis(fig[1, 1]; xlabel=L"T", ylabel=L"Y", title="data")
ax2 = Axis(fig[1, 2]; xlabel=L"(T - T_c)\,L^{1/\nu}", ylabel=L"Y\,L^{\beta/\nu}", title="collapse")
xs_scaled, ys_scaled, _, _ = scaled_data(sp)
for i in eachindex(Ls)
    scatterlines!(ax1, Ts, ys[i]; color=colors[i], markersize=6, label="L = $(Ls[i])")
    scatter!(ax2, xs_scaled[i], ys_scaled[i]; color=colors[i], markersize=7)
end
xlims!(ax2, -6, 6)
Legend(fig[1, 3], ax1)
fig
```

## Installation

ScalingCollapse.jl is a registered Julia package (Julia 1.9 or later):

```julia
using Pkg
Pkg.add("ScalingCollapse")
```

## Where to go next

- [Getting Started](@ref) walks through a first collapse step by step.
- The **Manual** explains the [input formats](@ref "Input data"), the
  [scaling functions](@ref "Scaling functions"), the [fit](@ref "Fitting"), the
  [error analysis](@ref "Error analysis") and the
  [quality functions](@ref "Quality functions").
- The **Examples** apply the package to real data: the [2D Ising model](@ref "2D Ising model"),
  [site percolation](@ref "Site percolation with error bars"), a
  [custom scaling form](@ref "Specific heat: a custom scaling function") and more.
- The [API Reference](@ref) documents every exported function.

!!! note "Development status"
    ScalingCollapse.jl is still young. Expect breaking changes between minor versions, and
    please [open an issue](https://github.com/maltepuetz/ScalingCollapse.jl/issues) if you
    run into a problem.

## Citing

If you use ScalingCollapse.jl in your research, please cite it with the metadata in
[`CITATION.cff`](https://github.com/maltepuetz/ScalingCollapse.jl/blob/main/CITATION.cff)
(GitHub's "Cite this repository" button exports it as BibTeX or APA) and state the
version you used. Every release is archived on Zenodo:
[10.5281/zenodo.10712426](https://doi.org/10.5281/zenodo.10712426) always resolves to the
latest version, and each version's record carries its own DOI for citing the exact release.
