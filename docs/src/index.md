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
      text: API Reference
      link: /api
    - theme: alt
      text: View on GitHub
      link: https://github.com/maltepuetz/ScalingCollapse.jl

features:
  - icon: "📉"
    title: One call, full fit
    details: Pass your data and system sizes; a parameter scan and Nelder–Mead refinement find the best collapse.
    link: /getting-started
  - icon: "🧩"
    title: Flexible scaling functions
    details: Use the :x, :xy and :xny presets, fix any parameter by name, or write your own scaling function.
  - icon: "📏"
    title: Parameter errors
    details: Profile, conditional or bootstrap error estimates for every fitted parameter.
---
```

```@meta
CurrentModule = ScalingCollapse
```

```@raw html
<div style="height: 3rem" aria-hidden="true"></div>
```

# What is ScalingCollapse.jl?

ScalingCollapse.jl automates finite-size scaling. Finite-size scaling
determines the critical parameters of a phase transition, such as the critical
temperature ``T_c`` and the exponents ``\nu`` and ``\beta``, from the way an
observable changes with the system size ``L``. Given measurements for several
sizes, the package finds the parameters for which all curves collapse onto a
single master curve.

!!! note "Development status"
    ScalingCollapse.jl is still in early development. Expect bugs and breaking
    changes, and feel free to open an issue if you run into a problem.

## Installation

```julia
using Pkg
Pkg.add("ScalingCollapse")
```

Continue with [Getting Started](@ref) for a first collapse.
