# Scaling

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://maltepuetz.github.io/ScalingCollapse.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://maltepuetz.github.io/ScalingCollapse.jl/dev/)
[![Build Status](https://github.com/maltepuetz/ScalingCollapse.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/maltepuetz/ScalingCollapse.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/maltepuetz/ScalingCollapse.jl/branch/main/graph/badge.svg?token=G2BD929KV0)](https://codecov.io/gh/maltepuetz/ScalingCollapse.jl)
[![Aqua](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)


This is a package for automatic finite size scaling. Finite size scaling is a method to determine the critical parameters of a phase transition by exploiting, so called, finite size effects. This package implements an automatic optimization algorithm to find the best parameters for a finite size scaling analysis automatically.

It is still in early development, so expect bugs and breaking changes. Feel free to create issues in case you stumble upon any problems!

## Installation
For installation run the following command:
```julia
julia> using Pkg; Pkg.add("https://github.com/maltepuetz/Scaling.jl.git")
```

## Getting started
In order to find the best parameters for a finite size scaling analysis, we need to create a `ScalingProblem`. Lets say we have data for the binder cumulant of the 2D Ising model `binder` for different temperatures `Ts` and different system sizes `Ls`. We want to find the critical temperature and the critical exponent of the correlation length. Let's create a `ScalingProblem` as follows:
```julia
using Scaling
sp = ScalingProblem(Ts, binder, Ls;
    sf = ScalingFunction(:x),
    dx = [-1.0, 1.0],
)
```
The kwarg `sf` is used to set the scaling function. In this case we use the scaling function preset `:x` which only rescales the x-axis as `x -> (x - p1)/p1 * L^(1/p2)`. We can also give the parameters custom names by passing the kwarg `p_names::Vector{String}` to the `ScalingFunction` constructor.
We use the kwarg `dx` to limit the optimization to a fixed interval because powerlaw dependece is usually only accurate close to the critical point. The default value is `dx = [-Inf, Inf]`.
Other kwargs are:
- `p_space = [0.1:0.1:3.0 for _ in 1:n_parameters(sf)]]` This sets the parameter space that will be scanned for local minima of the cost function. The optimization algorithm will start from the determined local minima to find the global minimum. Note that the the optimal parameters should lie within the given parameter space. Also note, that the initial scan can take a while for too small stepsizes. If you have a good idea where the optimal parameters should lie, you can set the parameter space accordingly. You can also skip the local minima search completely, by setting the kwarg `starting_ps`.
- `starting_ps::Vector{Float64}` This kwarg can be used to skip the local minima search. If `starting_ps` is set to a specific parameter set, the optimization algorithm will start from this point. This can be useful if you already have a good guess for the optimal parameters. Note that you might not find the optimal parameters if you start from a bad guess.
- `error::Bool=false` This kwarg can be used to turn the error calculation on or off. If `error=true`, the error of the parameters will be calculated using a `2S` method. This means that we fix all but one parameter and vary the free paramter, such that we double the cost function value. The error is then the maximum of the variation in plus and minus direction.
- `verbose::Bool=false` This kwarg can be used to turn the verbose output on or off. If `verbose=true`, the algorithm will print some information about the optimization process.

There are also other preset ScalingFunctions:
- `:x` Rescales the x-axis as `x -> (x - p1)/p1 * L^(1/p2)`.
- `:xy` This is the default scaling function. It rescales the x-axis as `x -> (x - p1)/p1 * L^(1/p2)` and the y-axis as `y -> y * L^(p3/p2)`.
- `:xny` This scaling function is similar to `:xy` but it rescales the y-axis as `y -> y * L^(-p3/p2)`. In this case the critical exponent `p3` has a negative sign.

Note that you can also create your own scaling function! We can make our scaling function looking a little bit cooler by giving our parameters custom names with the kwarg `p_names::Vecor{String}`. 

Now that we know all that lets do another scaling problem. This time we know that the critical temperature is somewhere between `T_c = 2.0` and `T_c = 2.4` and the critical exponent is somewhere between `nu = 0.5` and `nu = 1.5`. We can set the parameter space accordingly:
```julia
using Scaling
sp = ScalingProblem(Ts, binder, Ls;
    sf = ScalingFunction(:x; p_names = ["T_c", "nu"]),
    dx = [-1.0, 1.0],
    p_space = [2.0:0.1:2.4, 0.5:0.1:1.5],
    error = true,
)
```
We limited the parameter space so the algorithm runs faster, calculated some error bars and made everything look cooler by giving the parameters custom names! 
### Fixing parameters
Let's say that we know from our statistical mechanics course that the critical exponent of the correlation length is `nu = 1`. We can fix this parameter by passing it as a kwarg to the `ScalingFunction` constructor.
```julia
using Scaling
sp = ScalingProblem(Ts, binder, Ls;
    sf = ScalingFunction(:x; p_names = ["T_c", "nu"], nu=1),
    dx = [-1.0, 1.0],
    p_space = [2.0:0.1:2.4],
    error = true,
)
```
Now the algorithm will only optimize the critical temperature `T_c` and the critical exponent of the correlation length will be fixed to `nu = 1`.

## Creating custom scaling functions
To be written...


## Have fun and let me know if you have any problems!

