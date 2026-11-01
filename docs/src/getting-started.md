```@meta
CurrentModule = ScalingCollapse
```

# Getting Started

To find the best parameters for a finite size scaling analysis, we need to create a `ScalingProblem`. Lets say we have data for the binder cumulant of the 2D Ising model `binder` for different temperatures `Ts` and different system sizes `Ls`. We want to find the critical temperature and the critical exponent of the correlation length. Let's create a `ScalingProblem` as follows:
```julia
using ScalingCollapse
sp = ScalingProblem(Ts, binder, Ls;
    sf = ScalingFunction(:x),
    dx = [-1.0, 1.0],
)
```
The kwarg `sf` is used to set the scaling function. In this case we use the scaling function preset `:x` which only rescales the x-axis as `x -> (x - p1)/p1 * L^(1/p2)`. We can also give the parameters custom names by passing the kwarg `p_names::Vector{String}` to the `ScalingFunction` constructor.
We use the kwarg `dx` to limit the optimization to a fixed interval because powerlaw dependece is usually only accurate close to the critical point. The default value is `dx = [-Inf, Inf]`.
Other kwargs are:
- `p_space = [0.1:0.1:3.0 for _ in 1:n_parameters(sf)]` This sets the parameter space that will be scanned for local minima of the cost function. The optimization algorithm will start from the determined local minima to find the global minimum. Note that the the optimal parameters should lie within the given parameter space. Also note, that the initial scan can take a while for too small stepsizes. If you have a good idea where the optimal parameters should lie, you can set the parameter space accordingly. You can also skip the local minima search completely, by setting the kwarg `starting_ps`.
- `starting_ps::Vector{Float64}` This kwarg can be used to skip the local minima search. If `starting_ps` is set to a specific parameter set, the optimization algorithm will start from this point. This can be useful if you already have a good guess for the optimal parameters. Note that you might not find the optimal parameters if you start from a bad guess.
- `error=true` This kwarg selects the error analysis: `true` or `:profile` (default), `:conditional`, `:bootstrap`, or `false` (no error analysis). With `:profile`, the error of each parameter is a local profile width: we displace the parameter from its optimum, re-fit all other free parameters at every displacement, and find where the quality reaches `error_threshold` times its minimum (default `error_threshold = 1.3`) when moving outward from the fit. The error is the maximum of the displacements in plus and minus direction. With `:conditional`, the other parameters are instead held fixed at their optimal values, which ignores correlations and gives smaller or equal errors (equal with one free parameter or uncorrelated parameters). If the profile finds a lower quality than the reported minimum, the fit is re-optimized and the analysis restarted (with a warning), so `optimal_ps` and `minimum` can change during the error analysis. If the restarts run out, all errors are `NaN` with a warning; a lower quality found with the parameter outside `p_space` gives `NaN` for that parameter and a warning to widen `p_space`. Note that these errors are widths defined by a quality threshold, not standard deviations with a known coverage: profiling accounts for correlations between parameters, but the threshold is not calibrated and correlated data points (e.g. curves built from the same Monte Carlo samples) are not accounted for.
  With `:bootstrap`, the errors are standard deviations that come from the error bars of the data: we draw `n_bootstrap` (default 100) synthetic datasets, in which every y value is replaced by a Gaussian random number centred on it with its error bar as width, re-fit each one starting from the optimal parameters, and report the standard deviation of the re-fitted parameters. This requires error bars and assumes they are Gaussian and independent (correlated data points are again not accounted for). Pass a seeded random number generator, e.g. `rng = Random.Xoshiro(42)`, for reproducible errors. The bootstrap costs about `n_bootstrap` times one fit. The bootstrap errors scatter between random seeds by a relative amount of about `1/sqrt(2(n_bootstrap - 1))`, i.e. about 7% for the default of 100, which is enough to quote errors to one significant digit; use about 500 or more for a reliable second digit.
- `verbose::Bool=false` This kwarg can be used to turn the verbose output on or off. If `verbose=true`, the algorithm will print some information about the optimization process.

There are also other preset ScalingFunctions:
- `:x` Rescales the x-axis as `x -> (x - p1)/p1 * L^(1/p2)`.
- `:xy` This is the default scaling function. It rescales the x-axis as `x -> (x - p1)/p1 * L^(1/p2)` and the y-axis as `y -> y * L^(p3/p2)`.
- `:xny` This scaling function is similar to `:xy` but it rescales the y-axis as `y -> y * L^(-p3/p2)`. In this case the critical exponent `p3` has a negative sign.

Note that you can also create your own scaling function! We can make our scaling function looking a little bit cooler by giving our parameters custom names with the kwarg `p_names::Vector{String}`. 

Now that we know the basics, let us solve another scaling problem. This time we know that the critical temperature is somewhere between `T_c = 2.0` and `T_c = 2.4` and the critical exponent is somewhere between `nu = 0.5` and `nu = 1.5`. We set the parameter space accordingly:
```julia
using ScalingCollapse
sp = ScalingProblem(Ts, binder, Ls;
    sf = ScalingFunction(:x; p_names = ["T_c", "nu"]),
    dx = [-1.0, 1.0],
    p_space = [2.0:0.1:2.4, 0.5:0.1:1.5],
    error = true,
)
```
We limited the parameter space so the algorithm runs faster, calculated some error bars and made everything look cooler by giving the parameters custom names! 
## Fixing parameters
Let's say that we know from our statistical mechanics course that the critical exponent of the correlation length is `nu = 1`. We can fix this parameter by passing it as a kwarg to the `ScalingFunction` constructor.
```julia
using ScalingCollapse
sp = ScalingProblem(Ts, binder, Ls;
    sf = ScalingFunction(:x; p_names = ["T_c", "nu"], nu=1),
    dx = [-1.0, 1.0],
    p_space = [2.0:0.1:2.4],
    error = true,
)
```
Now the algorithm will only optimize the critical temperature `T_c` and the critical exponent of the correlation length will be fixed to `nu = 1`.
