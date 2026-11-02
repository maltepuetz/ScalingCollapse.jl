```@meta
CurrentModule = ScalingCollapse
```

# Finite-size scaling

This page summarizes the idea behind a scaling collapse and how it maps onto the package.

## The scaling hypothesis

At a continuous phase transition, the correlation length diverges as
``\xi \sim |x - x_c|^{-\nu}``, where ``x`` is the control parameter (temperature, coupling,
probability, ...) and ``x_c`` its critical value. In a finite system of linear size ``L``,
``\xi`` is cut off by ``L``, and close to ``x_c`` an observable ``Y`` depends on ``x`` and
``L`` only through the ratio ``L/\xi``. This leads to the scaling form

```math
Y(x, L) = L^{\kappa/\nu} \, \mathcal{F}\big((x - x_c)\, L^{1/\nu}\big),
```

where ``\kappa`` is the critical exponent of ``Y`` (e.g. ``-\beta`` for the order parameter,
``\gamma`` for the susceptibility) and ``\mathcal{F}`` is a universal function. Rescaling
``x \to (x - x_c) L^{1/\nu}`` and ``Y \to Y L^{-\kappa/\nu}`` maps the data of all system sizes
onto ``\mathcal{F}``, a **scaling collapse**. The collapse only works for the right
parameters, so finding the best collapse determines ``x_c``, ``\nu`` and ``\kappa``.

## Choosing the scaling function

The presets of [`ScalingFunction`](@ref) cover the common observables. Each preset has two or
three parameters `p1`, `p2`, `p3`, which you can name with `p_names`:

| observable                                     | behaviour at ``x_c``          | preset | parameters                           |
|:-----------------------------------------------|:------------------------------|:-------|:-------------------------------------|
| dimensionless: Binder cumulant, ``\xi/L``, spanning probability | independent of ``L`` | `:x`   | ``x_c``, ``\nu``                     |
| order parameter (magnetization, cluster density) | ``\sim L^{-\beta/\nu}``     | `:xy`  | ``x_c``, ``\nu``, ``\beta``          |
| susceptibility, correlation time               | ``\sim L^{\gamma/\nu}``       | `:xny` | ``x_c``, ``\nu``, ``\gamma``         |

All presets rescale x as ``(x - p_1) L^{1/p_2}``. `:xy` multiplies y by ``L^{p_3/p_2}``,
`:xny` by ``L^{-p_3/p_2}``, so that the third parameter is a positive exponent in both cases.
Other scaling forms, e.g. with logarithms or corrections to scaling, need a
[custom scaling function](@ref "Custom scaling functions").

## Measuring the collapse

To compare the rescaled curves, a **quality function** interpolates them and measures their
spread; lower values mean a better collapse. The default, [`MultipleSplines`](@ref),
interpolates every system size separately and computes the variance between the curves on a
grid of points (weighted by the error bars, if given). The curves are only compared where all
system sizes overlap. The keyword `dx` restricts this further to a window of rescaled x
values: the scaling form holds best close to ``x_c``, i.e. for small
``|x - x_c| L^{1/\nu}``.

## Things that limit the accuracy

- **Corrections to scaling.** The scaling form is the leading term for large ``L`` and small
  ``|x - x_c|``. Small system sizes and wide windows `dx` bias the parameters. Compare fits
  with the smallest sizes removed and with a narrower window; see the
  [percolation example](@ref "Site percolation with error bars").
- **Correlated parameters.** ``x_c`` and ``\nu`` (and ``\nu`` and the y exponent) are usually
  correlated. The default profile error analysis takes this into account, see
  [Error analysis](@ref).
- **Correlated data.** Curves computed from the same Monte Carlo samples (e.g. histogram
  reweighting) are correlated. None of the error methods takes this into account.
