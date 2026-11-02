# Changelog

All notable changes to ScalingCollapse.jl are documented in this file.

Every pull request with a user-facing change adds an entry under `## Unreleased`.
A pull request that genuinely needs no entry — a typo fix, a CI-only change, a
dependency bump of a test or docs environment — carries the `skip-changelog` label
instead.

## Unreleased

## 0.3.0 - 2026-11-02

### Breaking changes

- The presets `:x`, `:xy` and `:xny` rescale x as `(x - p1) * L^(1/p2)` instead of
  `(x - p1)/p1 * L^(1/p2)` (#58). The new keyword `relative=true` selects the old form.
  A relative window `dx` equals `dx .* p1` in the new form, so to reproduce an old fit
  pass `relative=true` or multiply `dx` by the critical point.
- The default `quality` is `MultipleSplines()` instead of `Linear()`, and
  `MultipleSplines` interpolates with pchip (monotone piecewise-cubic Hermite) instead of
  cubic splines, also in the default `quality_scan` (#63). The old defaults are
  `quality=Linear()` and `quality_scan=MultipleSplines(scan_mode=true, interp=:cubic)`.
  pchip cuts peaks and dips off at the highest (lowest) measured point; for such data
  `interp=:cubic` is recommended.
- The default error analysis gives profile errors, which re-fit the other free
  parameters at every displacement and so account for correlations (#47). They are never
  smaller than the conditional errors of the same fit, the previous default, which
  `error=:conditional` still gives. The keyword `error` takes `true`/`:profile`,
  `:conditional`, `:bootstrap` or `false`.
- A custom scaling function that does not return a `Data` throws in the quality
  functions (#49). Returning `Data` was already required.

### Added

- Bootstrap errors with `error=:bootstrap` (#48): the fit is repeated on `n_bootstrap`
  (default 100) data sets resampled within the error bars, drawn from `rng`, and the
  error is the standard deviation of the refitted parameters. Needs error bars.
- `parameter_errors(sp; method, ...)` computes the errors of a solved problem without
  changing it, and `parameter_errors!(sp; ...)` stores them in `sp` (#59).
- `starting_ps` accepts several starting vectors; the final fit runs from each and keeps
  the best (#55).
- `StepSimplexer(steps)` sets the initial Nelder–Mead step of each free parameter, for
  use as `optimizer=Optim.NelderMead(initial_simplex=StepSimplexer(steps))` (#57).
- `MultipleSplines(interp=:pchip | :cubic)` (#63). pchip needs 2 points per system size,
  cubic splines 4.
- Docstrings for `QualityFunction`, `Data`, `n_parameters`, `scaled_p_names` and
  `fixed_p_names` (#64).

### Changed

- `SingleSpline` fits its spline only to the points in the overlap interval, with knots
  at quantiles of those points (one per 5 points). This removes an upward bias in `nu`
  and errors that were 10–100× too large (#50).
- With the default (affine) initial simplex, the inner fits of the profile errors are
  refined with steps of the conditional errors, so narrow valleys of the other
  parameters are resolved (#60). Before, the profile error could equal the conditional
  one; on the 2D Ising Binder data σ(ν) goes from 0.019 to 0.032.
- The quality functions and the error analysis are faster, with bit-identical results
  (#49): `MultipleSplines` with error bars about 18×, `Linear` 3.5–5×, `Houdayer` and
  `SingleSpline` about 2×.
- `scaled_data(...; splines=true)` uses the interpolation of the quality function (#63).
- A custom scaling function without `x_scale`/`y_scale` is shown as
  `custom scaling (set x_scale)` (#64).
- New dependency: `Random` (#48).

### Fixed

- `Linear` returned `NaN` when a point and both its neighbours have the same scaled x,
  e.g. on data uniform in `(x - x_c) L^(1/nu)` at `p1 = x_c`, so the fit missed the true
  minimum (#45). Coincident neighbours are now averaged.
- A `NaN` parameter no longer causes a `BoundsError` in the profile error analysis; all
  quality functions reject NaN intervals (#63).

### Documentation

- New documentation site with a tutorial, a manual (finite-size scaling, input data,
  scaling functions, fitting, error analysis, quality functions), worked examples with
  plots, and a grouped API reference (#56, #62, #64).
- Rewritten README and docstrings (#64).
- The error analysis page quantifies what the default `error_threshold = 1.3` means,
  from a benchmark on synthetic data and the cases with exactly known parameters (#65).

## 0.2.6 - 2026-10-29

### Fixed

- The Nelder–Mead optimization no longer stops at iteration 0 when the initial simplex
  reaches outside `p_space` (#42).
- The Nelder–Mead optimization no longer stalls with Optim ≥ 1.12 when initial simplex
  vertices have a non-finite quality (#43).
- The error analysis no longer hangs for parameters whose optimum is exactly zero; errors
  that cannot be determined are `NaN`, with a warning.
- Quality functions return `Inf` instead of `NaN` or negative values when too few points
  remain in the optimization interval.

### Changed

- Unknown keyword arguments, invalid options (`p_space`, `starting_ps`, `dx`,
  `error_threshold`), invalid data and custom scaling functions with a wrong signature
  throw clear `ArgumentError`s (#41).
- The parameter scan throws a clear error when the quality is infinite everywhere or has
  no isolated minimum.
- A warning is issued when Nelder–Mead stops without converging.
- Allow Optim 2 (#43).

## 0.2.5 - 2026-08-14

### Added

- The Nelder–Mead optimizer and its `Optim.Options` are configurable through the
  `optimizer` and `optimizer_options` keywords (#40).

### Changed

- Faster power scaling and a refactored data layout (#36, #37).

## 0.2.4 - 2025-11-26

### Added

- Custom presets for `ScalingFunction` by extending `ScalingCollapse.handle_preset` (#35).

## 0.2.3 - 2025-11-06

### Fixed

- Removed the unfinished InteractiveScaling extension, which caused compilation errors.

## 0.2.2 - 2025-04-15

### Fixed

- The error analysis gave a wrong error estimate when the quality was `NaN` at the
  initial parameters.

## 0.2.1 - 2024-04-16

### Fixed

- The interval boundaries in the `SingleSpline` quality function were mirrored.

## 0.2.0 - 2024-04-16

### Breaking changes

- Renamed the quality functions `Spline` to `MultipleSplines` and `SingleMasterCurve` to
  `Linear`.

### Added

- `SingleSpline` quality function, based on KissSmoothing.jl (#29).

## 0.1.1 - 2024-03-19

### Changed

- The optimization works with an arbitrary number of parameters, previously at most 3
  (#26).

## 0.1.0 - 2024-02-26

- Initial release: automatic finite-size scaling with a parameter scan and Nelder–Mead
  optimization, quality functions, error analysis, Measurements.jl support and plotting
  utilities.
