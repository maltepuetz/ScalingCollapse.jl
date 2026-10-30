# Profile parameter errors (`error = :profile`, new default)

Draft PR description. Design: `docs/superpowers/specs/2026-10-30-profile-errors-design.md`.

## Summary

The default parameter errors were conditional: the error of parameter *i* was the displacement at
which the quality reaches `error_threshold` times its minimum while all other free parameters stay
fixed at their optimum. That ignores correlations between parameters and underestimates errors
(for example E3, `autoscale_site_order_exact`: sigma(p_c) = 2.6e-5 against autoScale's 7.5e-5).

This PR adds profile errors and makes them the default.

- For parameter *i* and displacement *x*, the profile quality is the quality minimized over all
  other free parameters with `p[i] = optimal_ps[i] + x` (a Nelder-Mead fit with the user's
  `optimizer` and `optimizer_options`; the nuisance parameters stay inside their `p_space`).
- The error is the larger of the two displacements, one on each side of the fit, at which the
  profile first reaches `error_threshold * minimum` when moving outward. The threshold rule and
  the default 1.3 are unchanged.
- The search starts at the conditional roots and grows outward by a factor of 1.5, then bisects to
  a relative precision of 1e-3. The profile at every point is started from the better of the
  conditional candidate and the nearest earlier point, so the profile is never above the conditional
  quality: a profile error is never smaller than the conditional error of the same fit (up to the
  bisection precision).
- Only a verified finite threshold crossing gives a finite error. Otherwise the error is NaN with a
  warning that names the reason (the profile never reaches the threshold within 60 growth steps; the
  rescaled data lose their overlap within `dx`; the quality at the optimum is not below the
  threshold).
- The profiled parameter is not bounded by `p_space` (as before); nuisance parameters are.
- With one free parameter the profile is the plain quality and both modes give identical errors.
- If a profile evaluation finds a quality lower than `sp.minimum` (the original fit had not
  converged), the fit is re-optimized from that point, a warning is issued, and the analysis is
  restarted (at most 3 restarts). `optimal_ps` and `minimum` can therefore change during the
  error analysis. If the restarts are exhausted, `optimal_ps` and `minimum` take the best point
  found, all errors are NaN, and a warning says the minimum could not be confirmed. A lower
  quality found outside `p_space` gives NaN for that parameter and a warning to widen `p_space`.
- The result is a local profile width: it measures the valley that contains the fit, not the
  projection of the whole region `quality <= threshold * minimum`.

These are still widths defined by a quality threshold, not standard deviations with a known
coverage. Profiling accounts for correlations between parameters; it does not calibrate the
threshold and does not account for correlated data points. This is stated in the docstring and
the README.

## API

`ScalingProblem(...; error = ...)` accepts a `Bool` or a `Symbol`:

| value | meaning |
|---|---|
| `true` (default) | same as `:profile` |
| `:profile` | local profile width (new default) |
| `:conditional` | previous method, same numbers (other parameters fixed) |
| `false` | no error analysis |

Any other value throws an `ArgumentError` listing the accepted values. `sp.error` stores the
normalized method. `:conditional` gives the same numbers as before: same bracket start and
precision 1e-10, same treatment of Inf as a crossing. The code was restructured into helpers shared
with the profile method (`_conditional_errors!` is no longer byte-for-byte the old text). Four
fixtures frozen before the change (commit 0cd93d5) check it: `==` on one machine (verified), rtol
1e-8 in CI.

## Rulings made while building

- A lower quality counts only if `q < minimum - max(1e-6 * minimum, optimizer_options.g_abstol)`.
  A purely relative 1e-6 is below Nelder-Mead's stopping noise and exhausted the restarts on real
  data (Susceptibility fits with `SingleSpline`).
- The first outward bracket step is `max(|conditional root|, 1e-6 * max(|p_i|, 1))`. A conditional
  root of about 1e-16 (a minimum sitting on a jump of `Linear`) otherwise made 60 growth steps cover
  only about 4e-6 and produced a spurious NaN.
- After a re-optimization, if the new minimum is above the lowest profile point already seen, that
  point is adopted: a worse minimum than one already found is never reported.
- `test/ising.jl`, Magnetization: the nu assertion was loosened from `atol=0.1` to `atol=0.15`. The
  restart finds a deeper minimum (quality 8.33e-6 to about 4.1e-6) separated from the original by a
  barrier, at nu of about 0.89 instead of 0.91. The error analysis is allowed to move the optimum,
  and this is such a case; the maintainer may prefer a different expectation or fit setup.
- Inner-fit retry: if the first inner fit did not converge and already started from the conditional
  candidate, the retry starts from where that fit stopped (retrying the identical fit is pointless).
- One warning per parameter per reason, not per side (less detail on which side failed).
- `_conditional_errors!` was refactored into shared helpers, accepted on the condition that the four
  Task 1 fixtures stay `==` on the same machine (verified).
- `test/error_methods.jl` loads the bond_strength and synthetic data with small inline HDF5 loaders
  instead of including `test/collapses/cases_*.jl`, to avoid duplicate `CASES` definitions (about
  10 lines of loader duplication). The `error` docstring entry was first reduced to the four accepted
  values and fully explained later.
- Restart and inner-fit warnings in the test log are left visible, because they report real
  non-convergence of those fits (list below).

## One-off checks

All three were run with a scratch script (not committed) over the whole collapse suite (211 cases,
603 parameter errors; one case has a single free parameter).

### 1. Profile versus conditional

Ratio sigma_profile / sigma_conditional over all 603 parameters:

| min | 10 % | median | 90 % | max |
|---|---|---|---|---|
| 0.625 | 1.06 | 1.41 | 2.21 | 7.28 |

No profile error is below conditional * (1 - 1e-3) when both are measured at the same fit.
The raw minimum (0.625) is one parameter, `learning_potts_vertical_g28` nu: the profile run
re-optimized the fit (quality 2.646 to 2.021), so its conditional counterpart was measured at a
different, worse fit. The conditional error at the re-optimized fit is 0.00380 against the
profile's 0.00530 (ratio 1.39). The same holds for the other three cases whose minimum moved (conditional errors rerun at the
profile run's final optimum with `starting_ps`; log in the scratchpad, `samefit.log`):

| case | minimum (conditional run to profile run) | profile / conditional at the same fit |
|---|---|---|
| learning_potts_vertical_g28 | 2.646 to 2.021 | 2.02, 1.39, 1.87 |
| synthetic_noise_err_SingleSpline | 90.36 to 71.25 | 2.23, 2.00, 1.18 |
| synthetic_corrections_SingleSpline (and `_sanity`, identical fit) | 88.35 to 70.97 | 2.87, 2.22, 1.26 |

Over the remaining 591 parameters the smallest ratio is 0.9995 and only 5 are below 1, all within the
bisection precision of 1e-3.

### 2. Stability

The profile errors were rerun with `optimizer_options = Optim.Options(g_tol = 1e-12)` for every
fit (original fit and inner fits) and the bisection precision at 1e-5 (a local edit in a separate
copy of the package, so `src/` of this branch is unchanged). Relative change of each error:

| median | 90 % | > 1 % | > 5 % | > 10 % | max |
|---|---|---|---|---|---|
| 0.02 % | 1.5 % | 67 of 603 | 26 | 8 | 89 % |

The largest change is `synthetic_three_sizes_MultipleSplines` x_c (2.2e-6 to 4.2e-6, 89 %). That
case is noise-free data whose quality minimum is about 5e-10; the tighter fit moves the minimum
itself (5.4e-10 to 4.3e-10), and at a minimum this close to zero the error follows the fit's
numerical noise. The other seven changes above 10 % are 16 % (`synthetic_noise_noerr_Linear` nu),
14 % (the same synthetic_three_sizes_MultipleSplines case, nu), 14 % and 11 %
(`synthetic_three_sizes_Linear`), 14 % (`curie_weiss_binder_exact` nu, 0.131 to 0.150), 11 %
(`synthetic_uneven_grid_Houdayer` p3) and 10 % (`synthetic_uneven_grid_MultipleSplines` x_c).
Errors are therefore stable to a few percent for most parameters, but not better than about 10 to
15 % for noise-free fits with a near-zero minimum, `Linear`/`Houdayer` landscapes and weakly
constrained parameters. The changes include the shift of the original fit, not only of the
profile search, so they overstate the profile-only sensitivity. No error became NaN or finite
differently between the two runs.

### 3. E3

`autoscale_site_order_exact`: sigma(p_c) = 7.60e-5 against autoScale's 7.5e-5 (+1.4 %, limit 20 %).
Conditional: 2.58e-5.

## Test-suite changes

- Collapse suite recalibrated: `sigmamax` of all exact cases re-derived from the profile errors
  (round(3 sigma, 2 significant digits)); `delta_abs` unchanged (every parameter passes
  `|p - p*| <= 3 sigma + delta_abs`, most now at 3 sigma alone); `[last]` in `references.toml`
  regenerated. The `[published]` check is now value-only (`check_errors=false`), because the stored
  published errors come from the conditional method.
- **Removed `[published]` entry:** `learning_potts_vertical_g28`. The profile run finds a lower
  quality (2.646 at the published optimum to 2.021 at the new minimum). Shifts relative to the
  published sigma: x_c +0.92 sigma, nu -0.62 sigma, p3 +0.85 sigma. This was the only `[published]`
  case that failed the value check. Comment in `test/collapses/cases_learning.jl`.
- **`bond_strength_exact`:** now fits L = 64 to 256 (was L = 16 to 256) and is no longer marked
  broken. Fit: p_c = 0.500301 +- 2.0e-4 (1.5 sigma), nu = 1.3173 +- 0.0102 (1.6 sigma),
  beta = 0.1267 +- 0.0076 (1.6 sigma). The twin `bond_strength` (L = 16 to 256) is unchanged.
- New `test/error_methods.jl` covers the methods, fixtures, lower-quality restarts, `p_space`
  bounds and warnings.

### Warnings in the test log

Left visible on purpose (11 in the last full run):

- Restart after a lower quality (7): `learning_potts_vertical_g28` (2.646 to 2.021),
  `synthetic_noise_err_SingleSpline` (90.36 to 71.25), `synthetic_corrections_SingleSpline` and
  `synthetic_corrections_SingleSpline_sanity` (88.35 to 70.97 each), and 3 restart warnings from 2
  fits in `test/ising.jl`: Magnetization (2 restarts, 8.33e-6 to 4.10e-6 to 4.06e-6) and one further
  fit (1 restart, 1.218e-6 to 1.205e-6; Susceptibility judging by its parameter values).
- "may be underestimated" (inner fits did not converge): `learning_ising_vertical_g1` (p3),
  `learning_potts_vertical_g1` (p2, p3).
- `learning_potts_vertical_g1` also emits one Nelder-Mead 1000-iteration warning.

Not covered by a dedicated test: the branch of `_reoptimize!` that adopts the profile's lowest point
when the re-optimized minimum is above it.

## Known limitation

A lower quality is only detected at points the profile evaluates, which lie beyond the conditional
roots. A non-converged fit is therefore caught when the better minimum lies along a tilted valley
that the profile walks into, not when it is a small offset in a round bowl around the reported
optimum.

## Runtime

The full `Pkg.test()` takes 9 min 17 s (1313 tests pass), of which the collapse suite takes
522 s; the suite took about 100 s before. This increase (profile errors refit the other
parameters at every bisection step) was accepted in the design.

## Not in this PR

- No version bump; `Project.toml` is unchanged. The bump (0.3.0) follows separately.
- Bootstrap errors, a covariance output, any change to `error_threshold`.

Not yet run: CI on all targets (Julia lts, 1, nightly; Linux, macOS aarch64, Windows including
32-bit). Only a local run on macOS with Julia 1.13 was done.

🤖 Generated with [Claude Code](https://claude.com/claude-code)
