# Profile parameter errors: design

Date: 2026-10-30

## Why

Package finding 2 of PR #44: `error_analysis!` moves one parameter while holding all others at
their optimum, and reports the displacement at which `quality / minimum` reaches
`error_threshold`. That is the width of the quality minimum along one axis. When parameters are
correlated, the minimum is a tilted valley, and the axis width can be far smaller than how far a
parameter can really move. The collapse suite shows it:

- `bond_strength`: β moves 5–10σ across size sets, in lockstep with p_c.
- E3 (`autoscale_site_order_exact`): σ(p_c) = 2.6e-5, about 3× smaller than autoScale's 7.5e-5
  for the same data.

## Goal

Make the default error a **local profile width**. For free parameter *i*, the profile quality
at displacement *x* is the quality minimized over all other free parameters with
`p[i] = optimal_ps[i] + x`. The error is the larger of the two displacements, one on each side
of the fit, at which the profile quality **first** reaches `error_threshold · minimum` when moving
outward from the fit.

This is deliberately local: it measures the valley that contains the fit, not the projection of
the whole region `{quality ≤ threshold · minimum}`, which can include separate valleys. The
threshold rule itself (default 1.3) is unchanged.

What these numbers are and are not (also stated in the docstring and README): they are widths
defined by a quality threshold, not standard deviations with a known coverage. Profiling
accounts for correlations between parameters. It does not calibrate the threshold 1.3, and it
does not account for correlated data points (e.g. curves whose points come from the same Monte
Carlo samples).

## Non-goals

- The version bump to 0.3.0 (separate PR, after this one).
- Bootstrap errors (a later, opt-in feature; see "Future").
- A covariance or correlation output.
- Changing the meaning or default of `error_threshold`.
- Fixing the Inf-as-crossing behavior of the conditional method (see "Valid crossings"); it
  stays bit for bit as today.
- Findings 1, 3 and 4 of PR #44 (finding 1 is already fixed in #45).

## Anti-goals

- A profile error smaller than the conditional error for the same fit (beyond the bisection
  precision). The algorithm enforces this explicitly (see "Profile quality").
- A finite error that does not come from a verified finite threshold crossing.
- Numerically valid-looking errors around a point that could not be confirmed as a minimum.

## Spike evidence (2026-10-29, on `main` at 4ea8837, before #45)

A local patch ran the whole collapse suite (211 cases) with both methods:

| | conditional (today) | profile |
|---|---|---|
| E3 σ(p_c) | 2.6e-5 | 7.6e-5 (autoScale: 7.5e-5) |
| `bond_strength` σ(p_c), σ(β) | 1.1e-4, 0.0031 | 4.8e-4, 0.012 |
| ratio profile / conditional, all parameters | – | min 1.00, median 1.42, 90 % < 2.2, max 7.3 |
| one free parameter | – | identical |
| collapse suite runtime | 103 s | 478 s |

Two failure modes of a naive implementation showed up:

1. **Basin hopping.** Today's bracket search starts at 10 % of the parameter value. With
   re-fitting, that first step can land past a barrier in a different valley whose quality is
   below threshold, so the bracket keeps growing. This gave σ(x_c) = 0.14 for x_c = 0.7
   (`synthetic_noise_err_Linear`, 570× today's error). Starting from the conditional root and
   growing slowly fixed it (4.3×).
2. **The optimum is not a minimum.** In 5 cases (mostly SingleSpline, also
   `synthetic_three_sizes_MultipleSplines` and `learning_potts_vertical_g28`) the profile found
   quality below `sp.minimum`: the final optimization had not converged. The threshold may then
   never be reached, giving NaN.

The spike bounded the profiled parameter by `p_space` (it used `check_bounds=true` on the full
vector), which this design does not; see "Bounds".

The runtime increase is accepted by the maintainer.

## Design

### API

`error` accepts a `Bool` or a `Symbol`:

| value | meaning |
|---|---|
| `true` (default) | same as `:profile` |
| `:profile` | local profile width (new default) |
| `:conditional` | today's method, bit for bit (other parameters fixed) |
| `false` | no error analysis |

Any other value throws an `ArgumentError` that lists the accepted values. The field `sp.error`
stores the normalized method; the representation is the builder's call (e.g. a `Symbol` with
`:none` for `false`), but `Base.show` and every reader of `sp.error` must keep working.

`:conditional` must reproduce today's errors exactly (same code path, same bracket start, same
precision 1e-10), including today's treatment of Inf as a crossing.

### One free parameter

With one free parameter, the profile is the plain quality. `:profile` then returns the
conditional result unchanged (no second bisection), so both modes give identical errors.

### Bounds

- The profiled parameter *i* is **not** bounded by `p_space`, as in today's conditional analysis.
- The other ("nuisance") parameters are kept inside their `p_space` bounds by the inner fit:
  the inner objective returns Inf when any nuisance coordinate is outside its bounds, and
  otherwise evaluates `quality(sp, p; check_bounds=false)`. (`check_bounds=true` cannot be used,
  because it also checks parameter *i*.) The initial simplex still goes through
  `BoundedSimplexer` on the nuisance bounds.
- If the box truncates the valley, the profile width comes out smaller than the unbounded one,
  but never smaller than the conditional width.

### Profile quality

`prof(x)` for parameter *i* at displacement *x*:

1. Evaluate the **conditional candidate**: the quality with the nuisance parameters at their
   optimal values.
2. Take the **continuation candidate**: the re-fitted nuisance values of the previously
   evaluated displacement closest to *x* on the same side (initially the optimal values).
3. Start the inner Nelder–Mead (`sp.optimizer.parameters`, `sp.optimizer_options`) from
   whichever candidate has the lower quality.
4. `prof(x)` = the minimum of the inner fit's result and both candidates' qualities. This
   enforces `prof(x) ≤ conditional quality at x`, so the profile width is ≥ the conditional
   width up to the bisection precision.
5. **Inner-fit failure.** If the inner fit does not converge, retry it once from the
   conditional candidate. If it still does not converge, use the best value found and count the
   failure; when any inner fit of a parameter failed, warn once for that parameter (its width
   may be underestimated). It is not set to NaN.

Each evaluation is classified as **valid finite**, **invalid** (Inf: no nuisance values with
overlapping data within `dx`), or NaN (treated as invalid).

### Bracketing and bisection (per side)

1. **Conditional roots first.** Run today's bracket and bisection on the conditional quality to
   get `d₋ < 0 < d₊`. If a conditional root fails, the parameter's error is NaN with today's
   warning.
2. **Bracket outward from the conditional root.** Starting at `x = d±`, multiply by 1.5 until
   the first valid finite `prof(x)` with `prof(x) / minimum ≥ error_threshold[i]`. The last
   point below threshold is the inner end. Moving outward step by step from the fit, and
   stopping at the first crossing, is what keeps the result local. It cannot guarantee that no
   barrier narrower than one step is skipped; the 1.5 growth factor keeps that risk small (it
   is no worse than today's doubling).
3. **Valid crossings.** If the outward step reaches an invalid point before a valid crossing,
   bisect between the inner end and the invalid point to find a valid finite point at or above
   the threshold. If none is found within the bisection budget, the error is NaN with a warning
   that the quality became undefined (no overlap within `dx`) before reaching the threshold.
   A finite error is reported only from a verified finite crossing.
4. **Bisect** between the inner and outer ends, keeping only valid finite points as ends, with
   a relative precision of about 1e-3 (the exact value is the builder's call).
5. Give up after a bounded number of outward steps (NaN, warning), as `_bracket` does today.
6. The error is the larger of the two absolute widths, as today.

### Lower quality found during the profile

Every `prof` evaluation is compared with `sp.minimum`. If one is lower by more than a small
relative tolerance (the spike used 1e-6), the original fit was not optimal:

- **The full parameter vector is inside `p_space`** (parameter *i* included): re-run the final
  optimization from it (`optimize_starting_ps!` with that point as start), update
  `optimal_ps` and `minimum`, and restart the error analysis from scratch for all parameters
  (no cached profile points, brackets or partial errors survive a restart). Each restart warns
  with the old and new minimum.
- **Parameter *i* is outside `p_space`:** this does not show that the bounded fit failed. Do
  not restart. Set this parameter's error to NaN and warn that a lower quality exists outside
  `p_space` and that `p_space` may need widening.
- **Restarts exhausted** (a small fixed number, builder's call, e.g. 3): keep the best feasible
  fit found (lowest quality inside `p_space`) as `optimal_ps` / `minimum`, set **all** errors
  to NaN, and warn that the minimum could not be confirmed.

The docstring states that the error analysis can change `optimal_ps` and `minimum`.

### Docs

- `ScalingProblem` docstring (`src/scalingproblem.jl`) and README: the `error` values, the
  local-profile definition, that the analysis may move the optimum, and the "what these numbers
  are" paragraph from "Goal". Fix the README's stale "2×" wording to the actual default 1.3.
- `CLAUDE.md`, Architecture step 4.

### Tests

Unit tests (extend the "error analysis" testset in `test/error_handling.jl`, which already uses a
quadratic quality with error `√0.3`):

- **Correlated quadratic:** quality `1 + a² + b² + 2ρab`, threshold 1.3: conditional error `√0.3`
  for both parameters, profile error `√(0.3 / (1 − ρ²))`, compared with `rtol ≈ 1e-3`.
- **One free parameter:** `:profile` and `:conditional` give identical (`==`) errors.
- **Two valleys:** a quality with a second, deeper-than-threshold valley beyond a barrier on one
  side: the profile width stops at the first crossing and does not extend into the second
  valley.
- **Lower minimum inside `p_space`:** a reported optimum that is not the minimum triggers a
  restart and a warning and ends at the true minimum.
- **Lower minimum outside `p_space`:** NaN for that parameter, the `p_space` warning, no
  restart.
- **Restart exhaustion:** with the restart limit lowered to 0 inside the test, all errors are
  NaN, with the warning, and `optimal_ps` is the best feasible point.
- **Invalid region:** a quality that becomes Inf before reaching the threshold gives NaN and the
  "undefined" warning in `:profile` mode.
- **API:** an invalid `error` value throws `ArgumentError`; `error=false` skips the analysis.
- **`:conditional` compatibility:** before changing `src/`, record today's errors for a handful
  of fits as literal numbers: the Ising MC Binder fit (`test/ising.jl`), CDW, `bond_strength`
  and one three-parameter synthetic case. Test `error=:conditional` against them with `==`.
  No parallel conditional references for the whole suite.

Collapse suite (`test/collapses/`):

- Regenerate `[last]` with `update_references.jl`.
- Recalibrate exact-case `σmax` to ≈ 3× the new fitted error (the rule from PR #44). Keep `δabs`
  as is unless a case needs it under the new errors; any change is re-measured and recorded in
  the case comment.
- `[published]` stays frozen, but its stored errors come from the conditional method. The check
  against `[published]` keeps the value check (`|p − p_ref| ≤ 0.25σ_ref`) and drops the
  error-ratio check; the `[last]` check keeps both.
- **A `[published]` fit that moves (maintainer's decision).** If re-optimization moves a case with
  a `[published]` entry so that its value check fails (e.g. `learning_potts_vertical_g28`),
  remove that case's `[published]` entry and record the old and new quality and the shift in σ
  in a comment. Word it as "a lower quality was found", which shows the published fit was not
  optimal; it does not show that the new fit is the global minimum. Keep checking the case
  against `[last]` only. This is the precedent of R5 cut 40 (`learning_potts_horizontal_b40`)
  in PR #44. List every removed entry in the PR description.
- **`bond_strength_exact` (maintainer's decision):** fit only L = 64–256 and remove `broken`.
  Recalibrate its `σmax` like the other exact cases; `δabs` = 0. The regression twin
  `bond_strength` keeps L = 16–256. Rewrite the case comment with the measured cause. Size
  windows with profile errors (spike, before #45):

  | sizes | p_c (0.5) | ν (4/3) | β (5/36) | quality |
  |---|---|---|---|---|
  | 16–256 | 0.49975 ± 0.00048 | 1.284 ± 0.015 (3.2σ) | 0.1545 ± 0.0122 (1.3σ) | 2.13 |
  | 32–256 | 0.49904 ± 0.00056 | 1.272 ± 0.018 (3.4σ) | 0.177 ± 0.016 (2.4σ) | 1.40 |
  | 64–256 | 0.50030 ± 0.00020 (1.5σ) | 1.317 ± 0.010 (1.6σ) | 0.1267 ± 0.0076 (1.6σ) | 0.07 |

  The β miss on 16–256 was the conditional error ignoring the p_c–β correlation (σ(β) 4× too
  small). The ν miss is corrections to scaling at L = 16–32 (open-boundary largest-cluster
  fraction): the quality of ~2 shows that the curves do not collapse within their error bars
  there. The spanning probability from the same samples has smaller corrections. Caveat for the
  comment: the 81 points of one curve come from the same samples and are strongly correlated,
  which the quality treats as independent (hence quality ≪ 1 on 64–256).

### One-off checks for the PR (not permanent tests)

- No profile error is smaller than the conditional error for the same fit, over the whole suite.
- Stability: rerun the suite's profile errors with tighter inner tolerances and bisection
  precision; report the largest relative change in the PR description.
- E3's profile σ(p_c) is within 20 % of autoScale's 7.5e-5.

### Done when

- Full `Pkg.test()` passes on all CI targets (including 32-bit Windows).
- The one-off checks above are run and reported in the PR description.
- Docs updated.

## Left to the builder

- Inner Nelder–Mead tuning (initial simplex size, tolerances) for speed, as long as the
  stability check passes.
- The exact bisection precision, the lower-quality tolerance and the restart limit.
- How `sp.error` is stored, and how the restart limit is made lowerable in tests.
- Code structure inside `src/error_analysis.jl`.

## Future (not this PR)

- `error = :bootstrap` with `n_bootstrap` and an RNG keyword: resample `y + N(0, e)` and refit
  from the optimum. Only valid with error bars. The symbol-valued `error` keyword leaves room
  for it.
- The conditional method's Inf-as-crossing behavior, if `:conditional` is ever allowed to
  change.
