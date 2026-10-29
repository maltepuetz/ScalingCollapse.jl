# Collapse benchmark suite and the Optim ≥ 1.12 fix: design

Date: 2026-10-29

## Why

Finite-size scaling is subtle and very data-dependent, so it is hard to tell when a change to
ScalingCollapse.jl breaks something. Today the tests contain a handful of collapses (2D Ising Monte
Carlo, the Holstein CDW demo). That is far too few to notice when results quietly get worse.

This already happened. While surveying data for this design we found a live bug. With Optim ≥ 1.12,
which is what a fresh install resolves to today, the final Nelder–Mead step stops at iteration 0 on
real data. It then returns the best grid point of the parameter scan instead of the optimum, with
error bars about 15× too large. The current test suite passes anyway.

Here is the reproduction. It uses Fig. 3 of arXiv:2505.22720 (qubit loss), with the paper's own
`ScalingProblem` call:

| Package version | Optim | θ = 0 (published: r_c = 0.336814(21), ν = 1.5316(40)) | θ = π/2 (exact: ν = 4/3) |
|---|---|---|---|
| v0.2.2, current `main` | 1.11.0 | 0.336814(21), 1.5316(40) ✓ | 0.499984(43), 1.3370(59) ✓ |
| v0.2.2, current `main` | 1.12.0 / 1.13.3 | 0.337(35), 1.5(59) ✗ | 0.5(19), 1.3775(883) ✗ |

The merged PR #42 (`BoundedSimplexer`) keeps the initial simplex inside `p_space`. That does not fix
this case: the vertex `[0.575, 1.3775]` lies inside `p_space`, but the quality there is `Inf`,
because the overlap interval intersected with `dx` is empty. Optim ≥ 1.12 then reports a non-finite
convergence measure and stops. Optim 1.11 kept iterating.

## Goals

1. **Many collapses that cover a lot of variety.** Different physics, numbers of sizes (3 to 7),
   noise levels, presets (`:x`, `:xy`, `:xny`), quality functions, and both clean and messy data.
2. **Correctness where the answer is known.** If the true parameters are known exactly, or known
   very precisely from the literature, the test checks the fit against them.
3. **Regression detection where it is not.** Otherwise the test compares against a stored previous
   result and flags any drift.
4. **Fix the Optim ≥ 1.12 failure**, with this suite as its regression test.
5. **Fits in CI:** the whole test run stays under about 10–15 min per CI job. Up to tens of MB of
   test data is acceptable.

## Non-goals

- Performance benchmarking. Timing is only watched to stay inside the CI budget.
- Comparing our numbers against other FSS tools' outputs. We only reuse their input data.
- Plotting or visual checks.
- Pinning Optim to 1.11 as the fix. We fix the root cause instead (see below).

## Anti-goals (failure even if it "works")

- **Tolerances so loose that a broken fit passes.** The Optim bug moved ν by about 0.04 and
  inflated the errors 15×. The pass criteria must be tight enough to catch a shift or error
  inflation of that size in every case. This is checked by auditing the tolerances, not by making
  every case reproduce the bug: many cases (noise-free, wide `dx`) never hit a non-finite simplex
  vertex. *(Claude's call.)* Only E1 and R1 θ = 0 must be shown to fail on unfixed code with
  Optim 1.13. Any other case that happens to fail there is a bonus and is recorded in the PR.
- **Flaky tests.** All data is committed, any generation uses fixed seeds, and nothing is downloaded
  at test time.
- **Data in the repo that we are not allowed to redistribute** (see Licensing).
- **Reference values that get silently overwritten.** Updating stored results must be an explicit,
  reviewable step.

## Part 1: Fix the Nelder–Mead start (do this first)

The cause is an initial simplex vertex where the quality is non-finite. A vertex can be inside
`p_space` and still be `Inf`, for example when the overlap interval inside `dx` is empty.

**Approach (Claude's call).** Generalize PR #42's idea from "inside `p_space`" to "finite quality".
When building the initial simplex, evaluate the quality at each vertex. If it is not finite, move
the vertex toward the starting point by repeatedly halving the step, up to a fixed number of times
(for example 20). If a vertex is still non-finite after that, keep the current warning path. This
requires the simplex construction to know the objective. `_bounded_optimizer(sp)` may therefore
become something like `_initial_simplex(sp, f, x0)`, passed to Optim through a custom `Simplexer`
that closes over `f`. The exact mechanics are left to the builder. The user-supplied
`sp.optimizer.initial_simplex` must stay the base simplex that gets shrunk, so #39's configuration
option keeps working.

**Acceptance:**
- With Optim 1.13.x, the qubit-loss θ = 0 and θ = π/2 cases reproduce the Optim 1.11 numbers above:
  equal to within 1% of the error bar, or as close as Optim's changes allow, with tolerances
  documented in the test.
- All existing tests still pass.
- **Ships first, on its own:** a separate PR and a 0.2.6 patch release, before the suite
  (decided by the user).
- **Widen compat to Optim 2:** in the same PR, change `Optim = "1.6.2"` to `Optim = "1.6.2, 2"`.
  This was checked on 2026-10-29 with current `main` and the compat widened in a scratch copy:
  Optim 2.3.2 is API-compatible (all 164 existing tests pass), and it fails on the qubit-loss cases
  exactly like 1.12+. So the same fix must make E1 and R1 θ = 0 pass on both 1.13.x and 2.x.
- CI must run with the newest Optim the compat allows (this is already the default). Optionally,
  add one CI job that uses `julia-actions/julia-downgrade-compat` so that the lowest compatible
  versions are tested too. This part is the builder's call.

## Part 2: The suite

### Layout

```
test/
  collapses/
    runcollapses.jl          # included from runtests.jl; loops over all cases
    cases_exact.jl           # case definitions with known answers
    cases_regression.jl      # case definitions checked against stored references
    references.toml          # stored results for regression cases (generated, reviewed)
    update_references.jl     # explicit script to regenerate references.toml; never run in CI
  testdata/
    README.md                # source, citation, license for every file
    <source>/*.h5            # converted fixtures (HDF5; HDF5 is already a test dependency)
  generators/                # scripts that produce the self-generated datasets (fixed seeds);
                             # their outputs are committed, the scripts are not run in CI
```

Exact names are the builder's call. The separation between case definitions, data, references,
generators and the reference update script is required.

### A case

Each case is a small record:

- `name`: a unique, readable name, for example `qubitloss_theta_pi2_percolation`.
- **Data loader:** reads the fixture and returns `xs, ys, es, Ls`.
- `ScalingProblem` kwargs: exactly as published where the case comes from a paper.
- `kind`: `:exact` or `:regression`.
- For `:exact`, the known values for some or all free parameters. Unlisted parameters are not
  checked against truth.
- `source`: citation or URL. Its license appears in `testdata/README.md`.

Cases are run with `@testset` per case. On failure the test prints a compact table of the
parameter, expected value, fitted value, deviation in σ, and error ratio.

### Pass criteria

**Exact cases**, for each parameter with a known value `p*`:
- `|p_fit − p*| ≤ k · σ_fit + δ_abs`, with `k = 3` by default. `δ_abs` is a small per-case slack
  for known systematic bias, for example corrections to scaling at small `L`. It defaults to 0 and
  is justified in a comment wherever it is non-zero.
- `σ_fit ≤ σ_max`. This is a per-case ceiling, set at about 3× the error obtained with the fixed
  code. It catches a "grid point with huge error" failure even when the grid point happens to sit
  near the truth (θ = π/2 above: r_c = 0.5 is exact, but σ is 4× too large).

**Regression cases**, compared against `references.toml` (`optimal_ps`, `optimal_ps_error`,
`minimum`):
- `|p_fit − p_ref| ≤ 0.25 · σ_ref` for every parameter.
- `0.5 ≤ σ_fit / σ_ref ≤ 2`.
- These numbers are the default. They tolerate small numerical drift from dependency updates, and
  they fail on the Optim bug by a wide margin. They can be overridden per case with a comment.

`update_references.jl` reruns all regression cases and rewrites `references.toml`. The resulting
diff is what a reviewer looks at when a change is supposed to move the numbers.

**Published values are a separate, frozen anchor** *(Claude's call)*. Where a case comes from a
paper, `references.toml` holds two sections for it:
- `[published.<case>]`: the paper's numbers (value ± error). These are written by hand and never
  touched by `update_references.jl`.
- `[last.<case>]`: the last accepted output of the code, including `minimum`. This section is
  regenerable.

Paper cases are checked against both, with the criteria above. A case is only added once it has
been shown to reproduce its published numbers. A mismatch means either the fix is incomplete or
the case needs a written explanation; that explanation goes next to the case, and the case may
then check only `last`.

### Initial cases

Datasets marked "generate" are produced once by a script in `generators/` and committed.

**Exact or known answer**

| # | Case | Data | Known values | Notes |
|---|---|---|---|---|
| E1 | Qubit loss, θ = π/2 (2505.22720, Fig. 3) | reduced `I_c` curves, L = 32–256 | ν = 4/3, r_c = 1/2 | `:x`, MultipleSplines, the paper's call |
| E2 | Ising ↔ toric code, γ = 0.3 vertical cut (2504.12385, Fig. 3 inset) | `sweep_pd_vertical`, d = 8–64 | β_c'/β_c = 1, ν = 1 | `:xy`; y-exponent ≈ 0 (not asserted exactly). The published fit 1.00058(16) is 3.6σ from 1 (small-size bias), so this case needs `δ_abs ≈ 1e-3` on β_c'/β_c |
| E3 | 2D site percolation order parameter (autoScale.py, MIT) | 7 text files, L = 8–512 | p_c = 0.5927460, 1/ν = 3/4, β/ν = 5/48 | `:xy` with a negative y-exponent; also run with `:xny` |
| E4 | 2D bond percolation, Newman–Ziff | generate, L = 8–512 | p_c = 1/2, ν = 4/3; spanning prob. y-exp 0; P∞ y-exp −5/48 | Replaces the pyfssa file (unclear license) |
| E5 | Exact 2D Ising strip ξ_L/L (transfer matrix) | generate, noise-free | T_c = 2/ln(1+√2), ν = 1, y-exp 0 | Tight tolerance; tests the no-error-bars path |
| E6 | 2D Ising Binder (existing `IsingMonteCarlo.h5`) | existing | T_c = 2.269, ν = 1 | Move the existing checks under the new framework |
| E7 | Curie–Weiss, finite N exact | generate, noise-free | T_c = 1; in N: x-exp 1/2, Binder y-exp 0, ⟨m²⟩^½ y-exp −1/4 | Trivial and fast |
| E8 | Synthetic family | generate | chosen per variant | See below |

E8 variants (synthetic `y = L^a · F((x − x_c)/x_c · L^(1/ν)) · (1 + c·L^(−ω)) + noise`):
- clean;
- Gaussian noise with error bars;
- noise without error bars;
- corrections to scaling (ω = 1, BSA-style);
- only 3 sizes;
- uneven or non-overlapping x-grids per L;
- a monotonic F (crossing) and a peaked F (susceptibility-like);
- each preset (`:x`, `:xy`, `:xny`);
- one variant with a fixed parameter passed as a kwarg.

Each variant runs with `Linear`, `SingleSpline`, `MultipleSplines` and `Houdayer`. A quality
function may get its own tolerance when it is known to be weaker (as `SingleSpline` already is in
`test/ising.jl`).

**Regression (published estimates)**

| # | Case | Data | Reference |
|---|---|---|---|
| R1 | Qubit loss, all 15 θ (2505.22720, Fig. 3) | reduced `I_c` curves (θ = 0 includes L = 512) | published ν per θ (below); r_c published only for θ = 0, so the other r_c are `last` only |
| R2 | 2504.12385, Fig. 2 inset (β = β_c) | `betac.jld2`, d = 8–128 | γ_c = 0.59805(205), ν = 2.658(54), p3 = 0.0417(33) |
| R3 | 2504.12385, phase-diagram cuts (Fig. 4/9) | `sweep_pd_horizontal` / `_vertical` | stored fits `FSS_horizontal.jld2` / `FSS_vertical.jld2` |
| R4 | 2608.20268, Potts N⁽²⁾ and N⁽¹⁾ insets (Fig. 4) | `Potts_Ic_horizontal/3/averaged_results.jld2` | γ_c = 0.58067(129), ν = 2.476(33); γ_c = 0.76383(85), ν = 1.453(33) |
| R5 | 2608.20268, cuts incl. Fig. 10 (custom `AffineSimplexer(a=1e-6, b=1e-6)`) | `Potts_Ic_*/averaged_results.jld2` | stored `FSS.jld2` |
| R6 | Holstein CDW (existing `test/cdw.jl`) | existing | current result |

R1 published ν for θ index 1…15 (θ = 0 … π/2 in 15 equal steps), from the `fig3.ipynb` output.
All use `Ls[3:end]`, `:x`, `dx=[-1.25, 1.25]`, `p_space=[0.3:0.001:0.6, 1.2:0.1:1.8]`,
`error_threshold=1.3` and `MultipleSplines()`:
1.5316(40), 1.5349(44), 1.5320(48), 1.5328(52), 1.5325(60), 1.5298(75), 1.5287(71), 1.5295(66),
1.5267(57), 1.5231(57), 1.5215(51), 1.5122(49), 1.4886(70), 1.4206(83), 1.3370(59).
θ = 0 additionally has r_c = 0.336814(21).

**Provenance of R5 (verified).** The LearningPotts Zenodo `code/Manifest.toml` pins
ScalingCollapse **v0.2.5** (its tree hash matches the `v0.2.5` tag) with **Optim 1.13.3**. Re-running
the Fig. 10 notebook call with exactly that environment reproduces all 30 cuts × 3 parameters of
`Potts_Ic_vertical/4/FSS.jld2` exactly. The tiny `a = b = 1e-6` simplex never produces a
non-finite vertex, which is why these fits are correct despite Optim 1.13. Current `main` (with
#42) reproduces them to within 0.12σ, which is inside the 0.25σ regression tolerance.

Even on these vertical cuts, the fitted β_c'/β_c is about 1.0049 rather than the exact 1, which is
tens of σ away (all sizes from L = 8 are used). So R5 stays regression-only. Do not promote it to
an exact case.

For the other paper cases, the authoritative sources are each deposit's own notebook and
`Manifest.toml`. Each case records which package and Optim versions produced its published
numbers. A case whose published numbers came from Optim ≥ 1.12 with a large default simplex must
first be checked for the bug before its numbers are used as an anchor.

For R3 and R5 (about 140 cuts in total), start with all of them and measure the runtime. If the
suite exceeds the CI budget, keep an evenly spaced subset that includes each paper's critical
points, and document which cuts are included.

### Data conversion and provenance

- **Paper data** (Zenodo, CC-BY-4.0) is converted once into compact HDF5 fixtures: per size,
  `x`, `y`, `y_err`, plus `L` and attributes for the source, DOI and figure. The conversion scripts
  read the original JLD2 files. Those inputs are not in the repo: the user provides them locally
  (the Zenodo zips are 1–18 GB, but single members can be pulled by HTTP range requests). The
  scripts therefore live in `generators/` (or a `dev/` folder) and are documented, not run in CI.
  JLD2 and Measurements are not added to `test/Project.toml`.
- **Qubit loss:** the reduction from raw `expect_ZZs` samples to `I_c` mean ± SEM follows
  `fig3.ipynb` exactly: `Ic(κ)`, then the mean, then `std/√n`. An HDF5 file with all 15 θ, about
  185 KB, has already been produced and checked. Its layout is `theta{i}/L{L}/{r, Ic, Ic_err,
  nsamples}`, plus top-level `Ls` and `thetas`. It currently lives only in a temporary scratchpad,
  so the builder regenerates it with the reduction script, which takes about a minute. The script
  is committed to `generators/` and its input is the user's local `NishimoriQubitLoss/data`
  (`sweep_pd.jld2`, `nishimori_line.jld2`).
- `testdata/README.md` lists every file with its source, citation, DOI or URL, license, and the
  transformation applied.
- Total committed data target: under 20 MB.

### Licensing

- **The user's Zenodo deposits (CC-BY-4.0):** fine, with attribution in `testdata/README.md`.
- **autoScale.py data (MIT):** fine. Include the MIT notice next to the files.
- **pyfssa data (unclear license) and Kenji Harada's BSA data (no license):** not copied. Equivalent
  data is generated ourselves (E4, E6/E8).

### Runtime

Single fits on this data take 0.1–2.5 s. The budget allows a few hundred fits per CI job, including
the error analysis. The suite prints its total time. If it gets close to the budget, reduce the R3
and R5 subsets first, and never the exact cases.

## Order of work

1. Part 1 (the fix), as its own PR and 0.2.6 release, with the Optim 2 compat widened. Add E1 and
   the θ = 0 case from R1 as plain tests (the qubit-loss fixture is committed in this PR), and
   confirm they fail on the current `main` with Optim 1.13 and with 2.x before the fix. The suite
   framework later adopts these two tests.
2. The framework (case record, runner, pass criteria, references file and update script), with the
   existing Ising and CDW tests moved in.
3. Self-generated datasets (E4, E5, E7, E8).
4. Third-party data (E3).
5. The user's paper data (E2, R1–R5). This needs the local Zenodo files from the user.

## Left to the builder

- File names, the case-record representation (`NamedTuple` vs struct), and helper APIs.
- The exact mechanics of the finite-quality simplex shrink.
- Generator implementation details (lattice code, sample counts, grids), within the stated sizes
  and reference values.
- Whether the downgrade-compat CI job is added.

## Confirmed by the user

- E1: r_c = 1/2 exactly at θ = π/2 (bond percolation). The fit gives 0.499984(43).
