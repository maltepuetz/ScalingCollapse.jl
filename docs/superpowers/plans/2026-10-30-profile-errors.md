# Profile Parameter Errors Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the default parameter error a local profile width (other free parameters re-fitted), keep today's method as `error=:conditional`, and recalibrate the collapse suite.

**Architecture:** `error_analysis!` dispatches on a normalized `sp.error::Symbol`. `:conditional` is today's code moved unchanged into `_conditional_errors!`. `:profile` first computes the conditional roots, then brackets outward from them on a profile quality that re-fits the nuisance parameters with Nelder–Mead. A lower quality found on the way triggers re-optimization and a restart of the whole analysis.

**Tech Stack:** Julia ≥ 1.9, Optim.jl (Nelder–Mead), Test, HDF5/TOML (test only).

**Spec:** `docs/superpowers/specs/2026-10-30-profile-errors-design.md`. Read it first. The plan argues from it and does not repeat its rationale.

## Global Constraints

- Work on a new branch `profile-errors` from `origin/main`, in a git worktree. The main checkout has untracked copies of tracked files (`CLAUDE.md`, `docs/superpowers/`), so don't switch branches there. First commit: the spec and this plan.
- No version bump in this PR. `Project.toml` `version` stays unchanged.
- `:conditional` uses today's code path, bracket start and precision 1e-10, unchanged.
- Default `error_threshold` stays 1.3.
- CI runs on Julia `lts`, `1` and `nightly`, on Linux, macOS (aarch64) and Windows, including 32-bit Windows. Don't write `Int64` where `Int` is meant.
- The runtime increase of the suite (~100 s → ~500 s) is accepted.
- Warnings go through `@warn`. Their messages contain the regex anchors the tests use (listed per task).
- Run tests with `julia --project=. -e 'using Pkg; Pkg.test()'`. For single files, use the TestEnv REPL route in `CLAUDE.md`, with the working directory `test/`.

## Review Focus

1. **`starting_ps` given:** `p_space` is the sentinel `[-1_000_000:1_000_000 ...]`, so nuisance bounds are effectively absent. Profile errors must still be finite and ≥ conditional. Covered in Task 3 (`profile with starting_ps`).
2. **Fixed parameters in the `ScalingFunction`:** the profile runs over free parameters only. A fit with two free parameters and one fixed must work. One free parameter must take the identity shortcut. Covered in Task 3 (`fixed parameter`) and Task 3 (`one free parameter`).
3. **Per-parameter `error_threshold` vectors:** each parameter uses its own threshold. Covered in Task 3 (`per-parameter thresholds`).
4. **Nuisance box truncates the valley:** the width lies between the conditional width and the unbounded profile width, and equals the analytic clamped value. Covered in Task 3 (`nuisance bounds`).
5. **Optimum at exactly zero:** today's relative step `0.1·|p|` becomes absolute 0.1. The profile path must not loop or divide by zero there. Covered in Task 3 (`optimum at zero`).

---

### Task 1: Freeze today's errors as `:conditional` fixtures

**Files:**
- Create: `test/error_methods.jl`
- Modify: `test/runtests.jl` (add `include("error_methods.jl")` after `include("error_handling.jl")`)

**Interfaces:**
- Produces: `test/error_methods.jl` with a testset `"error methods"`. Later tasks add testsets to this file.

- [ ] **Step 1: Record today's errors** on `origin/main`, before any `src/` change. Use a scratch script, not committed, that runs these four fits with the default `error=true` and prints `repr.(sp.optimal_ps_error)`:
  - (a) the README Ising fit in `test/ising.jl:37-41` (`p_space=[2.0:0.1:2.4, 0.5:0.1:1.5]`, `dx=[-1.0, 1.0]`, `sf=ScalingFunction(:x; p_names=["T_c","nu"])`, data from `testdata/IsingMonteCarlo.h5`);
  - (b) the second CDW fit in `test/cdw.jl` (`:xny`, `nu=1` fixed, gamma free);
  - (c) `bond_strength` via `load_percolation_h5("strength")` and `BOND_STRENGTH_KWARGS` from `test/collapses/cases_percolation.jl`;
  - (d) `synthetic_clean_xy_Linear` from `test/collapses/cases_synthetic.jl`.

- [ ] **Step 2: Write the failing test** `"conditional errors unchanged"` in `test/error_methods.jl`. It reruns the four fits with `error=:conditional` and asserts `isapprox(sp.optimal_ps_error, RECORDED; rtol=1e-8)`, with the recorded literals.

  *Deviation from the spec's `==`:* the literals come from one machine. CI runs on x86-64, aarch64 and 32-bit, so last-bit floating-point differences between platforms are possible. `rtol=1e-8` still catches any change in method or precision. Same-machine identity is checked once in Task 2 Step 5.

- [ ] **Step 3: Run it.** Expected: FAIL (`error=:conditional` is not accepted yet, e.g. a `TypeError` converting `Symbol` to `Bool`).

- [ ] **Step 4: Commit** `test/error_methods.jl` and `test/runtests.jl` together with the spec and plan: "Freeze today's parameter errors as conditional-mode fixtures".

### Task 2: `error` keyword accepts a method; `:conditional` path

**Files:**
- Modify: `src/scalingproblem.jl`: field `error::Bool` → `error::Symbol` (line ~133); kwarg parsing (line ~222); `solve!` (line ~398); docstring of `error` (lines ~24-25).
- Modify: `src/error_analysis.jl`
- Modify: `test/error_handling.jl` (the `harness` NamedTuple gains `error=:conditional`)
- Test: `test/error_methods.jl`

**Interfaces:**
- Produces: `_error_method(x)::Symbol` in `src/scalingproblem.jl`. It maps `true`/`:profile` → `:profile`, `:conditional` → `:conditional`, `false` → `:none`, and throws `ArgumentError` otherwise. The message lists `true, false, :profile, :conditional`.
- Produces: `error_analysis!(sp; max_restarts::Int=3)`, which dispatches on `sp.error` (`:conditional` → `_conditional_errors!(sp)`, `:profile` → `_profile_errors!(sp; max_restarts)`).
- Produces: `_conditional_errors!(sp)`: today's `error_analysis!` body, moved verbatim.
- `solve!` calls `error_analysis!(sp)` unless `sp.error === :none`.

- [ ] **Step 1: Write failing tests** in `test/error_methods.jl`, testset `"error keyword"`:
  - `@test_throws ArgumentError ScalingProblem(Ts, binder, Ls; sf=ScalingFunction(:x), dx=[-1.0, 1.0], p_space=[2.0:0.1:2.4, 0.5:0.1:1.5], error=:bogus)`
  - `error=false` gives `sp.error === :none` and `all(iszero, sp.optimal_ps_error)`.
  - `error=true` gives `sp.error === :profile`, and `error=:conditional` gives `sp.error === :conditional`.

- [ ] **Step 2: Run them.** Expected: FAIL.

- [ ] **Step 3: Implement** `_error_method`, the field type change, `solve!`, the dispatch and `_conditional_errors!`. Until Task 3 lands, `_profile_errors!(sp; max_restarts)` calls `_conditional_errors!(sp)` and has a `# TODO Task 3` comment. Add `error=:conditional` to the `harness` in `test/error_handling.jl` so its existing assertions keep testing today's method.

- [ ] **Step 4: Run the full test suite.** Expected: PASS, including Task 1's fixtures.

- [ ] **Step 5: Same-machine identity check (one-off).** On the machine that recorded Task 1, assert `==` instead of `rtol` in a REPL. Expected: all four equal.

- [ ] **Step 6: Commit:** "Accept error=:profile/:conditional/true/false".

### Task 3: Profile width

**Files:**
- Modify: `src/error_analysis.jl`
- Test: `test/error_methods.jl` (testset `"profile width"`)

**Interfaces:**
- Consumes: `_conditional_errors!`; the helpers `_bracket`, `_root`, `_one_vector`; `BoundedSimplexer` and `_nan_to_inf` (`src/optimization.jl`).
- Produces:
  - `_profile_errors!(sp; max_restarts::Int=3)`. In this task it has no restart handling. A detected lower quality is just recorded and ignored until Task 4.
  - `_profile_value(sp, i::Int, x::Float64, warm::Vector{Float64}) -> (q::Float64, nuisance::Vector{Float64}, converged::Bool)` implements spec "Profile quality" steps 1–5. `q` is `Inf` for invalid points; NaN is mapped to Inf.
  - `_profile_side(sp, i::Int, d::Float64, f0::Float64) -> Union{Float64,Nothing}` implements spec "Bracketing and bisection" steps 2–5 for one side, given conditional root `d`. It returns the width or `nothing` (→ NaN).
- Constants: growth factor `1.5`; bisection relative precision `1e-3`; at most `60` outward steps.

Test harness: a `Base.@kwdef mutable struct ErrorHarness` in `test/error_methods.jl`, with the fields `error_analysis!` and `optimize_starting_ps!` read:
- `quality`: called as `(sp, ps; check_bounds=false)`; the test qualities ignore `check_bounds` unless stated;
- `optimal_ps`, `optimal_ps_error`, `minimum`, `error_threshold`, `p_space`, `starting_ps`, `verbose=false`, `error=:profile`, `sf=ScalingFunction(:x; p_names=["a","b"])`;
- `optimizer=Optim.NelderMead()`, `optimizer_options=Optim.Options()`.

`p_space` defaults to `[-10.0:0.1:10.0, -10.0:0.1:10.0]`.

- [ ] **Step 1: Write failing tests** (threshold 1.3 and minimum 1.0 unless stated):
  - `correlated quadratic`: `q = 1 + a² + b² + 2ρab`, ρ = 0.8, optimum `[0,0]`. Expected: both errors `≈ sqrt(0.3/(1-0.64))` (= 0.912871), `rtol=1e-3`. Under `error=:conditional` both are `≈ sqrt(0.3)`.
  - `one free parameter`: the Ising fit with `sf=ScalingFunction(:x; nu=1)` gives `==` errors under `:profile` and `:conditional`.
  - `two valleys`: `q = 1 + b² + 25a²(a−1)²`, optimum `[0,0]`. Expected: `error[1] ≈ 0.125224` (`rtol=1e-3`), not ≈ 1.1.
  - `invalid region`: `q = abs(a) < 0.5 ? 1 + a² + b² : Inf`. Expected: `error[1]` is NaN, with `@test_logs (:warn, r"undefined")`; `error[2] ≈ sqrt(0.3)`. With the cutoff at `0.6` instead: `error[1] ≈ sqrt(0.3)`, `rtol=1e-3`.
  - `nuisance bounds`: the correlated quadratic with ρ = 0.8 and `p_space = [-10.0:0.1:10.0, -0.1:0.01:0.1]`. Expected: `error[1] ≈ 0.62443` (root of `a² − 0.16a − 0.29 = 0`), `rtol=1e-3`. The quality callable must return `Inf` if `check_bounds` is ever `true`; this checks that the profile never asks the quality to bound-check.
  - `per-parameter thresholds`: `q = 1 + a² + b²`, `error_threshold=[1.3, 2.0]`. Expected: `≈ [sqrt(0.3), 1.0]`.
  - `optimum at zero`: `q = 1 + a² + b²` with optimum `[0,0]` finishes and gives `≈ sqrt(0.3)` for both.
  - `profile ≥ conditional`: for the README Ising fit (Task 1 (a)), every profile error ≥ the conditional error × (1 − 1e-3).
  - `profile with starting_ps`: the same fit with `starting_ps=[2.27, 1.0]` instead of `p_space`. Expected: finite errors, each ≥ conditional × (1 − 1e-3).
  - `fixed parameter`: the CDW fit of Task 1 (b), but with gamma free and nu fixed, gives finite profile errors ≥ conditional × (1 − 1e-3).

- [ ] **Step 2: Run them.** Expected: the profile-specific ones FAIL (they get conditional values).

- [ ] **Step 3: Implement** `_profile_errors!`, `_profile_value` and `_profile_side` following the spec's "One free parameter", "Bounds", "Profile quality" and "Bracketing and bisection" sections.
  - The nuisance objective is `ps -> any(out of p_space bounds for nuisance coords) ? Inf : _nan_to_inf(sp.quality(sp, full_ps; check_bounds=false))`.
  - The inner optimizer is `Optim.NelderMead(BoundedSimplexer(sp.optimizer.initial_simplex, lower, upper, f), sp.optimizer.parameters)` with `sp.optimizer_options`. Here `lower` and `upper` are the nuisance bounds.
  - Warm start: keep per side a list of `(x, nuisance)` and use the entry with the closest `x`.
  - Inner failure: `!Optim.converged(result)` → one retry from the conditional candidate. If it still fails, count it, and once per parameter emit `@warn` matching `r"may be underestimated"`.
  - Lower-quality detection: track the lowest `q` together with its full parameter vector. Use it in Task 4.

- [ ] **Step 4: Run** `test/error_methods.jl` and `test/error_handling.jl`. Expected: PASS.

- [ ] **Step 5: Commit:** "Compute profile parameter errors by default".

### Task 4: Lower quality found: restart, outside `p_space`, exhaustion

**Files:**
- Modify: `src/error_analysis.jl`
- Test: `test/error_methods.jl` (testset `"lower quality during the profile"`)

**Interfaces:**
- Consumes: `optimize_starting_ps!(sp)` (`src/optimization.jl`), which reads `sp.starting_ps` and sets `sp.optimal_ps` and `sp.minimum`.
- Produces: the final `_profile_errors!(sp; max_restarts)` behavior from the spec section "Lower quality found during the profile". The relative tolerance is `1e-6`.

- [ ] **Step 1: Write failing tests** with `ErrorHarness`:
  - `restart`: a narrow tilted valley `q = 1 + 100(a−b)² + 0.01(b−0.6)²`, reported `optimal_ps = [0.0, 0.0]`, `minimum = q(0,0) = 1.0036`.
    - Why it triggers: the conditional root of `a` is ≈ 0.055. Re-fitting `b` there gives `q ≈ 1.0030`, below the reported minimum, so the restart fires. (A plain isotropic bowl would not trigger it, because profile points are only evaluated beyond the conditional roots.)
    - Expected: a `@warn` matching `r"lower quality"`, then `sp.minimum ≈ 1.0` (`atol=1e-4`), `sp.optimal_ps ≈ [0.6, 0.6]` (`atol=1e-2`), and both errors `≈ sqrt(0.3 * 100.01)` (= 5.4775) with `rtol=1e-2`.
  - `outside p_space`: `q = 1 + (a−b)² + (abs(b) ≥ 0.5 ? -0.01 : 0.0)`, with `p_space = [-0.12:0.01:0.12, -1.0:0.01:1.0]`, optimum `[0,0]` and minimum 1. Expected: `error[1]` NaN, with `@warn` matching `r"outside p_space"`; `error[2]` finite; `sp.optimal_ps == [0.0, 0.0]` and `sp.minimum == 1.0` (no restart).
  - `restarts exhausted`: the `restart` setup with `error_analysis!(sp; max_restarts=0)`. Expected: all errors NaN, with `@warn` matching `r"could not be confirmed"`; `sp.minimum < 1.0036`, and `sp.optimal_ps` is inside `p_space`.

- [ ] **Step 2: Run them.** Expected: FAIL.

- [ ] **Step 3: Implement** the three branches.
  - A restart clears every per-analysis cache simply by calling `_profile_errors!` again from the top. No state survives across calls.
  - When restarts are exhausted, `optimal_ps` and `minimum` are set to the best feasible point seen, and all errors are NaN.

- [ ] **Step 4: Run** `test/error_methods.jl` and `test/error_handling.jl`. Expected: PASS.

- [ ] **Step 5: Commit:** "Re-optimize when the profile finds a lower quality".

### Task 5: Docs

**Files:**
- Modify: `src/scalingproblem.jl` docstring (the `error` and `error_threshold` entries, lines ~24-28)
- Modify: `README.md:37`
- Modify: `CLAUDE.md` (Architecture step 4)

- [ ] **Step 1: Update the docstrings and README.**
  - Document the four `error` values and the local profile definition.
  - State that the analysis may update `optimal_ps` and `minimum`.
  - Add the spec's "what these numbers are and are not" paragraph from its Goal section, in your own words.
  - Make `error_threshold`'s description say "with the other free parameters re-fitted (`:profile`) or fixed (`:conditional`)".
  - Check the README and docstrings for any stale "2×" threshold wording and fix it to 1.3.

- [ ] **Step 2: Update `CLAUDE.md` step 4** to describe `:profile` (the default) and `:conditional`. Drop its "not the 2×" remark if the README no longer says 2×.

- [ ] **Step 3: Build the docs.** Run `julia --project=docs docs/make.jl`. Expected: it builds without docstring errors.

- [ ] **Step 4: Commit:** "Document profile errors".

### Task 6: Collapse suite under profile errors

**Files:**
- Modify: `test/collapses/framework.jl` (`check_regression`, `run_cases`)
- Modify: `test/collapses/framework_tests.jl`
- Modify: `test/collapses/cases_percolation.jl` (`bond_strength_exact`)
- Modify: `test/collapses/cases_*.jl` (`σmax` of exact cases, and `SYNTHETIC_SIGMAMAX` in `cases_synthetic.jl`)
- Modify: `test/collapses/references.toml` (`[last]` regenerated; `[published]` entries removed only per the spec rule)

**Interfaces:**
- Produces: `check_regression(name, sp, ref; check_errors::Bool=true)`. `run_cases` passes `check_errors=false` for the `[published]` group.

- [ ] **Step 1: Write the failing framework test** in `framework_tests.jl`: with `check_errors=false`, the case `fake([1.0, 2.0], [0.21, 0.1])` (σ ratio 2.1) passes, and the value-off case `fake([1.03, 2.0], [0.1, 0.1])` still fails.

- [ ] **Step 2: Run it.** Expected: FAIL. Then implement the keyword: with `check_errors=false`, skip the `0.5 ≤ ratio ≤ 2` note. Run again. Expected: PASS.

- [ ] **Step 3: Switch `bond_strength_exact` to L = 64–256.** Its `load` slices `load_percolation_h5("strength")` from the index of `L == 64`. Remove `broken=true` and leave `δabs` empty. Rewrite the case comment with the size-window table and cause from the spec. The twin `bond_strength` keeps L = 16–256.

- [ ] **Step 4: Measure all exact cases.** Use a scratch script (not committed) that runs every `kind === :exact` case and prints each free parameter's fit, σ and |p − p*|/σ. For each case set `σmax[i]` = 3 × the new σ, rounded to 2 significant figures, which is the PR #44 rule. Update the "measured:" numbers in each case comment. Change `δabs` only if a parameter now fails at 3σ + δabs, and re-measure and comment it if you do.

- [ ] **Step 5: Regenerate `[last]`.** Run `julia -e 'using Pkg; Pkg.activate(; temp=true); Pkg.develop(path="."); Pkg.add(["HDF5", "TOML"]); include("test/collapses/update_references.jl")'`. Expected: only `[last]` changes in `git diff test/collapses/references.toml`.

- [ ] **Step 6: Run the collapse suite.** For every case that fails its `[published]` value check, apply the spec rule:
  - Confirm that the new `minimum` is lower than the published fit's quality. If it isn't, stop and report; that is not this rule's case.
  - Remove the `[published]` entry.
  - Add a comment to the case: old and new quality, the shift in σ, and the wording "a lower quality was found".
  - Collect the list for the PR description.

- [ ] **Step 7: Run the full test suite.** Expected: PASS on all non-broken cases; `bond_strength_exact` is a plain pass.

- [ ] **Step 8: Commit:** "Recalibrate the collapse suite for profile errors".

### Task 7: One-off checks and PR notes

**Files:**
- Create: `docs/superpowers/notes/2026-10-30-profile-errors-pr.md` (PR description draft)

- [ ] **Step 1: Profile ≥ conditional check.** Use a scratch script that runs every case under both methods. Expected: no profile error below conditional × (1 − 1e-3). Record the ratio quantiles.

- [ ] **Step 2: Stability check.** Rerun the profile errors with `optimizer_options=Optim.Options(g_tol=1e-12)` and the bisection precision temporarily at `1e-5` (local edit, reverted afterwards). Record the largest relative change in any error.

- [ ] **Step 3: E3 check.** Expected: `autoscale_site_order_exact` σ(p_c) is within 20 % of 7.5e-5.

- [ ] **Step 4: Write the PR description draft.** Include:
  - a summary of the method;
  - the API table;
  - the results of Steps 1–3;
  - removed `[published]` entries;
  - the `bond_strength_exact` change;
  - the suite runtime;
  - the note that the version bump follows separately.

  End it with `🤖 Generated with [Claude Code](https://claude.com/claude-code)`.

- [ ] **Step 5: Commit** the notes: "Add PR notes for profile errors". Don't push or open a PR unless the maintainer asks.
