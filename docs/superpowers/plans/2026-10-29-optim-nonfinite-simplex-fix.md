# Optim ≥ 1.12 Nelder–Mead Fix (v0.2.6) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the Nelder–Mead step work on Optim ≥ 1.12 (including 2.x) by keeping every
initial-simplex vertex at a finite quality. Pin the fix with two real-data tests from
arXiv:2505.22720, widen the Optim compat to 2, and release it as v0.2.6.

**Architecture:** PR #42's `BoundedSimplexer` already reflects or clamps vertices into `p_space`.
It gains the objective `f`, and after the bounds step it halves any vertex with non-finite `f`
toward the starting point. `_bounded_optimizer(sp)` becomes `_bounded_optimizer(sp, f)`, so both
call sites pass the same closure they hand to `optimize`. A committed HDF5 fixture of the reduced
qubit-loss `I_c` curves feeds a new test file.

**Tech Stack:** Julia ≥ 1.9, Optim.jl (1.x and 2.x), HDF5.jl (already a test dependency), Test.

**Spec:** `docs/superpowers/specs/2026-10-29-collapse-benchmark-suite-design.md` (Part 1 only; Part 2,
the suite, gets its own plan).

## Global Constraints

- Branch from `origin/main` (currently `fbe327e`, which includes #42). Branch name:
  `fix-nonfinite-simplex`.
- Ship as its own PR and a 0.2.6 patch release. `Project.toml` `version = "0.2.6"`.
- Optim compat: `Optim = "1.6.2, 2"`. Do not pin Optim to 1.11.
- No new test dependencies. JLD2 and Measurements must not enter `test/Project.toml`.
- A user-supplied `optimizer.initial_simplex` remains the base simplex. #39's configuration
  option must keep working (`test/nelder_mead_configuration.jl` stays green, unmodified).
- Tests run with working directory `test/`. Fixture paths are relative to it.
- CI covers 32-bit Windows: no `Int64` literals or assumptions in new code.
- Never loosen the tolerances below to get a pass. If a pass criterion fails after the fix,
  investigate first.
- Each commit message ends with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.

## Review Focus

1. **Starting point on the edge of `p_space`.** This happens, for example, in the README Ising
   case, where the scan minimum sits at ν = 1.5, the upper bound. Reflection has to run before the
   halving: halving toward an edge point from outside never gets back inside. The existing README
   test in `test/ising.jl` must still pass.
2. **A vertex that never becomes finite within the halving limit.** The simplexer must terminate
   and return the vertex. `optimize_starting_ps!`'s existing non-convergence warning stays the
   user-facing signal. A test covers this in Task 2.
3. **Halving collapses a vertex onto the starting point**, which gives a degenerate simplex. The
   halving is capped, so a vertex gets at most 2⁻²⁰ of its step and never exactly x0. The unit test
   in Task 2 asserts all vertices are distinct.
4. **Custom simplexers** (the Fig. 10 `AffineSimplexer(a=1e-6, b=1e-6)` style and
   `RecordingSimplexer` in the tests) must be wrapped, not replaced. They are already covered by
   `test/nelder_mead_configuration.jl`, which must stay unmodified and green.
5. **The parameter-scan stage** uses `sp.quality_scan`, not `sp.quality`, as `f`. Passing the wrong
   closure would shrink against the wrong function. Task 2 makes both call sites pass their own
   closure.

---

### Task 1: Qubit-loss fixture and failing regression tests

**Files:**
- Create: `test/testdata/qubitloss/qubitloss_Ic.h5`
- Create: `test/generators/qubitloss_reduce.jl` (documentation of provenance; not run in CI)
- Create: `test/testdata/README.md`
- Create: `test/qubitloss.jl`
- Modify: `test/runtests.jl` (add `include("qubitloss.jl")` after `include("cdw.jl")`)

**Interfaces:**
- Produces: the fixture layout `theta{i}/L{L}/{r, Ic, Ic_err, nsamples}` (Float64 vectors;
  `nsamples` Int), with top-level `Ls = [8, 16, 32, 64, 128, 256]` and `thetas` (15 values,
  0 … π/2). `theta1/L512` also exists. Part 2 of the spec reuses this file unchanged.

- [ ] **Step 1: Put the fixture and reduction script in place.**
  - The verified file and its script exist at
    `/private/tmp/claude-501/-Users-malte-github-ScalingCollapse-jl/4aae5966-721c-4e31-ac8f-b6598c9c0ad2/scratchpad/qubitloss/qubitloss_Ic.h5`
    (185281 bytes) and `.../scratchpad/qubitloss/reduce.jl`. Copy them to the two paths above.
  - Add a header comment to the script: purpose, inputs (`sweep_pd.jld2` and `nishimori_line.jld2`
    from Zenodo 10.5281/zenodo.15530142, `NishimoriQubitLoss.zip`), and usage
    `julia --project=<env with JLD2, HDF5> qubitloss_reduce.jl <data dir> <out.h5>`.
  - If the scratchpad copy is gone, ask the user to copy `NishimoriQubitLoss/data` to a folder
    outside `~/Desktop` (macOS blocks Desktop access) and regenerate. The reduction follows
    `fig3.ipynb` exactly:
    - `Ic(κ) = -(1+κ)/2·log2(|1+κ|/2) - (1-κ)/2·log2(|1-κ|/2)`, retried in `BigFloat` when
      non-finite, with `NaN` mapped to 0;
    - per `(θ, L, r)`: the mean of `Ic.(expect_ZZs)`, and `std/√n` as the error.

- [ ] **Step 2: Write `test/testdata/README.md`.** One entry per file:
  - `IsingMonteCarlo.h5`: existing; source unknown, so leave a TODO for the maintainer.
  - `qubitloss/qubitloss_Ic.h5`:
    - Pütz et al., "Flow to Nishimori universality in weakly monitored quantum circuits with qubit
      loss", arXiv:2505.22720;
    - data from Zenodo 10.5281/zenodo.15530142, CC-BY-4.0;
    - transformation: raw `expect_ZZs` → `I_c` mean ± SEM via `test/generators/qubitloss_reduce.jl`.

- [ ] **Step 3: Write the failing tests in `test/qubitloss.jl`.**
  - Add a loader `load_qubitloss(theta_index, Ls) -> (xs, ys, es, Ls)` that reads the fixture.
  - Add a fit helper that calls `ScalingProblem(xs, ys, es, Ls; ...)` with the paper's kwargs:
    `sf=ScalingFunction(:x)`, `dx=[-1.25, 1.25]`, `p_space=[0.3:0.001:0.6, 1.2:0.1:1.8]`,
    `error_threshold=1.3`, `quality=MultipleSplines()`.

```julia
@testset "Qubit loss (arXiv:2505.22720, Fig. 3)" begin
    @testset "θ = π/2: bond percolation (exact)" begin
        sp = fit_qubitloss(15, [32, 64, 128, 256])
        (rc, nu), (σrc, σnu) = sp.optimal_ps, sp.optimal_ps_error
        @test abs(rc - 0.5) <= 3σrc
        @test abs(nu - 4 / 3) <= 3σnu
        # ceilings at ~3× the error of a converged fit (4.33e-5, 5.92e-3); a scan grid point
        # returned without Nelder-Mead progress has 1.9e-4 and 8.8e-2
        @test σrc <= 1.3e-4
        @test σnu <= 1.8e-2
    end
    @testset "θ = 0: Nishimori (published, v0.2.2 + Optim 1.11)" begin
        sp = fit_qubitloss(1, [32, 64, 128, 256, 512])
        ref = [0.33681399873460816, 1.5316126320454644]
        σref = [2.0879977065466418e-5, 0.004006908264828457]
        @test all(abs.(sp.optimal_ps .- ref) .<= 0.25 .* σref)
        @test all(0.5 .<= sp.optimal_ps_error ./ σref .<= 2)
    end
end
```

- [ ] **Step 4: Run the tests on unfixed `main` and confirm they fail.**
  - Run: `cd test && julia --project=.. -e 'using Pkg; Pkg.test()'`. Optim resolves to 1.13.x; check
    this with `Pkg.status("Optim"; mode=PKGMODE_MANIFEST)` in the test environment.
  - Expected: FAIL. The θ = π/2 σ ceilings fail (σ about 1.9e-4 and 8.8e-2), and the θ = 0
    deviation fails (0.337 vs 0.336814). All other test sets pass.
  - Record the actual numbers for the PR description.
  - The spec also requires these tests to fail before the fix on Optim 2.x. That was verified on
    2026-10-29 (`main` with compat widened, Optim 2.3.2: θ = 0 gave `[0.337, 1.5]`, θ = π/2 gave
    ν = 1.3775 ± 0.088). Quote it in the PR; no need to re-run.

- [ ] **Step 5: Commit** (local only; do not push a red branch).
  - Stage `test/testdata/qubitloss/qubitloss_Ic.h5`, `test/generators/qubitloss_reduce.jl`,
    `test/testdata/README.md`, `test/qubitloss.jl` and `test/runtests.jl`.
  - Message: `Add qubit-loss collapse regression tests (fail on Optim >= 1.12)`.

### Task 2: Shrink non-finite simplex vertices toward the start

**Files:**
- Modify: `src/optimization.jl`: `BoundedSimplexer`, `Optim.simplexer(::BoundedSimplexer, x)`,
  `_bounded_optimizer`, the two `optimize(...)` call sites, and the non-convergence hint text
- Test: `test/nelder_mead_configuration.jl` (append a new testset only; existing tests unchanged)

**Interfaces:**
- Consumes: the Task 1 tests as the acceptance check.
- Produces:
  - `BoundedSimplexer(simplexer::Optim.Simplexer, lower::Vector{Float64}, upper::Vector{Float64}, f)`,
    which is internal and not exported;
  - `_bounded_optimizer(sp, f) -> Optim.NelderMead`;
  - `const _MAX_SIMPLEX_HALVINGS = 20`.

- [ ] **Step 1: Write the failing unit test** by appending to `test/nelder_mead_configuration.jl`:

```julia
@testset "initial simplex avoids non-finite quality" begin
    # the AffineSimplexer vertex 0.6125 leaves p_space and is reflected to 0.3875, which is
    # inside p_space but has Inf quality -> must be halved toward 0.5
    f(x) = x[1] < 0.45 ? Inf : sum(abs2, x .- [0.5, 1.4])
    S = ScalingCollapse.BoundedSimplexer(
        ScalingCollapse.Optim.AffineSimplexer(a=0.025, b=0.1), [0.3, 1.2], [0.6, 1.8], f)
    simplex = ScalingCollapse.Optim.simplexer(S, [0.5, 1.4])
    @test all(v -> isfinite(f(v)), simplex)
    @test all(v -> all([0.3, 1.2] .<= v .<= [0.6, 1.8]), simplex)
    @test length(unique(simplex)) == 3

    # never finite away from x0: terminates, keeps vertices distinct
    g(x) = x == [0.5, 1.4] ? 0.0 : Inf
    S = ScalingCollapse.BoundedSimplexer(
        ScalingCollapse.Optim.AffineSimplexer(a=0.025, b=0.1), [0.3, 1.2], [0.6, 1.8], g)
    simplex = ScalingCollapse.Optim.simplexer(S, [0.5, 1.4])
    @test length(unique(simplex)) == 3
end
```

- [ ] **Step 2: Run it and confirm it fails.**
  - Run: `cd test && julia --project=.. -e 'using Pkg; Pkg.test()'`.
  - Expected: FAIL with a `MethodError` from the 4-argument `BoundedSimplexer` constructor.

- [ ] **Step 3: Implement.**
  - Add the field `f` to `BoundedSimplexer` (type-parametrize it as `F`).
  - In `Optim.simplexer(S::BoundedSimplexer, initial_x)`, keep the existing reflect/clamp loop
    unchanged. Then, for each vertex, while `!isfinite(S.f(vertex))` and the halving count is below
    `_MAX_SIMPLEX_HALVINGS`, set `vertex .= initial_x .+ (vertex .- initial_x) ./ 2`. Skip vertices
    equal to `initial_x`.
  - Change `_bounded_optimizer(sp)` to `_bounded_optimizer(sp, f)`.
  - At each call site, bind the closure once and pass it both to `_bounded_optimizer` and to
    `optimize`. In `parameter_scan!` that is
    `f = ps -> _nan_to_inf(sp.quality_scan(sp, ps; check_bounds=true))`; in
    `optimize_starting_ps!` the same with `sp.quality`.
  - Update the comment block above `BoundedSimplexer` to say it handles both out-of-bounds and
    non-finite-quality vertices, and mention that Optim ≥ 1.12 stops at iteration 0 otherwise.
  - In the non-convergence hint, replace "reaches outside of p_space" with "has vertices where the
    quality function is not finite (e.g. no overlap of the rescaled data within dx)". Keep the
    `AffineSimplexer(b=0.01)` suggestion.

- [ ] **Step 4: Run the full suite and confirm it passes.**
  - Run: `cd test && julia --project=.. -e 'using Pkg; Pkg.test()'` (Optim 1.13.x).
  - Expected: all tests pass, including both qubit-loss test sets and the README case in
    `test/ising.jl`.
  - If θ = 0 misses the 0.25σ criterion: the published run used a different base simplex
    (v0.2.2 hardcoded `AffineSimplexer(b=0.1)`, now `a=0.025, b=0.1`). Report the observed
    deviation in σ to the maintainer instead of widening the tolerance.

- [ ] **Step 5: Commit.**
  - Stage `src/optimization.jl` and `test/nelder_mead_configuration.jl`.
  - Message: `Shrink non-finite initial simplex vertices toward the start point`.

### Task 3: Optim 2 compat, version 0.2.6, cross-version verification, PR

**Files:**
- Modify: `Project.toml` (`version = "0.2.6"`; `[compat] Optim = "1.6.2, 2"`)

- [ ] **Step 1: Edit `Project.toml`** as above.

- [ ] **Step 2: Run the suite on Optim 2.x.**
  - Run: `cd test && julia --project=.. -e 'using Pkg; Pkg.test()'`.
  - Expected: Optim resolves to 2.3.x (check this as in Task 1), and all tests pass.

- [ ] **Step 3: Run the suite on Optim 1.13.x and 1.11.0.**
  - Use a scratch environment that `Pkg.develop`s this checkout plus Aqua, HDF5 and Test, and pins
    Optim with `Pkg.add(name="Optim", version="1.13.3")`, then `"1.11.0"`.
  - Run `test/runtests.jl` from `test/` in each.
  - Expected: all tests pass on both versions.

- [ ] **Step 4: Commit.**
  - Stage `Project.toml`.
  - Message: `Allow Optim 2 and bump version to 0.2.6`.

- [ ] **Step 5: Push and open the PR against `main`.**
  - Title: `Fix Nelder-Mead stalling on Optim >= 1.12; allow Optim 2`.
  - The body includes:
    - the before/after table (θ = 0 and θ = π/2 on Optim 1.11 / 1.13 / 2.3, from Task 1 Step 4
      and the Task 3 runs);
    - the root cause in two sentences;
    - the note that the existing suite passed with the bug present;
    - the data attribution;
    - a note that registering v0.2.6 after merge is the maintainer's step.

  Then bind the PR with the ccd_pr tools.
