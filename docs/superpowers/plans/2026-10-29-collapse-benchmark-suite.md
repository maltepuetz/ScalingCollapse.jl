# Collapse Benchmark Suite Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Turn the test suite into a broad collapse benchmark. Every case is checked against its
exact answer where one is known, and against stored (and, for papers, published) results
otherwise.

**Architecture:** A small framework lives in `test/collapses/`:
- a `CollapseCase` record;
- one runner;
- two pure checker functions that return result rows;
- `references.toml`, with a frozen `[published]` section and a regenerable `[last]` section;
- an explicit `update_references.jl`.

Case definitions are grouped by source. Data is committed under `test/testdata/<source>/`.
Generators and converters live in `test/generators/`, are run once by hand, and their outputs are
committed.

**Tech Stack:** Julia ≥ 1.9, Test, HDF5, TOML (stdlib). The generators also use JLD2 and
Measurements in their own scratch environment, never in `test/Project.toml`.

**Spec:** `docs/superpowers/specs/2026-10-29-collapse-benchmark-suite-design.md` (Part 2. Part 1
shipped as maltepuetz/ScalingCollapse.jl#43).

## Global Constraints

- **Branch:** `collapse-benchmark-suite`, cut from `fix-nonfinite-simplex`, because it needs #43's
  fix and fixture. If #43 is merged when the PR opens, rebase onto `main` and target `main`.
  Otherwise open the PR stacked on `fix-nonfinite-simplex`.
- **Dependencies:** no new non-stdlib test dependencies. Adding the stdlib `TOML`
  (uuid `fa267f1f-6049-4f14-aa54-33bafae1ed76`) to `test/Project.toml` is allowed.
- **Data:** everything is committed, and nothing is downloaded at test time. Generators use fixed
  seeds (`Random.Xoshiro(20261029)`). Total new data stays under 20 MB.
- **CI time:** under 10–15 min per CI job. Locally, the whole `Pkg.test()` must stay under 4 min,
  because 32-bit Windows CI is ~3× slower. If it goes over, reduce R3/R5 cuts first (evenly
  spaced, always keeping each paper's critical-point cut), never the exact cases.
- **Exact criterion:** `|p − p*| ≤ 3σ_fit + δ_abs` (δ_abs defaults to 0), and `σ_fit ≤ σ_max`
  with `σ_max ≈ 3 ×` the error of the fixed code's fit.
- **Regression criterion:** `|p − p_ref| ≤ 0.25 σ_ref`, and `0.5 ≤ σ_fit/σ_ref ≤ 2`, against
  both `[published]` (where it exists) and `[last]`.
- **`δ_abs`:** non-zero only for data with known finite-size bias (E2, E3, E4, E6). Each comes
  with a comment giving the reason and the measured deviation. Never for E1, E5, E7 or the clean
  E8 variants. Never loosen a criterion to get a pass without a ledgered ruling.
- `update_references.jl` never writes `[published]`.
- **Commits** end with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
- **Tests** run with working directory `test/`. CI includes 32-bit Windows, so no `Int64`
  literals or assumptions.
- **Licensing:** CC-BY data carries attribution in `test/testdata/README.md`. autoScale data
  (MIT) ships with its license file. No pyfssa or BSA files.

## Facts established while planning (2026-10-29, fixed code, Optim 2.3.2)

All rows are compared with the stored fits from the deposits.

| Set | Source call | Reproduces stored fits? |
|---|---|---|
| R2 2504 Fig. 2 inset | `ScalingProblem(gammas, Ics, ds; sf=:xy, dx=[-2,2], p_space=[0.3:0.1:0.9, 1.0:0.1:5.0, 0.0:0.01:0.1], quality=MultipleSplines())`, `betac.jld2`, ds = 8–128 | exactly |
| E2 + R3v 2504 vertical cuts (31) | `ScalingProblem(betas[2:5] ./ βc, [Ics[k][gi, :] for k in 2:5], ds[2:5]; sf=:xy, p_space=[0.6:0.01:1.3, 0.5:0.1:1.5, -0.5:0.1:0.5], quality=MultipleSplines(N_steps=100), dx=[-1.25, 1.75])`, βc = ln(1+√2)/2 | max 0.002σ |
| R3h 2504 horizontal cuts (41) | `ScalingProblem(gammas[bi][2:end], Ics[bi][2:end], ds[2:end]; sf=:xy, dx=[-1.25, 1.25], p_space=[0.3:0.01:1.0, 1.0:0.1:2.5, -0.4:0.1:0.4], quality=MultipleSplines(N_steps=100))` | max 4e-5σ |
| R4 2608 Fig. 4 insets | `Potts_Ic_horizontal/3`, `ScalingProblem([γs[k][bi] for k], Ics[:, bi], ds; sf=:xy, dx = bi==2 ? [-2,2] : [-1.5,1.5], p_space=[0.4:0.01:1, 0.1:0.1:3, -0.5:0.1:0.5], quality=MultipleSplines())` | exactly (N⁽²⁾, bi = 2); N⁽¹⁾ matches the paper text |
| R5h 2608 horizontal cuts (40) | `Potts_Ic_horizontal/1`, same call with `dx=[-1.5,1.5]`, `p_space=[0.5:0.01:1, 0.1:0.1:3, -0.5:0.1:0.5]` | 39/40 exactly; **cut 40 was bug-affected** (stored at quality 2.032, fixed code finds 1.866, γ_c moves 1.19σ) |
| R5v 2608 Fig. 10 (30) | see the spec ("Provenance of R5"); custom `AffineSimplexer(a=1e-6, b=1e-6)` | within 0.12σ |
| E2 (γ index 16, γ = 0.3) | fit 1.000575(161), 0.9981(48), 0.00104(62) | |

- **Software versions:** the deposits were produced with ScalingCollapse v0.2.2 + Optim 1.11 (2504,
  2505) and v0.2.5 + Optim 1.13.3 (2608).
- **Runtimes locally:** R3v 9 s, R3h 15 s, R5h 41 s, R5v about 6 s.
- **Preset conventions:** `:xy` is `y·L^(p3/p2)`; `:xny` is `y·L^(-p3/p2)`.
- **What `p3` means:** for an order parameter `m ~ L^(-β/ν)`, use `:xy` with `p3 = β`. For a
  susceptibility `χ ~ L^(γ/ν)`, use `:xny` with `p3 = γ`.
- **Local copies of the deposit files:**
  `/private/tmp/claude-501/-Users-malte-github-ScalingCollapse-jl/4aae5966-721c-4e31-ac8f-b6598c9c0ad2/scratchpad/{lt/LearningTransitions,potts/LearningPotts}/data/`.
  If they are gone, re-fetch the zip members by HTTP range read:
  - `LearningTransitions.zip` on Zenodo record 15227834;
  - `LearningPotts.zip` on Zenodo record 21918134;
  - `.../scratchpad/scripts/rzipx.py <file url> <outdir> <member…>`, where the url comes from
    `https://zenodo.org/api/records/<id>`.

## Review Focus

1. **A case whose fit throws** (e.g. a quality function that cannot handle a no-error-bar
   variant). It must fail that one case, with its name, and not abort the whole collapse suite.
   The runner wraps each case in its own `@testset` and catches errors into a failing `@test`.
2. **A regression case missing from `references.toml`.** It must fail with a message naming the
   case and pointing to `update_references.jl`, not pass silently. Task 1 tests this.
3. **`update_references.jl` run on a file with `[published]` entries.** It must keep them
   byte-for-byte as values, and only `[last]` changes. Task 1 tests this.
4. **Parameter-count mismatch** between a case's reference and the fit (e.g. after fixing a
   parameter). The checker must report a failure, not throw a `BoundsError` or zip silently.
   Task 1 tests this.
5. **Runtime creep on CI** (x86 Windows). The suite prints its total time and per-group times, and
   Task 6 checks the local budget.

---

### Task 1: Framework, references, and migrating the existing collapses (E1, R1, E6, R6)

**Files:**
- Create: `test/collapses/framework.jl`. Holds `CollapseCase`, `run_case`, `check_exact`,
  `check_regression`, `load_references`, `print_rows`.
- Create: `test/collapses/runcollapses.jl`. Includes the framework and all `cases_*.jl`, runs
  every case, and prints the timing.
- Create: `test/collapses/cases_qubitloss.jl` (E1, R1), `test/collapses/cases_ising.jl` (E6),
  `test/collapses/cases_cdw.jl` (R6).
- Create: `test/collapses/references.toml` and `test/collapses/update_references.jl`.
- Create: `test/collapses/framework_tests.jl`, which holds the checker unit tests.
- Create: `test/testdata/cdw_holstein.jl`. Move the `X_L*`/`Y_L*`/`E_L*` constants here
  verbatim from `test/cdw.jl`.
- Modify:
  - `test/cdw.jl`: replace the constants with `include("testdata/cdw_holstein.jl")`;
  - `test/qubitloss.jl`: delete it, since its two checks become E1 and R1 θ = 0;
  - `test/runtests.jl`: drop `include("qubitloss.jl")`, add
    `include("collapses/framework_tests.jl")` and `include("collapses/runcollapses.jl")` last;
  - `test/Project.toml`: add TOML.

**Interfaces:**
- **Produces:**
  - `CollapseCase(; name::String, kind::Symbol, load, kwargs::NamedTuple, exact=Dict{Int,Float64}(), δabs=Dict{Int,Float64}(), σmax=Dict{Int,Float64}(), source::String)`.
    `kind ∈ (:exact, :regression)`. `load()` returns `(xs, ys, es, Ls)` with `es === nothing`
    when there are no error bars. Dict keys index the free parameters.
  - `run_case(c)::ScalingProblem`, which calls `ScalingProblem(xs, ys, [es,] Ls; c.kwargs...)`.
  - `check_exact(c, sp)::Vector{CheckRow}` and
    `check_regression(name, sp, ref::Dict{String,Any})::Vector{CheckRow}`, where
    `ref = Dict("optimal_ps"=>…, "optimal_ps_error"=>…)`.
  - `CheckRow = @NamedTuple{param::Int, expected::Float64, got::Float64, dev_sigma::Float64, err_ratio::Float64, pass::Bool, note::String}`.
  - `load_references(path) -> Dict{String,Any}`, with top-level keys `"published"` and `"last"`,
    each mapping a case name to its ref dict.
  - `CASES::Vector{CollapseCase}`, a global that each `cases_*.jl` `append!`s to.
  - `update_references.jl`, run from `test/` as `julia --project=.. collapses/update_references.jl`.
    It runs every `:regression` case and rewrites only `[last]`, with keys sorted, and stores
    `optimal_ps`, `optimal_ps_error` and `minimum`.

- [ ] **Step 1: Write failing checker tests in `test/collapses/framework_tests.jl`.**
  - Use fake fit results: `fake(ps, errs) = (optimal_ps=ps, optimal_ps_error=errs, minimum=1.0)`.
  - `check_exact` with `exact=Dict(1=>0.5, 2=>4/3)`:
    - `fake([0.50001, 1.34], [1e-5, 0.01])` passes;
    - `fake([0.5, 1.3775], [1.9e-4, 0.088])` with `σmax=Dict(1=>1.3e-4, 2=>1.8e-2)` fails both
      rows on σmax;
    - `fake([0.51, 4/3], [1e-3, 0.01])` fails row 1, and passes with `δabs=Dict(1=>0.01)`.
  - `check_regression`:
    - `ref=Dict("optimal_ps"=>[1.0, 2.0], "optimal_ps_error"=>[0.1, 0.1])` with
      `fake([1.02, 2.0], [0.1, 0.1])` passes (0.2σ);
    - `fake([1.03, 2.0], …)` fails (0.3σ);
    - `fake([1.0, 2.0], [0.04, 0.1])` fails (ratio 0.4);
    - `fake([1.0, 2.0], [0.21, 0.1])` fails (ratio 2.1);
    - `fake([1.0], [0.1])` returns one failing row with a note containing `"parameter count"`.
  - A missing reference: `check_regression("nope", sp, nothing)` returns one failing row with a
    note containing `"update_references.jl"`.
  - `update_references`: factor its writer as `merge_last(old::Dict, new_last::Dict) -> Dict`.
    Test that `merge_last(Dict("published"=>Dict("a"=>Dict("optimal_ps"=>[1.0])), "last"=>Dict()), Dict("a"=>Dict("optimal_ps"=>[2.0])))["published"]`
    equals the input `published`.

- [ ] **Step 2: Run the tests and confirm they fail** with `UndefVarError: check_exact`.
  - Run: `cd test && julia --project=.. -e 'using Pkg; Pkg.test()'`.

- [ ] **Step 3: Implement `framework.jl`** with the interfaces above.
  - **Runner:** for each case, open a `@testset "<name>"`. Inside, `try sp = run_case(c)` and test
    every row's `pass`; on error, `@test false` with the case name in the message. On any failing
    row, call `print_rows` (a table of param / expected / got / devσ / ratio / note).
  - **Which references a case is checked against:**
    - `:exact` cases are checked against `c.exact` only.
    - `:regression` cases are checked against `[last]` always, and against `[published]` when the
      case has an entry there.
  - **Timing:** print the total time of `runcollapses.jl` with `@info`.

- [ ] **Step 4: Write the first cases.**
  - **E1** (`qubitloss_theta_pi2_exact`): the θ = π/2 case from `test/qubitloss.jl`.
    `exact = Dict(1=>0.5, 2=>4/3)`, `σmax = Dict(1=>1.3e-4, 2=>1.8e-2)`.
  - **R1** (`qubitloss_theta{i}`, i = 1…15): the paper call. Use `Ls = [32, 64, 128, 256, 512]` for
    i = 1, otherwise `[32, 64, 128, 256]`.
  - **E6** (`ising_mc_binder_exact`): `IsingMonteCarlo.h5` binder, with `:x`, `dx=[-1, 1]`,
    `p_space=[0.1:0.1:3, 0.1:0.1:3]` and `exact = Dict(1=>2/log(1+√2), 2=>1.0)`. Measure first;
    set `δabs` and `σmax` per the Global Constraints, with a comment.
  - **R6** (`cdw_holstein`): the first `ScalingProblem` call of `test/cdw.jl`, verbatim.

- [ ] **Step 5: Seed `references.toml`.**
  - `[published.qubitloss_theta1]`: `optimal_ps = [0.33681399873460816, 1.5316126320454644]`,
    `optimal_ps_error = [2.0879977065466418e-5, 0.004006908264828457]`.
  - `[published.qubitloss_theta{i}]` for i = 2…15: ν only is not representable as a vector, so
    store r_c as well, using `NaN` for its value and error. The checker skips NaN reference entries
    and notes `"no published value"`. Add this NaN-skip to Step 1's tests first (red, then green).
    - Published ν (value, error) for i = 2…15:
      1.5349(44), 1.5320(48), 1.5328(52), 1.5325(60), 1.5298(75), 1.5287(71), 1.5295(66),
      1.5267(57), 1.5231(57), 1.5215(51), 1.5122(49), 1.4886(70), 1.4206(83), 1.3370(59).
    - These are rounded to the notebook's printing: the error has 2 significant digits and the
      value is printed to its error's last digit. Use those decimals as given.
  - Then run `update_references.jl` to create `[last]`. Review the diff: `[published]` unchanged,
    and every `[last]` R1 ν within 0.25σ of its published value.

- [ ] **Step 6: Run the full suite.**
  - Run: `cd test && julia --project=.. -e 'using Pkg; Pkg.test()'`.
  - Expected: all pass, including the moved CDW and Ising tests, unchanged.

- [ ] **Step 7: Commit:** `Add collapse benchmark framework with qubit-loss, Ising and CDW cases`.

### Task 2: Exact noise-free generators (E5 Ising strip, E7 Curie–Weiss)

**Files:**
- Create: `test/generators/ising_strip.jl` and `test/generators/curie_weiss.jl`. Each writes
  `test/testdata/exact/<name>.h5` (layout: `Ls`, `L{L}/{x, y}`) and needs only HDF5 (and
  LinearAlgebra for the cross-check).
- Create: `test/testdata/exact/ising_strip.h5` and `test/testdata/exact/curie_weiss.h5`.
- Create: `test/collapses/cases_exact_models.jl`.

**Interfaces:**
- Consumes: `CollapseCase` and `CASES` from Task 1.

- [ ] **Step 1: Write the E5 generator.**
  - Compute the correlation length of a periodic Ising strip of width L from Kaufman's
    transfer-matrix spectrum (J = 1). With `K = 1/T` and `K*` defined by `tanh K* = e^{-2K}`, let
    `cosh γ_q = cosh 2K* cosh 2K − sinh 2K* sinh 2K cos(πq/L)`, and take
    `ξ_L⁻¹ = ½ Σ_{r=0}^{L−1} (γ_{2r} − γ_{2r+1})`, where γ₀ carries the sign of `K − K*`.
  - Before writing any file, the generator must **assert** two things:
    1. `ξ_L` matches `1/log(λ₀/λ₁)` from brute-force diagonalisation of the 2^L × 2^L
       transfer matrix, to `1e-10` relative, for L = 4, 6, 8 at T ∈ {2.0, 2.269, 2.6}. If the
       sign convention above is wrong, fix it here.
    2. `ξ_L/L` at `T_c = 2/ln(1+√2)` tends to `4/π` (within `1e-3` at L = 256).
  - **Output:** `y = ξ_L/L` for L = 8, 16, 32, 64, on 61 points of T ∈ [2.1, 2.45].

- [ ] **Step 2: Write the E7 generator.**
  - Use `Z_N(T) = Σ_k binomial(N,k)·exp(M²/(2NT))`, with `M = 2k − N`, in log-space via
    `logabsbinomial` or `lgamma` sums.
  - Outputs, for N = 64, 128, 256, 512, 1024 on 61 points of T ∈ [0.8, 1.2]:
    - the Binder ratio `U = 1 − ⟨M⁴⟩/(3⟨M²⟩²)`;
    - `m_rms = sqrt(⟨M²⟩)/N`.

- [ ] **Step 3: Write the cases in `cases_exact_models.jl`.** Use no `es` and the default quality.
  - **E5** `ising_strip_xi_exact`: `:x`, `exact = Dict(1=>2/log(1+√2), 2=>1.0)`.
  - **E7a** `curie_weiss_binder_exact`: `:x` with the system size passed as `N`, so
    `exact = Dict(1=>1.0, 2=>2.0)` (x scales as `N^{1/2}`).
  - **E7b** `curie_weiss_mrms_exact`: `:xy`, `exact = Dict(1=>1.0, 2=>2.0, 3=>0.5)` (from
    `p3/p2 = 1/4`).
  - Choose `dx`/`p_space` so that the truth is interior. Run the cases, then set `σmax` to 3× the
    measured errors. `δabs = 0` for all of them. If a case fails at 3σ, that is a finding about
    the package or the generator, not about the tolerance: investigate (check the generator's
    asserts first) and ledger it.

- [ ] **Step 4: Run the suite and confirm it passes;** record the fitted values in the ledger.

- [ ] **Step 5: Commit:** `Add exact Ising strip and Curie-Weiss collapse cases`.

### Task 3: Synthetic family (E8)

**Files:**
- Create: `test/generators/synthetic.jl`, which writes `test/testdata/synthetic/synthetic.h5`
  with one group per variant.
- Create: `test/collapses/cases_synthetic.jl`.

- [ ] **Step 1: Generate the data.**
  - **Model:** `y = L^{s·p3/p2} · F((x − x_c)/x_c · L^{1/p2}) · (1 + c·L^{−ω}) + noise`, with
    `x_c = 0.7`, `p2 = ν = 1.2` and `p3 = 0.3`. `s = −1` for `:xy` variants, `+1` for `:xny`,
    and `p3` is absent for `:x`.
  - **F:** monotone `F(z) = 1 + tanh(z)`, or peaked `F(z) = exp(−z²) + 0.1`.
  - **Defaults:** `L = [16, 32, 64, 128]`, 41 equidistant x ∈ [0.55, 0.85], `c = 0`, no noise.
  - **Variants:**

    | Variant | Changes from the defaults |
    |---|---|
    | `clean_x` | `:x`, monotone F |
    | `clean_xy` | `:xy`, monotone F |
    | `clean_xny` | `:xny`, peaked F |
    | `noise_err` | `:xy`; Gaussian noise with σ = 1% of \|y\|; `es` stored |
    | `noise_noerr` | the same data as `noise_err`, without `es` |
    | `three_sizes` | `:xy`; L = [16, 32, 64] |
    | `uneven_grid` | `:xy`; per L, 41 sorted uniform-random x ∈ [0.55, 0.85] |
    | `corrections` | `:xy`; c = 0.5, ω = 1, plus 1% noise with `es` |
    | `fixed_nu` | `:xy`, with `ScalingFunction(:xy; p2=1.2)` in the case (the default names are `p1, p2, p3`) |

- [ ] **Step 2: Write the cases.**
  - Run each variant with each of `Linear()`, `SingleSpline()`, `MultipleSplines()` and
    `Houdayer()`; the case name is `synthetic_<variant>_<quality>`.
  - Use `dx = [-1.5, 1.5]` and `p_space = [0.5:0.01:0.9, 0.6:0.1:2.0, 0.0:0.05:0.6]`, dropping
    dimensions for `:x` and `fixed_nu`.
  - All variants are `:exact` with the true values, except `corrections`. That one is
    `:regression`, plus an extra exact sanity bound `|p − p*| ≤ 0.05·|p*|`, implemented as
    `δabs = 0.05·|p*|` in a second, exact case named `..._sanity` (Claude's call: the fit is
    biased by design).
  - Set `σmax` to 3× the measured error.

- [ ] **Step 3: Run and triage.**
  - Run: `cd test && julia --project=.. -e 'using Pkg; Pkg.test()'`.
  - **A combination that throws or fails:** check whether it is a limitation of that quality
    function by design (read its docstring and source). If so, exclude it with a comment and
    ledger it.
  - **Otherwise:** it is a package finding. Ledger it, and leave the case marked
    `@test_broken`-style (use `kind = :exact` with a `broken = true` field added to
    `CollapseCase`, defaulting to false) so it is visible but not red. Report every such finding
    to the user.
  - Add `broken` handling to the Task 1 framework tests (red, then green) before using it.

- [ ] **Step 4: Commit:** `Add synthetic collapse family across presets and quality functions`.

### Task 4: Percolation (E4 self-generated, E3 autoScale)

**Files:**
- Create: `test/generators/percolation.jl`, which writes `test/testdata/percolation/bond_square.h5`.
- Create: `test/testdata/autoscale/orderParam_L{8,16,32,64,128,256,512}.dat` and
  `test/testdata/autoscale/LICENSE`.
- Create: `test/collapses/cases_percolation.jl`.

- [ ] **Step 1: Write the E4 generator.**
  - **Method:** Newman–Ziff bond percolation on an L × L square lattice, open boundaries, using
    union-find with path compression and union by size.
  - **Per sample:** add bonds in random order and record after each bond
    - top-bottom spanning (0/1), and
    - the largest cluster fraction.
  - **Canonical average:** convolve each sample with binomial weights at
    p ∈ 0.45:0.001:0.55, computed in log-space. Report the mean and standard error over samples.
  - **Sizes and samples:** L = 16, 32, 64, 128, 256, with 10⁴ samples for L ≤ 128 and 2·10³ for
    L = 256. The seed is fixed.
  - **Assert** that the spanning probability at p = 0.5 lies in [0.4, 0.6] for every L, as a
    sanity check.

- [ ] **Step 2: Fetch E3.**
  - `curl` the seven `ORDER_PARAMETER/orderParam_L*.dat` files and `LICENSE` from
    `https://raw.githubusercontent.com/omelchert/autoScale/<default branch>/`. Look up the default
    branch with `gh api repos/omelchert/autoScale --jq .default_branch`.
  - Columns are p, y, error, and `#` lines are comments.
  - Record the source URL and commit SHA in `test/testdata/README.md`.

- [ ] **Step 3: Write the cases.** Set `δabs` and `σmax` as in the Global Constraints, with
  comments.
  - **E4a** `bond_spanning_exact`: `:x`, `exact = Dict(1=>0.5, 2=>4/3)`.
  - **E4b** `bond_strength_exact`: `:xy`, `exact = Dict(1=>0.5, 2=>4/3, 3=>5/36)`.
  - **E3** `autoscale_site_order_exact`: `:xy`,
    `exact = Dict(1=>0.59274605079210, 2=>4/3, 3=>5/36)`, with `dx` from autoScale's own fit
    window, mapped to the relative x convention.

- [ ] **Step 4: Run the suite and confirm it passes. Commit:**
  `Add bond and site percolation collapse cases`.

### Task 5: Paper data (E2, R2, R3, R4, R5)

**Files:**
- Create: `test/generators/convert_papers.jl`. It runs in a scratch environment with JLD2,
  Measurements and HDF5, and writes value/error arrays per curve.
- Create:
  - `test/testdata/learning_ising/{betac,vertical,horizontal}.h5`, for 2504.12385, with only the
    sizes the calls use: vertical `ds[2:5]`, horizontal `ds[2:end]`;
  - `test/testdata/learning_potts/{insets,horizontal,vertical}.h5`, for 2608.20268, from
    `Potts_Ic_horizontal/3`, `Potts_Ic_horizontal/1` and `Potts_Ic_vertical/4`.
- Create: `test/collapses/cases_learning.jl`.
- Modify: `test/collapses/references.toml` (`[published]` entries) and `test/testdata/README.md`.

- [ ] **Step 1: Write and run the converter.** Then check the conversion: refitting a few curves
  from the HDF5 file must give the same numbers as refitting from the JLD2 file.

- [ ] **Step 2: Write the cases**, using the calls in the "Facts" table verbatim:
  - **E2** `learning_ising_vertical_g16_exact`: `exact = Dict(1=>1.0, 2=>1.0)`.
    `δabs = Dict(1=>1e-3)`, with the comment "published fit 1.00058(16) is 3.6σ from 1;
    small-size bias".
  - **R2** `learning_ising_betac`.
  - **R3** `learning_ising_vertical_g{i}` (31 cuts) and `learning_ising_horizontal_b{i}`
    (41 cuts).
  - **R4** `learning_potts_inset_N1` and `learning_potts_inset_N2`.
  - **R5** `learning_potts_horizontal_b{i}` (40 cuts) and `learning_potts_vertical_g{i}`
    (30 cuts), the latter with the custom simplex.

- [ ] **Step 3: Fill `[published]`.**
  - **Where the values come from:** the deposits' stored fit files (`FSS_vertical.jld2`,
    `FSS_horizontal.jld2`, `Potts_Ic_horizontal/1/FSS.jld2`, `Potts_Ic_vertical/4/FSS.jld2`).
  - **Paper numbers:**
    - R2: `[0.5980465451677133, 2.6581251587998618, 0.041657264644830376] ±
      [0.0020456949924619735, 0.05404956842410173, 0.003291440393339273]`;
    - R4 N⁽²⁾: `[0.5806695158504818, 2.475730758904198, 0.06933984903666862] ±
      [0.0012880674176765387, 0.033199220533840834, 0.0032187643446133464]`;
    - R4 N⁽¹⁾: 0.76383(85), 1.453(33), −0.0004(4), from the paper text.
  - **Exclude `learning_potts_horizontal_b40` from `[published]`.** Add a comment on that case:
    "stored fit (v0.2.5 + Optim 1.13.3) was hit by the Optim ≥ 1.12 simplex bug: quality 2.032
    vs 1.866 at the true minimum; γ_c differs by 1.19σ. Checked against [last] only."
  - Then run `update_references.jl`. Review: `[last]` agrees with `[published]` within 0.25σ for
    every case that has both.

- [ ] **Step 4: Run the full suite, measure the time, and commit:**
  `Add paper collapse cases from arXiv:2504.12385 and 2608.20268`.

### Task 6: Tolerance audit, runtime, docs, PR

- [ ] **Step 1: Tolerance audit.**
  - For every exact case, list `|p − p*|/σ`, `σ/σmax` and `δabs` (with its justification) in the
    PR description.
  - Because `σmax ≈ 3σ_fit`, every exact case already catches 15× error inflation. Flag every
    exact case where `3σ_fit + δabs` on the ν-type parameter is ≥ 0.04: a shift of the size the
    Optim bug produced would pass there. Ledger those cases and explain why that is acceptable.

- [ ] **Step 2: Mutation check.** Revert #43's halving locally (do not commit), with Optim 2.3.2
  resolved. Record which cases fail; at least E1 and R1 θ = 0 must fail. Restore the code.

- [ ] **Step 3: Runtime.** Run the full `Pkg.test()` locally. If it takes over 4 min, thin R3/R5
  per the Global Constraints and ledger the cuts kept.

- [ ] **Step 4: Docs.**
  - Finish `test/testdata/README.md` with all sources, licenses and transformations.
  - Add a short "Collapse benchmark suite" section to `CLAUDE.md` covering:
    - where the cases live;
    - how to add one;
    - how to run `update_references.jl`;
    - that `[published]` is never auto-written.

- [ ] **Step 5: Commit, push and open the PR** (target per Global Constraints).
  - Its body includes:
    - the case inventory;
    - the audit table;
    - the mutation result;
    - the runtime;
    - the R5 cut-40 finding;
    - any `broken` synthetic findings.
  - The body ends with
    `🤖 Generated with [Claude Code](https://claude.com/claude-code)`.
  - Bind the PR with the ccd_pr tools.
