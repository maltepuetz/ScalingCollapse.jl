# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

ScalingCollapse.jl is a Julia package (Julia ≥ 1.9) for automatic finite-size scaling: given observables measured at several system sizes `L`, it finds the scaling parameters (e.g. `T_c`, `nu`, `beta`) that best collapse all curves onto one master curve.

## Commands

```bash
# Run the full test suite (includes Aqua.jl code-quality checks)
julia --project=. -e 'using Pkg; Pkg.test()'

# Build the docs locally (DocumenterVitepress; needs node/npm; output in docs/build/)
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
julia --project=docs docs/make.jl

# Live preview of the built site (from docs/)
npm run docs:dev
```

On Julia ≥ 1.11 `Pkg.develop` writes a `[sources]` entry with an **absolute** path into the tracked `docs/Project.toml`. Check `git diff docs/Project.toml` afterwards and revert it. A docs build also kills a running `vitepress dev` server, because it regenerates `docs/build/.documenter`; restart it with `npm run docs:dev`.

`docs/src/.vitepress/config.mts` is **vendored** from DocumenterVitepress's template (same setup as QuditClifford.jl), so it is frozen at the version it was copied from. DocumenterVitepress substitutes its own copy only when this file is absent, and fills the `REPLACE_ME_DOCUMENTER_VITEPRESS*` placeholders in ours. On a DocumenterVitepress bump, diff it against `template/src/.vitepress/config.mts` of the new version and re-apply the two local changes named in its header: the navbar `logo` (the `logo-light.svg`/`logo-dark.svg` pair) and the favicon links plus the head script that picks between `favicon-light.ico` and `favicon-dark.ico`. Afterwards confirm that no `REPLACE_ME` is left in `docs/build/.documenter/.vitepress/config.mts`. The logo and favicons live in `docs/src/assets/`.

Running a single test file: test files are `include`d from `test/runtests.jl`, which does `using ScalingCollapse, Test, Aqua, HDF5`. Test-only deps live in `test/Project.toml`. Data-driven tests open `testdata/IsingMonteCarlo.h5` via a path relative to `test/`, so the working directory must be `test/`. The easiest approach is a REPL using TestEnv.jl (`using TestEnv; TestEnv.activate()`), then `cd("test")`, `using ScalingCollapse, Test, HDF5`, `include("ising.jl")`. Alternatively, temporarily comment out includes in `runtests.jl`.

CI (`.github/workflows/CI.yml`) runs on Julia `lts`, `1` and `nightly`, across Linux/macOS/Windows, including 32-bit Windows (x86). Watch for `Int`/`Int64` assumptions, since a 32-bit vs 64-bit bug was fixed before.

## Architecture

Everything happens in the `ScalingProblem` constructor: `ScalingProblem(args...; kwargs...)` parses kwargs, builds `Vector{Data}`, and calls `solve!` immediately. There is no separate solve step for users. Results are stored on the mutable struct (`optimal_ps`, `optimal_ps_error`, `minimum`).

Pipeline (`src/scalingproblem.jl` → `src/optimization.jl` → `src/error_analysis.jl`):
1. **Data ingestion** (`src/data.jl`): `unzip_data` accepts many input shapes for `xs`/`ys`/`es`: a shared vector, a vector of vectors, or a matrix whose orientation is inferred from `length(Ls)`. Each size becomes a `Data(L, xs, ys, es)`, sorted by x. `Measurement` y-values are split into value and uncertainty. If any `es` is nonzero, `errors_defined=true`, which switches quality functions to their weighted variants.
2. **Parameter scan** (`parameter_scan!`): evaluates `quality_scan` (default `MultipleSplines(scan_mode=true)`) on the full Cartesian grid `p_space`, finds grid-local minima, runs Nelder–Mead from each one, and keeps the best as `starting_ps`. Passing `starting_ps` skips this step. In that case `p_space` is set to the sentinel `[-1_000_000:1_000_000 ...]`, which `Base.show` checks for.
3. **Final optimization** (`optimize_starting_ps!`): Nelder–Mead on `quality` (default `Linear()`) from `starting_ps`. Both stages use the user-configurable `optimizer::Optim.NelderMead` and `optimizer_options::Optim.Options`. Every fit wraps `optimizer.initial_simplex` in `BoundedSimplexer`. The exported `StepSimplexer(steps)` takes one step per free parameter; profile inner fits pass it through `_restrict_simplexer` so it keeps only the nuisance parameters' steps.
4. **Error analysis** (`error_analysis!`, on by default; `error` is `true`/`:profile`, `:conditional` or `false`): the error of a parameter is the larger of the left and right displacements at which `quality/minimum` reaches `error_threshold` (default 1.3). `:profile` (the default) re-fits the other free parameters inside `p_space` at every displacement and takes the first crossing moving outward from the fit (a local width); the profiled parameter is not bounded by `p_space`. A lower quality found during the profile re-optimizes the fit and restarts the analysis (up to 3 times), so `optimal_ps` and `minimum` can change; exhausted restarts give all-NaN errors. `:conditional` is the previous method: bisection with the other parameters held at their optimum. With one free parameter both agree. `parameter_errors(sp; method, ...)` (keywords default to the settings of `sp`) runs the analysis on a copy of a solved problem and throws instead of re-optimizing; `parameter_errors!` stores the errors and settings in `sp` and may re-optimize.

**Quality functions** (`src/quality.jl`, `src/quality_functions/*.jl`): subtypes of the abstract `QualityFunction` are callable structs with the signature `(qf)(sp, ps; check_bounds=false) -> Float64`. They apply `sp.sf.f` to every `Data`, restrict to the overlap interval intersected with `sp.dx`, and return `Inf` when out of bounds (`check_bounds=true` compares against `sp.p_space`) or when the interval is empty. Available functions: `Linear`, `SingleSpline`, `MultipleSplines` (Dierckx splines), and `Houdayer`. To add a new one, create a file in `quality_functions/`, `include` it in `quality.jl`, and export it from `src/ScalingCollapse.jl`.

**Scaling functions** (`src/scalingfunction.jl`): `ScalingFunction` wraps a user or preset function `f(d::Data, p1, p2, ...) -> Data`. Presets are `:x` (2 params), `:xy` (3, the default), and `:xny` (3, negative y exponent). They rescale x as `(x - p1) * L^(1/p2)`; `relative=true` selects the form `(x - p1)/p1 * L^(1/p2)`, which was the default up to 0.2 (a relative window `dx` equals `dx .* p1` in the absolute form). They dispatch through `handle_preset(::Val{preset}; kwargs...)`, which is the intended extension point for new presets. Fixing a parameter is done by passing its name as a kwarg (e.g. `nu=1.0`). `specialize_wrapper` then builds `sf.f`, which takes only the free parameters. So `sf.f`, `n_parameters(sf)`, `p_space`, `starting_ps`, and `optimal_ps` all refer to free parameters only. `p_names` lists all of them; use `scaled_p_names` for the free ones.

**Other public API**: `scaled_data(sp[, ps]; splines=false)` returns rescaled data for plotting. `residuals(sp; ...)` (`src/residuals.jl`) evaluates the quality landscape over a parameter grid.

`ext/InteractiveScaling.jl` is an entirely commented-out GLMakie plotting extension. It is not registered in `Project.toml` (no `[weakdeps]`/`[extensions]`) and is currently inactive.

## Collapse benchmark suite

`test/collapses/` holds about 200 end-to-end collapse cases, run by `Pkg.test()` through `collapses/runcollapses.jl`. `framework.jl` defines `CollapseCase`, the runner and the two checks. The cases live in `cases_*.jl`, each of which `push!`es onto `CASES`. Their data live in `test/testdata/`, with sources, licenses and transformations in `test/testdata/README.md`; generators and converters live in `test/generators/` and are not run by the tests.

- **Exact cases** (`kind=:exact`) compare with known values: `|p − p*| ≤ 3σ_fit + δabs` and `σ_fit ≤ σmax`, with `σmax` ≈ 3× the fitted error of the current code. `δabs` is only for data with a measured finite-size bias, with a comment giving the reason and the measured deviation. Never loosen a criterion to get a pass.
- **Regression cases** (`kind=:regression`) compare with `test/collapses/references.toml`: `|p − p_ref| ≤ 0.25σ_ref` and `0.5 ≤ σ_fit/σ_ref ≤ 2`, against `[last]` and, where it exists, `[published]`.
- `broken=true` documents a known failure: the case is reported as one `@test_broken` and flags an error as soon as it passes.

To add a case, `push!` a `CollapseCase` in the matching `cases_*.jl` (a new `cases_*.jl` file is picked up automatically), commit its data under `test/testdata/` and document the data in its README. For a regression case, then regenerate `[last]`.

`update_references.jl` rewrites the `[last]` table from the current code; review the diff before committing it. It needs the test dependencies (HDF5, TOML), which the root project lacks, so run it from the repo root in the test environment, either with TestEnv.jl (installed in your default environment):

```bash
julia --project=. -e 'using TestEnv; TestEnv.activate(); include("test/collapses/update_references.jl")'
```

or in a temporary environment:

```bash
julia -e 'using Pkg; Pkg.activate(; temp=true); Pkg.develop(path="."); Pkg.add(["HDF5", "TOML"]); include("test/collapses/update_references.jl")'
```

`[published]` (values from papers and their data deposits) is never written by `update_references.jl`; edit it by hand only, with the source in the case's comment.
