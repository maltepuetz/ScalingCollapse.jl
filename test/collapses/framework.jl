# Framework of the collapse benchmark suite: cases, runner and checkers.
using ScalingCollapse
using Test
using TOML

"""
    CollapseCase(; name, kind, load, kwargs, exact, δabs, σmax, broken, source)

One benchmark collapse. `kind` is `:exact` (compared with known exact values `exact`, keyed by
free-parameter index) or `:regression` (compared with `references.toml`). `load()` returns
`(xs, ys, es, Ls)` with `es === nothing` when there are no error bars; `kwargs` go to
`ScalingProblem`. A `broken` case documents a known failure: its rows are reported with
`@test_broken`, so it stays green while failing (a throw counts as failing) and is flagged as an
error as soon as it passes.
"""
Base.@kwdef struct CollapseCase
    name::String
    kind::Symbol
    load::Any
    kwargs::NamedTuple
    exact::Dict{Int,Float64} = Dict{Int,Float64}()
    δabs::Dict{Int,Float64} = Dict{Int,Float64}()
    σmax::Dict{Int,Float64} = Dict{Int,Float64}()
    broken::Bool = false
    source::String
end

const CheckRow = @NamedTuple{param::Int, expected::Float64, got::Float64, dev_sigma::Float64, err_ratio::Float64, pass::Bool, note::String}

# every cases_*.jl appends its cases here
const CASES = CollapseCase[]

function run_case(c::CollapseCase)
    xs, ys, es, Ls = c.load()
    es === nothing && return ScalingProblem(xs, ys, Ls; c.kwargs...)
    return ScalingProblem(xs, ys, es, Ls; c.kwargs...)
end

"""
Exact criterion: `|p - p*| <= 3σ_fit + δabs` and, where given, `σ_fit <= σmax`.
`err_ratio` is `σ_fit / σmax` (NaN without σmax).
"""
function check_exact(c::CollapseCase, sp)
    rows = CheckRow[]
    for i in sort!(collect(keys(c.exact)))
        expected = c.exact[i]
        got, σ = sp.optimal_ps[i], sp.optimal_ps_error[i]
        δ = get(c.δabs, i, 0.0)
        dev = abs(got - expected)
        notes = String[]
        dev <= 3σ + δ || push!(notes, "|p-p*| > 3σ+δabs")
        σmax = get(c.σmax, i, NaN)
        isnan(σmax) || σ <= σmax || push!(notes, "σ_fit > σmax")
        push!(rows, (param=i, expected=expected, got=got, dev_sigma=dev / σ,
            err_ratio=σ / σmax, pass=isempty(notes), note=join(notes, "; ")))
    end
    return rows
end

"""
Regression criterion against `ref = Dict("optimal_ps"=>…, "optimal_ps_error"=>…)`:
`|p - p_ref| <= 0.25σ_ref` and, if `check_errors`, `0.5 <= σ_fit/σ_ref <= 2`. NaN reference
entries mean "no published value" and pass. `[published]` is checked with `check_errors=false`,
since its stored errors come from the conditional method.
"""
function check_regression(name, sp, ref; check_errors::Bool=true)
    nan_row(note) = [(param=0, expected=NaN, got=NaN, dev_sigma=NaN, err_ratio=NaN, pass=false, note=note)]
    ref === nothing && return nan_row("no reference for $name; run collapses/update_references.jl")
    ps, σs = Float64.(ref["optimal_ps"]), Float64.(ref["optimal_ps_error"])
    length(ps) == length(sp.optimal_ps) ||
        return nan_row("parameter count $(length(sp.optimal_ps)) != reference $(length(ps))")
    rows = CheckRow[]
    for i in eachindex(ps)
        got, σ = sp.optimal_ps[i], sp.optimal_ps_error[i]
        if isnan(ps[i]) || isnan(σs[i])
            push!(rows, (param=i, expected=NaN, got=got, dev_sigma=NaN, err_ratio=NaN,
                pass=true, note="no published value"))
            continue
        end
        dev, ratio = abs(got - ps[i]) / σs[i], σ / σs[i]
        notes = String[]
        dev <= 0.25 || push!(notes, "|p-p_ref| > 0.25σ_ref")
        !check_errors || 0.5 <= ratio <= 2 || push!(notes, "σ_fit/σ_ref outside [0.5, 2]")
        push!(rows, (param=i, expected=ps[i], got=got, dev_sigma=dev, err_ratio=ratio,
            pass=isempty(notes), note=join(notes, "; ")))
    end
    return rows
end

function load_references(path)
    ref = TOML.parsefile(path)
    get!(ref, "published", Dict{String,Any}())
    get!(ref, "last", Dict{String,Any}())
    return ref
end

"Replace only the `last` table of `old`; `published` is never touched."
function merge_last(old::Dict, new_last::Dict)
    merged = copy(old)
    merged["last"] = new_last
    return merged
end

function print_rows(name, rows, label="")
    cell(x::Real) = lpad(string(round(x; sigdigits=6)), 12)
    cell(x::Integer) = lpad(string(x), 12)
    println("  ", name, isempty(label) ? "" : " [$label]")
    println("    ", join(lpad.(["param", "expected", "got", "devσ", "ratio"], 12)), "  note")
    for r in rows
        println("    ", cell(r.param), cell(r.expected), cell(r.got), cell(r.dev_sigma),
            cell(r.err_ratio), "  ", r.pass ? "ok " : "FAIL ", r.note)
    end
end

"""
Run `rows` through `@test`, printing the table when any row fails. For a `broken` case all rows
count as one `@test_broken` over their conjunction (so partly passing cases stay Broken), and no
table is printed.
"""
function report_rows(name, rows, label="", broken=false)
    if broken
        @test_broken all(r.pass for r in rows)
        return
    end
    for r in rows
        @test r.pass
    end
    all(r.pass for r in rows) || print_rows(name, rows, label)
end

function run_cases(cases, refs; runner=run_case)
    for c in cases
        @testset "$(c.name)" begin
            try
                sp = runner(c)
                groups = Tuple{Vector{CheckRow},String}[]
                if c.kind === :exact
                    push!(groups, (check_exact(c, sp), "exact"))
                else
                    push!(groups, (check_regression(c.name, sp, get(refs["last"], c.name, nothing)), "last"))
                    haskey(refs["published"], c.name) &&
                        push!(groups, (check_regression(c.name, sp, refs["published"][c.name]; check_errors=false), "published"))
                end
                if c.broken  # one result over the rows of all groups
                    report_rows(c.name, reduce(vcat, first.(groups)), "", true)
                else
                    for (rows, label) in groups
                        report_rows(c.name, rows, label)
                    end
                end
            catch err
                c.broken ? (@test_broken false) : (@test false)
                c.broken || @error "collapse case $(c.name) threw" exception = (err, catch_backtrace())
            end
        end
    end
end
