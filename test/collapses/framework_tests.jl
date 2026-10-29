# Unit tests of the collapse benchmark checkers, using fake fit results.
if !@isdefined(CollapseCase)
    include("framework.jl")
end

import Test: AbstractTestSet, record, finish
# collects results without printing or throwing
struct NoThrowTestSet <: AbstractTestSet
    description::String
    results::Vector{Any}
    NoThrowTestSet(desc; kw...) = new(desc, Any[])
end
record(ts::NoThrowTestSet, t) = (push!(ts.results, t); t)
record(ts::NoThrowTestSet, t::Test.AbstractTestSet) = (append!(ts.results, t.results); t)
finish(ts::NoThrowTestSet) = ts

@testset "Collapse benchmark framework" begin
    fake(ps, errs) = (optimal_ps=ps, optimal_ps_error=errs, minimum=1.0)
    dummy_load() = error("not used")
    exact_case(; kw...) = CollapseCase(; name="fake", kind=:exact, load=dummy_load,
        kwargs=NamedTuple(), source="test", kw...)

    @testset "check_exact" begin
        c = exact_case(exact=Dict(1 => 0.5, 2 => 4 / 3))
        @test all(r.pass for r in check_exact(c, fake([0.50001, 1.34], [1e-5, 0.01])))

        c = exact_case(exact=Dict(1 => 0.5, 2 => 4 / 3), σmax=Dict(1 => 1.3e-4, 2 => 1.8e-2))
        rows = check_exact(c, fake([0.5, 1.3775], [1.9e-4, 0.088]))
        @test length(rows) == 2
        @test !any(r.pass for r in rows)
        @test all(occursin("σmax", r.note) for r in rows)

        c = exact_case(exact=Dict(1 => 0.5, 2 => 4 / 3))
        rows = check_exact(c, fake([0.51, 4 / 3], [1e-3, 0.01]))
        @test [r.pass for r in rows] == [false, true]
        c = exact_case(exact=Dict(1 => 0.5, 2 => 4 / 3), δabs=Dict(1 => 0.01))
        @test all(r.pass for r in check_exact(c, fake([0.51, 4 / 3], [1e-3, 0.01])))
    end

    @testset "check_regression" begin
        ref = Dict{String,Any}("optimal_ps" => [1.0, 2.0], "optimal_ps_error" => [0.1, 0.1])
        @test all(r.pass for r in check_regression("x", fake([1.02, 2.0], [0.1, 0.1]), ref))
        @test !all(r.pass for r in check_regression("x", fake([1.03, 2.0], [0.1, 0.1]), ref))
        @test !all(r.pass for r in check_regression("x", fake([1.0, 2.0], [0.04, 0.1]), ref))
        @test !all(r.pass for r in check_regression("x", fake([1.0, 2.0], [0.21, 0.1]), ref))

        rows = check_regression("x", fake([1.0], [0.1]), ref)
        @test length(rows) == 1 && !rows[1].pass && occursin("parameter count", rows[1].note)

        rows = check_regression("nope", fake([1.0], [0.1]), nothing)
        @test length(rows) == 1 && !rows[1].pass && occursin("update_references.jl", rows[1].note)

        # NaN reference entries mean "no published value": skipped, never a failure
        ref = Dict{String,Any}("optimal_ps" => [NaN, 2.0], "optimal_ps_error" => [NaN, 0.1])
        rows = check_regression("x", fake([123.0, 2.0], [5.0, 0.1]), ref)
        @test all(r.pass for r in rows)
        @test rows[1].note == "no published value"
        rows = check_regression("x", fake([123.0, 2.1], [5.0, 0.1]), ref)
        @test [r.pass for r in rows] == [true, false]
    end

    @testset "merge_last" begin
        old = Dict{String,Any}("published" => Dict("a" => Dict("optimal_ps" => [1.0])),
            "last" => Dict())
        merged = merge_last(old, Dict("a" => Dict("optimal_ps" => [2.0])))
        @test merged["published"] == old["published"]
        @test merged["last"]["a"]["optimal_ps"] == [2.0]
    end

    @testset "broken cases" begin
        # a broken case that fails is reported as Broken, not as a failure; one that passes is flagged
        passing = [(param=1, expected=1.0, got=1.0, dev_sigma=0.0, err_ratio=0.1, pass=true, note="")]
        failing = [(param=1, expected=1.0, got=2.0, dev_sigma=9.0, err_ratio=0.1, pass=false, note="x")]
        @test !exact_case().broken
        @test exact_case(broken=true).broken
        ts = @testset NoThrowTestSet "failing broken" begin
            report_rows("fake", failing, "exact", true)
        end
        @test length(ts.results) == 1 && ts.results[1] isa Test.Broken
        ts = @testset NoThrowTestSet "passing broken" begin
            report_rows("fake", passing, "exact", true)
        end
        @test length(ts.results) == 1 && ts.results[1] isa Test.Error
        ts = @testset NoThrowTestSet "failing not broken" begin
            redirect_stdout(devnull) do
                report_rows("fake", failing, "exact", false)
            end
        end
        @test ts.results[1] isa Test.Fail
    end
end
