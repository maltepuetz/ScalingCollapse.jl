# Runs every collapse case (exact and regression) of the benchmark suite.
t0 = time()
@isdefined(CollapseCase) || include("framework.jl")
for f in sort(filter(f -> startswith(f, "cases_") && endswith(f, ".jl"), readdir(@__DIR__)))
    include(joinpath(@__DIR__, f))
end

@testset "Collapse benchmark cases" begin
    run_cases(CASES, load_references(joinpath(@__DIR__, "references.toml")))
end
@info "Collapse benchmark suite finished" seconds = round(time() - t0; digits=1)
