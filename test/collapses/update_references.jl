# Rewrites the [last] table of references.toml from the current code. Never touches [published].
# Needs the test dependencies (HDF5, TOML), which the root project lacks. From the repo root, with
# TestEnv.jl installed in the default environment:
#   julia --project=. -e 'using TestEnv; TestEnv.activate(); include("test/collapses/update_references.jl")'
# or in a temporary environment:
#   julia -e 'using Pkg; Pkg.activate(; temp=true); Pkg.develop(path="."); Pkg.add(["HDF5", "TOML"]); include("test/collapses/update_references.jl")'
include("framework.jl")
for f in sort(filter(f -> startswith(f, "cases_") && endswith(f, ".jl"), readdir(@__DIR__)))
    include(joinpath(@__DIR__, f))
end

path = joinpath(@__DIR__, "references.toml")
old = load_references(path)
new_last = Dict{String,Any}()
for c in CASES
    c.kind === :regression || continue
    sp = run_case(c)
    new_last[c.name] = Dict{String,Any}(
        "optimal_ps" => collect(Float64, sp.optimal_ps),
        "optimal_ps_error" => collect(Float64, sp.optimal_ps_error),
        "minimum" => Float64(sp.minimum),
    )
    @info "updated" c.name sp.optimal_ps
end
open(path, "w") do io
    TOML.print(io, merge_last(old, new_last); sorted=true)
end
