using Scaling
using Test
using Aqua
using HDF5

@testset "Scaling.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(Scaling; ambiguities=false,)
    end
    include("data.jl")
    include("ising.jl")
    include("cdw.jl")
end
