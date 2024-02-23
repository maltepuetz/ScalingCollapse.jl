using ScalingCollapse
using Test
using Aqua
using HDF5

@testset "ScalingCollapse.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(ScalingCollapse; ambiguities=false,)
    end
    include("data.jl")
    include("ising.jl")
    include("cdw.jl")
    include("scalingfunction.jl")
end
