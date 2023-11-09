using Scaling
using Test
using Aqua

@testset "Scaling.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(Scaling)
    end
    # Write your tests here.
end
