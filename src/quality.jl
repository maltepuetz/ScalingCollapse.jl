abstract type QualityFunction end

function Base.show(io::IO, qf::QualityFunction)
    print(io, typeof(qf), " <: QualityFunction")
end

include("quality_functions/Linear.jl")
include("quality_functions/SingleSpline.jl")
include("quality_functions/MultipleSplines.jl")
include("quality_functions/Houdayer.jl")
