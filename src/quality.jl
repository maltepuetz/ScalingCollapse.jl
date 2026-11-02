"""
    ScalingCollapse.QualityFunction

Abstract supertype of the quality functions, which measure how well the rescaled data of all
system sizes collapse onto one curve: [`MultipleSplines`](@ref) (the default),
[`Linear`](@ref), [`SingleSpline`](@ref) and [`Houdayer`](@ref). Lower values mean a better
collapse.

A quality function is a callable struct `qf(sp::ScalingProblem, ps; check_bounds=false)`
that returns a `Float64` for the free parameters `ps`. It rescales every `Data` in
`sp.data` with `sp.sf.f(d, ps...)`, compares the curves within the overlap of all sizes
intersected with `sp.dx`, and returns `Inf` if that interval is empty, or if
`check_bounds=true` and `ps` lies outside `sp.p_space`. A new quality function can be
defined in the same way and passed as `quality` to [`ScalingProblem`](@ref).
"""
abstract type QualityFunction end

function Base.show(io::IO, qf::QualityFunction)
    print(io, typeof(qf), " <: QualityFunction")
end

include("quality_functions/Linear.jl")
include("quality_functions/SingleSpline.jl")
include("quality_functions/MultipleSplines.jl")
include("quality_functions/Houdayer.jl")
