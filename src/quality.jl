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

##### Helpers shared by the quality functions #####

# true if ps lies outside p_space
function _out_of_bounds(sp, ps)
    for (i, p) in enumerate(ps)
        (p < sp.p_space[i][1] || p > sp.p_space[i][end]) && return true
    end
    return false
end

# the data of all sizes rescaled with the free parameters ps. sp.sf.f has the abstract type
# Function and ps has no static length, so they are passed through a function barrier as
# concrete types: only the call of _scale_all is dispatched at run time.
_scale_all(sp, ps) = _scale_all(sp.sf.f, sp.data, ntuple(i -> Float64(ps[i]), length(ps)))

function _scale_all(f::F, data, ps::NTuple{N,Float64}) where {F,N}
    return Data[f(d, ps...) for d in data]
end

# the window where all rescaled sizes overlap, intersected with dx
function _overlap_interval(scaled_data, dx)
    lo = max(maximum(d.xs[1] for d in scaled_data), dx[1])
    hi = min(minimum(d.xs[end] for d in scaled_data), dx[2])
    return lo, hi
end

# the points scaled_data[l] at the indices ranges[l] of all sizes, merged into one list sorted
# by x. Ties keep the order of the sizes, as a stable sort of the concatenated points would.
function _merge_sorted(scaled_data, ranges)
    m = sum(length, ranges)
    xs = Vector{Float64}(undef, m)
    ys = Vector{Float64}(undef, m)
    es = Vector{Float64}(undef, m)
    next = [first(r) for r in ranges]  # next unmerged index of every size
    for k in 1:m
        # the size with the smallest next x (the first of tied sizes)
        best = 0
        for l in eachindex(ranges)
            next[l] <= last(ranges[l]) || continue
            if best == 0 || scaled_data[l].xs[next[l]] < scaled_data[best].xs[next[best]]
                best = l
            end
        end
        d, i = scaled_data[best], next[best]
        xs[k], ys[k], es[k] = d.xs[i], d.ys[i], d.es[i]
        next[best] += 1
    end
    return xs, ys, es
end

include("quality_functions/Linear.jl")
include("quality_functions/SingleSpline.jl")
include("quality_functions/MultipleSplines.jl")
include("quality_functions/Houdayer.jl")
