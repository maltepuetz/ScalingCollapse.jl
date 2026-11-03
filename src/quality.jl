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

##### Workspace of the quality functions #####

# Buffers of the quality evaluations of one problem, so that the built-in quality functions do
# not allocate when they are called repeatedly (by the scan, the fits and the error analysis).
# Its type carries the number N of free parameters and the types of the scaling function, which
# are abstract in ScalingProblem: code that receives a Workspace is compiled for them. Every
# task that evaluates the quality needs its own workspace.
mutable struct Workspace{N,F,C}
    const f::F                      # sp.sf.f
    const coefficients::C           # sp.sf.coefficients: rescale in place (nothing: call f)
    const lower::Vector{Float64}    # bounds of p_space
    const upper::Vector{Float64}
    const data::Vector{Data}        # the rescaled data of every size
    yvals::Matrix{Float64}          # MultipleSplines: interpolated y (sizes × points) ...
    evals::Matrix{Float64}          # ... and errors
    const slopes::Vector{Float64}   # pchip slopes of one size
    const densities::Vector{Float64}  # data point densities of the sizes
    const xs::Vector{Float64}       # Linear: the points of all sizes, merged
    const ys::Vector{Float64}
    const es::Vector{Float64}
    const next::Vector{Int}         # next unmerged index of every size
    const ranges::Vector{UnitRange{Int}}  # indices of the points in the window, per size
end

function Workspace(sp)
    sf, data = sp.sf, sp.data
    n_points = sum(d -> length(d.xs), data)
    # buffers for the presets; a custom f replaces the elements with the Data it returns
    scaled = [Data(d.L, similar(d.xs), similar(d.ys), similar(d.es), Val(:sorted)) for d in data]
    return Workspace{n_parameters(sf),typeof(sf.f),typeof(sf.coefficients)}(
        sf.f,
        sf.coefficients,
        Float64[first(p) for p in sp.p_space],
        Float64[p[end] for p in sp.p_space],
        scaled,
        Matrix{Float64}(undef, 0, 0),
        Matrix{Float64}(undef, 0, 0),
        Vector{Float64}(undef, maximum(d -> length(d.xs), data)),
        Vector{Float64}(undef, length(data)),
        Vector{Float64}(undef, n_points),
        Vector{Float64}(undef, n_points),
        Vector{Float64}(undef, n_points),
        Vector{Int}(undef, length(data)),
        Vector{UnitRange{Int}}(undef, length(data)),
    )
end

# the workspace of sp; nothing for an object that only mimics a ScalingProblem (in the tests),
# whose quality function is then called without one
_workspace(sp) = nothing

# whether to split the work on sp into tasks; never for an object that only mimics a
# ScalingProblem (in the tests)
_threaded(sp) = false

# [f(state, x) for x in xs], where state = init(). If threaded, xs is split into at most
# Threads.nthreads() contiguous chunks that are worked through in parallel tasks, each with its
# own state (e.g. a Workspace), so f(state, x) must not share anything else between tasks. An
# exception in a task is rethrown as it is.
function _tmap(f, init, xs, threaded::Bool)
    n_tasks = threaded ? min(Threads.nthreads(), length(xs)) : 1
    n_tasks <= 1 && return _map_chunk(f, init(), xs)
    chunks = Iterators.partition(xs, cld(length(xs), n_tasks))
    tasks = [Threads.@spawn _map_chunk(f, init(), chunk) for chunk in chunks]
    return reduce(vcat, map(_fetch, tasks))
end

# a function barrier for the type of state
_map_chunk(f::F, state, xs) where {F} = [f(state, x) for x in xs]

function _fetch(task)
    try
        return fetch(task)
    catch e
        e isa TaskFailedException && throw(e.task.exception)
        rethrow()
    end
end

# A quality function without a method for a workspace (e.g. a custom one) is called without it
(qf::QualityFunction)(::Workspace, sp, ps; check_bounds=false) = qf(sp, ps; check_bounds)

# ps -> quality(ws, sp, ps; check_bounds) behind a function barrier: the closure has concrete
# types, so a loop or an optimizer that receives it as an argument is compiled for them and
# does not allocate. check_bounds=nothing omits the keyword (false for a QualityFunction).
function _quality_function(sp, quality, ws; check_bounds::Union{Nothing,Bool}=nothing)
    return _quality_function(quality, ws, sp, check_bounds)
end
function _quality_function(quality::Q, ws::W, sp, check_bounds::B) where {Q,W,B}
    return ps -> _call_quality(quality, ws, sp, ps, check_bounds)
end

# only a QualityFunction takes a workspace; any other callable (and any call without a
# workspace) gets quality(sp, ps; check_bounds), the keyword omitted for nothing
function _call_quality(quality, ws, sp, ps, check_bounds)
    return isnothing(check_bounds) ? quality(sp, ps) : quality(sp, ps; check_bounds)
end
function _call_quality(quality::QualityFunction, ws::Workspace, sp, ps, check_bounds)
    return quality(ws, sp, ps; check_bounds=something(check_bounds, false))
end

# the interpolation buffers of ws with n sizes and N_steps points
function _interpolation_buffers!(ws::Workspace, n, N_steps)
    if size(ws.yvals) != (n, N_steps)
        ws.yvals = Matrix{Float64}(undef, n, N_steps)
        ws.evals = Matrix{Float64}(undef, n, N_steps)
    end
    return ws.yvals, ws.evals
end

##### Helpers shared by the quality functions #####

# true if ps lies outside p_space
function _out_of_bounds(ws::Workspace, ps)
    for (i, p) in enumerate(ps)
        (p < ws.lower[i] || p > ws.upper[i]) && return true
    end
    return false
end

# the data of all sizes rescaled with the free parameters ps, in ws.data
function _scale_all!(ws::Workspace{N}, data, ps) where {N}
    return _scale_all!(ws.coefficients, ws, data, ntuple(i -> Float64(ps[i]), Val(N)))
end

# a preset: rescale into the buffers (::C makes Julia 1.10 specialize on coefficients)
function _scale_all!(coefficients::C, ws::Workspace, data, ps::Tuple) where {C}
    for l in eachindex(data)
        d = data[l]
        _linear_scaling!(ws.data[l], d, coefficients(d.L, ps...)...)
    end
    return ws.data
end

# a custom scaling function: store the Data it returns
function _scale_all!(::Nothing, ws::Workspace, data, ps::Tuple)
    for l in eachindex(data)
        ws.data[l] = ws.f(data[l], ps...)
    end
    return ws.data
end

# the window where all rescaled sizes overlap, intersected with dx
function _overlap_interval(scaled_data, dx)
    lo = max(maximum(d.xs[1] for d in scaled_data), dx[1])
    hi = min(minimum(d.xs[end] for d in scaled_data), dx[2])
    return lo, hi
end

# the points scaled_data[l] at the indices ranges[l] of all sizes, merged into the first
# entries of xs, ys and es, sorted by x; returns their number. Ties keep the order of the
# sizes, as a stable sort of the concatenated points would. next is a buffer of one index per
# size.
function _merge_sorted!(xs, ys, es, next, scaled_data, ranges)
    m = sum(length, ranges)
    for l in eachindex(ranges)
        next[l] = first(ranges[l])
    end
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
    return m
end

include("quality_functions/Linear.jl")
include("quality_functions/SingleSpline.jl")
include("quality_functions/MultipleSplines.jl")
include("quality_functions/Houdayer.jl")
