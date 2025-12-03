"""
    Data(L::Int, xs::Vector{Float64}, ys::Vector{Float64}, es::Vector{Float64})

A Data object stores the data for a single system size. Sorted by the x values.

# Fields
- `L::Int`: system size
- `xs::Vector{Float64}`: x values
- `ys::Vector{Float64}`: y values
- `es::Vector{Float64}`: y error values
"""
struct Data
    L::Int               # system size
    xs::Vector{Float64}  # x values
    ys::Vector{Float64}  # y values
    es::Vector{Float64}  # y error values
    function Data(L, xs, ys, es)
        # sort data by x values
        @assert length(xs) == length(ys) == length(es)
        perm = sortperm(xs)
        xsd = similar(xs)
        ysd = similar(ys)
        esd = similar(es)
        for i in eachindex(perm)
            p = perm[i]
            xsd[i] = xs[p]
            ysd[i] = ys[p]
            esd[i] = abs(es[p])
        end
        return new(L, xsd, ysd, esd)
    end
end

function Data(L, xs, ys)
    return Data(L, xs, ys, zeros(eltype(ys), length(ys)))
end

function Data(L, xs, ys::AbstractVector{M}) where {M<:Measurement}
    return Data(L, xs, Measurements.value.(ys), Measurements.uncertainty.(ys))
end

# overload method to check equality of two Data objects
function Base.isequal(d1::Data, d2::Data)
    return d1.L == d2.L && d1.xs == d2.xs && d1.ys == d2.ys && d1.es == d2.es
end

function _unzip_data_aux(ys::Matrix, ordering_length::Int)
    nrows, ncols = size(ys)
    (nrows == ncols) && (@warn "Data orientation is ambiguous, selecting column-first ordering")
    (ncols == ordering_length) && return eachcol(ys)
    (nrows == ordering_length) && return eachrow(ys)
    error("Unable to determine data orientation")
end

function _unzip_data_aux(xs::Vector{Vector{T}}, ordering_length::Int) where {T}
    (length(xs) == ordering_length) && return (x for x in xs)
    error("Expected input of length $(ordering_length), found $(length(xs)) instead")
end

function _unzip_data_aux(xs::Vector{T}, ordering_length::Int) where {T}
    return Iterators.repeated(xs, ordering_length)
end

function _unzip_data(xs, ys, Ls)
    map((L, x, y) -> Data(L, x, y), Ls, xs, ys)
end

function _unzip_data(xs, ys, es, Ls)
    map((L, x, y, e) -> Data(L, x, y, e), Ls, xs, ys, es)
end

function unzip_data(xs::Union{Vector,Matrix}, ys::Union{Vector,Matrix}, Ls::Vector{<:Integer})
    Lcount = length(Ls)
    xs_reordered = _unzip_data_aux(xs, Lcount)
    ys_reordered = _unzip_data_aux(ys, Lcount)
    _unzip_data(xs_reordered, ys_reordered, Ls)
end

function unzip_data(xs::Union{Vector,Matrix}, ys::Union{Vector,Matrix}, es::Union{Vector,Matrix}, Ls::Vector{<:Integer})
    Lcount = length(Ls)
    xs_reordered = _unzip_data_aux(xs, Lcount)
    ys_reordered = _unzip_data_aux(ys, Lcount)
    es_reordered = _unzip_data_aux(es, Lcount)
    _unzip_data(xs_reordered, ys_reordered, es_reordered, Ls)
end
