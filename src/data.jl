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
        perm = sortperm(xs)
        xs = xs[perm]
        ys = ys[perm]
        es = es[perm] .|> abs
        return new(L, xs, ys, es)
    end
end


# overload method to check equality of two Data objects
function Base.isequal(d1::Data, d2::Data)
    return d1.L == d2.L && d1.xs == d2.xs && d1.ys == d2.ys && d1.es == d2.es
end

# unzip data from arguments
function unzip_data(
    xs::Vector{T1},
    ys::Array{T2,2},
    Ls::Vector{T4}
) where {T1,T2<:Real,T4<:Integer}

    # get ordering of ys and check for consistency
    correct_ordering = true
    size(ys, 1) == length(Ls) && (correct_ordering = false)
    @assert correct_ordering ? length(xs) == size(ys, 1) : length(xs) == size(ys, 2)
    @assert correct_ordering ? length(Ls) == size(ys, 2) : length(Ls) == size(ys, 1)

    # unzip data
    data = Vector{Data}(undef, length(Ls))
    for (l, L) in enumerate(Ls)
        if correct_ordering
            data[l] = Data(L, xs, ys[:, l], zeros(Float64, length(xs)))
        else
            data[l] = Data(L, xs, ys[l, :], zeros(Float64, length(xs)))
        end
    end

    return data
end

function unzip_data(
    xs::Vector{T1},
    ys::Array{T2,2},
    es::Array{T3,2},
    Ls::Vector{T4}
) where {T1,T2,T3<:Real,T4<:Integer}

    # get ordering of ys and check for consistency
    correct_ordering = true
    size(ys, 1) == length(Ls) && (correct_ordering = false)
    @assert correct_ordering ? length(xs) == size(ys, 1) : length(xs) == size(ys, 2)
    @assert correct_ordering ? length(Ls) == size(ys, 2) : length(Ls) == size(ys, 1)

    # unzip data
    data = Vector{Data}(undef, length(Ls))
    for (l, L) in enumerate(Ls)
        if correct_ordering
            data[l] = Data(L, xs, ys[:, l], es[:, l])
        else
            data[l] = Data(L, xs, ys[l, :], es[l, :])
        end
    end

    return data
end

function unzip_data(
    xs::Array{T1,2},
    ys::Array{T2,2},
    Ls::Vector{T4}
) where {T1,T2<:Real,T4<:Integer}

    # get ordering of ys and check for consistency
    correct_ordering = true
    size(ys, 1) == length(Ls) && (correct_ordering = false)
    @assert correct_ordering ? length(Ls) == size(ys, 2) : length(Ls) == size(ys, 1)

    # unzip data
    data = Vector{Data}(undef, length(Ls))
    for (l, L) in enumerate(Ls)
        if correct_ordering
            data[l] = Data(L, xs[:, l], ys[:, l], zeros(Float64, size(xs, 1)))
        else
            data[l] = Data(
                L,
                xs[l, :],
                ys[l, :],
                zeros(Float64, size(ys, correct_ordering ? 1 : 2))
            )
        end
    end

    return data
end

function unzip_data(
    xs::Array{T1,2},
    ys::Array{T2,2},
    es::Array{T3,2},
    Ls::Vector{T4}
) where {T1,T2,T3<:Real,T4<:Integer}

    # get ordering of ys and check for consistency
    correct_ordering = true
    size(ys, 1) != length(Ls) ? nothing : (correct_ordering == false)
    @assert size(xs) == size(ys)
    @assert correct_ordering ? length(Ls) == size(ys, 2) : length(Ls) == size(ys, 1)#
    @assert size(ys) == size(es)

    # unzip data
    data = Vector{Data}(undef, length(Ls))
    for (l, L) in enumerate(Ls)
        if correct_ordering
            data[l] = Data(L, xs[:, l], ys[:, l], es[:, l])
        else
            data[l] = Data(L, xs[l, :], ys[l, :], es[l, :])
        end
    end

    return data
end

function unzip_data(
    xs::Vector{Vector{T1}},
    ys::Vector{Vector{T2}},
    Ls::Vector{T4}
) where {T1,T2<:Real,T4<:Integer}

    # check for consistency
    @assert length(xs) == length(ys)
    @assert length(xs) == length(Ls)

    # unzip data
    data = Vector{Data}(undef, length(Ls))
    for (l, L) in enumerate(Ls)
        data[l] = Data(L, xs[l], ys[l], zeros(Float64, length(xs[l])))
    end

    return data
end

function unzip_data(
    xs::Vector{Vector{T1}},
    ys::Vector{Vector{T2}},
    es::Vector{Vector{T3}},
    Ls::Vector{T4}
) where {T1,T2,T3<:Real,T4<:Integer}

    # check for consistency
    @assert length(xs) == length(ys)
    @assert length(xs) == length(Ls)
    @assert length(xs) == length(es)

    # unzip data
    data = Vector{Data}(undef, length(Ls))
    for (l, L) in enumerate(Ls)
        data[l] = Data(L, xs[l], ys[l], es[l])
    end

    return data
end

# unzip Measurements data
function unzip_data(
    xs::Vector{T1},
    ys::Array{Measurement{T2},2},
    Ls::Vector{T4}
) where {T1,T2<:Real,T4<:Integer}

    # get ordering of ys and check for consistency
    correct_ordering = true
    size(ys, 1) != length(Ls) ? nothing : (correct_ordering == false)
    @assert correct_ordering ? length(xs) == size(ys, 1) : length(xs) == size(ys, 2)
    @assert correct_ordering ? length(Ls) == size(ys, 2) : length(Ls) == size(ys, 2)

    # unzip data
    data = Vector{Data}(undef, length(Ls))
    for (l, L) in enumerate(Ls)
        if correct_ordering
            data[l] = Data(
                L,
                xs,
                Measurements.value.(ys[:, l]),
                Measurements.uncertainty.(ys[:, l])
            )
        else
            data[l] = Data(
                L,
                xs,
                Measurements.value.(ys[l, :]),
                Measurements.uncertainty.(ys[l, :])
            )
        end
    end

    return data
end

function unzip_data(
    xs::Array{T1,2},
    ys::Array{Measurement{T2},2},
    Ls::Vector{T4}
) where {T1,T2<:Real,T4<:Integer}

    # get ordering of ys and check for consistency
    correct_ordering = true
    size(ys, 1) != length(Ls) ? nothing : (correct_ordering == false)
    @assert size(xs) == size(ys)
    @assert correct_ordering ? length(Ls) == size(ys, 2) : length(Ls) == size(ys, 1)

    # unzip data
    data = Vector{Data}(undef, length(Ls))
    for (l, L) in enumerate(Ls)
        if correct_ordering
            data[l] = Data(
                L,
                xs[:, l],
                Measurements.value.(ys[:, l]),
                Measurements.uncertainty.(ys[:, l])
            )
        else
            data[l] = Data(
                L,
                xs[l, :],
                Measurements.value.(ys[l, :]),
                Measurements.uncertainty.(ys[l, :])
            )
        end
    end

    return data
end

function unzip_data(
    xs::Vector{Vector{T1}},
    ys::Vector{Vector{Measurement{T2}}},
    Ls::Vector{T4}
) where {T1,T2<:Real,T4<:Integer}

    # check for consistency
    @assert length(xs) == length(ys)
    @assert length(xs) == length(Ls)

    # unzip data
    data = Vector{Data}(undef, length(Ls))
    for (l, L) in enumerate(Ls)
        data[l] = Data(
            L,
            xs[l],
            Measurements.value.(ys[l]),
            Measurements.uncertainty.(ys[l])
        )
    end

    return data

end
