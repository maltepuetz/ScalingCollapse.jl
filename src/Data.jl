"""
    Data(L::Int, xs::Vector{Float64}, ys::Vector{Float64}, es::Vector{Float64})

A Data object stores the data for a single system size.

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
end


# unzip data from arguments
function unzip_data(xs::Vector{T1}, ys::Array{T2,2}, Ls::Vector{Int}) where {T1,T2<:Real}

    # get ordering of ys and check for consistency
    correct_ordering = true
    size(ys, 1) != length(Ls) ? nothing : (correct_ordering == false)
    @assert correct_ordering ? length(xs) == size(ys, 1) : length(xs) == size(ys, 2)
    @assert correct_ordering ? length(Ls) == size(ys, 2) : length(Ls) == size(ys, 2)

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
    ys::Vector{T2},
    es::Vector{T3},
    Ls::Vector{Int}
) where {T1,T2,T3<:Real}

    # get ordering of ys and check for consistency
    correct_ordering = true
    size(ys, 1) != length(Ls) ? nothing : (correct_ordering == false)
    @assert correct_ordering ? length(xs) == size(ys, 1) : length(xs) == size(ys, 2)
    @assert correct_ordering ? length(Ls) == size(ys, 2) : length(Ls) == size(ys, 2)
    @assert size(ys) == size(es)

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

function unzip_data(xs::Array{T1,2}, ys::Array{T2,2}, Ls::Vector{Int}) where {T1,T2<:Real}

    # get ordering of ys and check for consistency
    correct_ordering = true
    size(ys, 1) != length(Ls) ? nothing : (correct_ordering == false)
    @assert size(xs) == size(ys)
    @assert correct_ordering ? length(Ls) == size(ys, 2) : length(Ls) == size(ys, 1)

    # unzip data
    data = Vector{Data}(undef, length(Ls))
    for (l, L) in enumerate(Ls)
        if correct_ordering
            data[l] = Data(L, xs[:, l], ys[:, l], zeros(Float64, size(xs, 1)))
        else
            data[l] = Data(L, xs[l, :], ys[l, :], zeros(Float64, size(xs, 1)))
        end
    end

    return data
end

function unzip_data(
    xs::Array{T1,2},
    ys::Array{T2,2},
    es::Array{T3,2},
    Ls::Vector{Int}
) where {T1,T2,T3<:Real}

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

function unzip_data(xs::Vector{Vector}, ys::Vector{Vector}, Ls::Vector{Int})

    # check for consistency
    @assert length(xs) == length(ys)
    @assert length(xs) == length(Ls)

    # unzip data
    data = Vector{Data}(undef, length(Ls))
    for (l, L) in enumerate(Ls)
        data[l] = Data(L, xs[l], ys[l], zeros(Float64, length(xs)))
    end

    return data

end

function unzip_data(
    xs::Vector{Vector},
    ys::Vector{Vector},
    es::Vector{Vector},
    Ls::Vector{Int})

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
