function quality_houdayer(sp, parameters; check_bounds=false)

    # check if parameters are in bounds -- if not return Inf
    if check_bounds
        for (i, p) in enumerate(parameters)
            if p < sp.p_space[i][1] || p > sp.p_space[i][end]
                return Inf
            end
        end
    end

    scaled_data = [sp.sf.f(sp.data[i], parameters...) for i in eachindex(sp.data)]

    # set interval in which we want to optimize
    interval = Vector{Float64}(undef, 2)
    interval[1] = max(
        maximum(scaled_data[i].xs[1] for i in eachindex(scaled_data)),
        sp.dx[1]
    )
    interval[2] = min(
        minimum(scaled_data[i].xs[end] for i in eachindex(scaled_data)),
        sp.dx[2]
    )

    # check whether the constructed optimization interval is valid
    interval[1] > interval[2] && return Inf

    N = 0  # counts number of terms in the sum
    S = 0.0  # quality of the scaling

    for (i, d) in enumerate(scaled_data)

        j = searchsortedfirst(d.xs, interval[1])  # index of first element in interval

        while (j <= length(d.xs)) && (d.xs[j] <= interval[2])

            Y_ij, dY2_ij = fit_mastercurve(scaled_data, i, j)
            Y_ij == Inf && (j += 1; continue)  # no contribution

            S += (d.ys[j] - Y_ij)^2 / (d.es[j]^2 + dY2_ij)

            j += 1
            N += 1
        end
    end
    S /= N

    return S
end

function fit_mastercurve(scaled_data, i, j)

    # set arrays to store set of data used for fitting
    xs = Float64[]
    ys = Float64[]
    ws = Float64[]

    for (i_prime, d_prime) in enumerate(scaled_data)
        i_prime == i && continue
        j_prime = searchsortedlast(d_prime.xs, scaled_data[i].xs[j])
        ((j_prime == 0) || (j_prime == length(d_prime.xs))) && continue  # no contribution
        push!(xs, d_prime.xs[j_prime])
        push!(xs, d_prime.xs[j_prime+1])
        push!(ys, d_prime.ys[j_prime])
        push!(ys, d_prime.ys[j_prime+1])
        push!(ws, d_prime.es[j_prime])
        push!(ws, d_prime.es[j_prime+1])
    end

    length(xs) == 0 && return Inf, Inf  # no contribution

    # check whether all elements in ws are zeros
    all(ws .== 0.0) && (ws .= 1.0)
    ws .= 1 ./ ws .^ 2

    # if there are single zero errors, set weigth to max(ws \ Inf)
    for i in eachindex(ws)
        isinf(ws[i]) && (ws[i] = maximum(ws[isfinite.(ws)]))
    end

    K = sum(ws)
    Kx = sum(ws .* xs)
    Ky = sum(ws .* ys)
    Kxx = sum(ws .* xs .^ 2)
    Kxy = sum(ws .* xs .* ys)
    Delta = K * Kxx - Kx^2

    a = (Kxx * Ky - Kx * Kxy) / Delta
    b = (K * Kxy - Kx * Ky) / Delta

    Y_ij = a + scaled_data[i].xs[j] * b
    dY2_ij = 1 / Delta * (Kxx - 2 * scaled_data[i].xs[j] * Kx + scaled_data[i].xs[j]^2 * K)

    return Y_ij, dY2_ij
end
