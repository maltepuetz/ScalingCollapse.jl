
"""
    Linear()

Quality function that interpolates all scaled data points to a single master curve. For
every data point, we take the closest data points, to the left and to the right to
interpolate a linear function. The quality is then calculated as the sum of the squared
differences between the data points and the interpolated values.
"""
struct Linear <: QualityFunction end



##### Linear Quality Function #####
function (::Linear)(sp, parameters; check_bounds=false)

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

    # collect all data points in the interval
    xs = Float64[]
    ys = Float64[]
    es = Float64[]
    for d in scaled_data
        i_l = max(1, searchsortedlast(d.xs, interval[1]) + 1)
        i_r = min(length(d.xs), searchsortedfirst(d.xs, interval[2]) - 1)
        for i in i_l:i_r
            push!(xs, d.xs[i])
            push!(ys, d.ys[i])
            push!(es, d.es[i])
        end
    end

    # sort data points
    perm = sortperm(xs)
    xs = xs[perm]
    ys = ys[perm]
    es = es[perm]

    # calculate quality of the scaling
    S = 0.0
    if sp.errors_defined
        for i in eachindex(xs)[2:end-1]
            # interpolate linearly
            ybar = ((xs[i+1] - xs[i]) * ys[i-1] - (xs[i-1] - xs[i]) * ys[i+1]) / (xs[i+1] - xs[i-1])
            expected_var = es[i]^2 + ((xs[i+1] - xs[i]) / (xs[i+1] - xs[i-1]))^2 * es[i-1]^2 + ((xs[i-1] - xs[i]) / (xs[i+1] - xs[i-1]))^2 * es[i+1]^2
            S += (ys[i] - ybar)^2 / expected_var
        end
    else
        for i in eachindex(xs)[2:end-1]
            # interpolate linearly
            ybar = ((xs[i+1] - xs[i]) * ys[i-1] - (xs[i-1] - xs[i]) * ys[i+1]) / (xs[i+1] - xs[i-1])
            S += (ys[i] - ybar)^2
        end
    end
    S /= (length(xs) - 2)
    return S
end
