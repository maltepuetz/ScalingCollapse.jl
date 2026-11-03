"""
    Houdayer()

Quality function of Houdayer & Hartmann, Phys. Rev. B 70, 014418 (2004). It is biased on
coarse data grids (see below) and therefore not recommended.

# How it works
Following Houdayer & Hartmann, Phys. Rev. B 70, 014418 (2004), for every scaled data point
`(x_ij, y_ij)` of system size `i` in the overlap interval, the master curve is approximated
by a weighted linear least-squares fit to the two data points of every other system size
that bracket `x_ij`. The overlap interval is the window covered by all sizes, intersected
with `dx`. With the fit value `Y_ij` and its variance `dY²_ij`, the quality is

    S = 1/N Σ_ij (y_ij - Y_ij)² / (e_ij² + dY²_ij).

If no error bars are given, all weights are set to 1.

# Known bias
The linear fit misses the curvature of the master curve between the bracketing points. The
error this causes grows with the square of the spacing of the scaled data points, which is
largest for the biggest system size. A larger `ν` compresses the scaled x axis, shrinks this
spacing, and thereby lowers `S`, so on coarse data grids the fitted `ν` comes out too
large. For example, on noise-free synthetic data with 4 sizes and 41 points each (only
about 5 of the largest size lie in the interval), `ν` came out 2-3% too large; with 161
points per size, 0.1%.
The bias vanishes as the grid is refined. [`Linear`](@ref) and [`MultipleSplines`](@ref) do
not show it.

# Ruggedness
The bracketing data points change discontinuously as scaled points of different sizes pass
each other, so the quality jumps and has many shallow local minima on uniformly gridded data,
like `Linear` (see its docstring). `MultipleSplines` gives a smooth quality.
"""
struct Houdayer <: QualityFunction end

##### Houdayer Quality Function #####
function (qf::Houdayer)(sp, parameters; check_bounds=false)
    return qf(Workspace(sp), sp, parameters; check_bounds)
end

function (qf::Houdayer)(ws::Workspace, sp, parameters; check_bounds=false)

    # check if parameters are in bounds -- if not return Inf
    check_bounds && _out_of_bounds(ws, parameters) && return Inf

    scaled_data = _scale_all!(ws, sp.data, parameters)

    # set interval in which we want to optimize
    interval = _overlap_interval(scaled_data, sp.dx)

    # check whether the constructed optimization interval is valid (also rejects NaN)
    interval[1] <= interval[2] || return Inf

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
    N == 0 && return Inf
    S /= N

    return S
end

# the index j' of the data point of d just left of x, if d brackets x with j' and j' + 1,
# else 0
function _bracket(d, x)
    j_prime = searchsortedlast(d.xs, x)
    return j_prime == length(d.xs) ? 0 : j_prime
end

# weighted linear least-squares fit through the two data points of every other size that
# bracket scaled_data[i].xs[j]; the weights are 1/e² (all 1 if all errors are zero, and the
# largest finite weight for single zero errors). Returns the fit value and its variance at x.
function fit_mastercurve(scaled_data, i, j)

    x = scaled_data[i].xs[j]

    # first pass: whether any point contributes, and the largest finite weight
    any_point = false
    all_zero = true
    w_max = 0.0
    for (i_prime, d_prime) in enumerate(scaled_data)
        i_prime == i && continue
        j_prime = _bracket(d_prime, x)
        j_prime == 0 && continue  # no contribution
        any_point = true
        for e in (d_prime.es[j_prime], d_prime.es[j_prime+1])
            all_zero &= e == 0.0
            w = 1 / e^2
            isfinite(w) && (w_max = max(w_max, w))
        end
    end

    any_point || return Inf, Inf  # no contribution

    # second pass: the weighted sums
    K = Kx = Ky = Kxx = Kxy = 0.0
    for (i_prime, d_prime) in enumerate(scaled_data)
        i_prime == i && continue
        j_prime = _bracket(d_prime, x)
        j_prime == 0 && continue
        for k in (j_prime, j_prime + 1)
            x_k, y_k = d_prime.xs[k], d_prime.ys[k]
            w = all_zero ? 1.0 : 1 / d_prime.es[k]^2
            isinf(w) && (w = w_max)
            K += w
            Kx += w * x_k
            Ky += w * y_k
            Kxx += w * x_k^2
            Kxy += w * x_k * y_k
        end
    end
    Delta = K * Kxx - Kx^2

    a = (Kxx * Ky - Kx * Kxy) / Delta
    b = (K * Kxy - Kx * Ky) / Delta

    Y_ij = a + x * b
    dY2_ij = 1 / Delta * (Kxx - 2 * x * Kx + x^2 * K)

    return Y_ij, dY2_ij
end
