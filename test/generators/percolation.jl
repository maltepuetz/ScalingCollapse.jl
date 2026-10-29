# Newman–Ziff bond percolation on the L × L square lattice with open boundaries. Writes
# test/testdata/percolation/bond_square.h5. Not run in CI.
#
# Each sample adds the 2L(L − 1) bonds in random order, merging clusters with union-find (path
# compression, union by size), and records after every bond n whether a cluster connects the top
# and bottom rows (spanning, 0/1) and the size of the largest cluster divided by L² (strength).
# Each sample's microcanonical curve Q(n) is turned into its canonical value
# Q(p) = Σₙ B(n; M, p) Q(n), with binomial weights B computed in log-space (by the ratio
# recursion from the mode, truncated where B < e^{-60} of the mode, then normalised). The stored
# value is the mean over samples, and the error is the standard error of the per-sample canonical
# values (sample standard deviation / √nsamples).
#
# Exact values: p_c = 1/2, ν = 4/3, β = 5/36.
#
# Layout: `Ls`, `ps`, `<obs>/L{L}/{x, y, e}` for obs ∈ (spanning, strength); the attribute
# `nsamples` of `<obs>/L{L}` holds the number of samples.
#
# Usage (in an environment with HDF5):
#     julia percolation.jl [out.h5]

using HDF5
using Random

const SIZES = [16, 32, 64, 128, 256]
const NSAMPLES = Dict(16 => 10_000, 32 => 10_000, 64 => 10_000, 128 => 10_000, 256 => 2_000)
const PS = collect(0.45:0.001:0.55)

# bonds as (site, site) pairs; sites are 1:L² in row-major order, row 1 is the top row
function bonds(L)
    site(r, c) = (r - 1) * L + c
    out = Tuple{Int,Int}[]
    for r in 1:L, c in 1:L
        c < L && push!(out, (site(r, c), site(r, c + 1)))
        r < L && push!(out, (site(r, c), site(r + 1, c)))
    end
    return out
end

function findroot!(parent, i)
    root = i
    while parent[root] != root
        root = parent[root]
    end
    while parent[i] != root  # path compression
        parent[i], i = root, parent[i]
    end
    return root
end

"""
One Newman–Ziff sample. Returns the number of bonds at which top and bottom first connect and
fills `largest[n + 1]` with the largest cluster size after `n` bonds.
"""
function sample!(largest, rng, L, bondlist, order, parent, csize, top, bottom)
    N = L^2
    for i in 1:N
        parent[i], csize[i] = i, 1
        top[i], bottom[i] = i <= L, i > N - L
    end
    shuffle!(rng, order)
    nspan = typemax(Int)
    maxsize = 1
    largest[1] = 1
    for (n, b) in enumerate(order)
        a, c = findroot!(parent, bondlist[b][1]), findroot!(parent, bondlist[b][2])
        if a != c
            csize[a] < csize[c] && ((a, c) = (c, a))  # union by size: attach c to a
            parent[c] = a
            csize[a] += csize[c]
            top[a] |= top[c]
            bottom[a] |= bottom[c]
            maxsize = max(maxsize, csize[a])
            nspan == typemax(Int) && top[a] && bottom[a] && (nspan = n)
        end
        largest[n + 1] = maxsize
    end
    return nspan
end

"Normalised binomial weights B(n; M, p) for n ∈ range, computed in log-space."
function binomial_weights(M, p)
    mode = clamp(floor(Int, (M + 1) * p), 0, M)
    lr = log(p / (1 - p))
    logw = Dict(mode => 0.0)
    lo, hi = mode, mode
    while lo > 0  # log B(n − 1) − log B(n) = log(n / (M − n + 1)) − log(p / (1 − p))
        l = logw[lo] + log(lo / (M - lo + 1)) - lr
        l < -60 && break
        lo -= 1
        logw[lo] = l
    end
    while hi < M
        l = logw[hi] + log((M - hi) / (hi + 1)) + lr
        l < -60 && break
        hi += 1
        logw[hi] = l
    end
    w = [exp(logw[n]) for n in lo:hi]
    return lo:hi, w ./ sum(w)
end

function simulate(rng, L, nsamples)
    bondlist = bonds(L)
    M, N = length(bondlist), L^2
    order = collect(1:M)
    parent, csize = zeros(Int, N), zeros(Int, N)
    top, bottom = falses(N), falses(N)
    largest = zeros(Int, M + 1)
    weights = [binomial_weights(M, p) for p in PS]
    # running sums of the per-sample canonical values and their squares
    s_span, s2_span = zeros(length(PS)), zeros(length(PS))
    s_str, s2_str = zeros(length(PS)), zeros(length(PS))
    for _ in 1:nsamples
        nspan = sample!(largest, rng, L, bondlist, order, parent, csize, top, bottom)
        for (k, (ns, w)) in enumerate(weights)
            span, str = 0.0, 0.0
            for (j, n) in enumerate(ns)
                n >= nspan && (span += w[j])
                str += w[j] * largest[n + 1]
            end
            str /= N
            s_span[k] += span
            s2_span[k] += span^2
            s_str[k] += str
            s2_str[k] += str^2
        end
    end
    # mean and standard error of the mean (sample variance with n − 1)
    mean_err(s, s2) = (s ./ nsamples, sqrt.(max.(s2 .- s .^ 2 ./ nsamples, 0.0) ./ (nsamples - 1) ./ nsamples))
    return mean_err(s_span, s2_span), mean_err(s_str, s2_str)
end

function write_percolation(path)
    rng = Random.Xoshiro(20261029)
    k = findfirst(≈(0.5), PS)
    results = map(SIZES) do L
        t = @elapsed (span, span_e), (str, str_e) = simulate(rng, L, NSAMPLES[L])
        println("L = $L: $(NSAMPLES[L]) samples in $(round(t; digits=1)) s, ",
            "R(1/2) = $(span[k]) ± $(span_e[k]), P(1/2) = $(str[k]) ± $(str_e[k])")
        # sanity check: the spanning probability at p_c is close to 1/2 on the square
        0.4 <= span[k] <= 0.6 || error("spanning probability at p = 1/2 is $(span[k]) for L = $L")
        return (spanning=(span, span_e), strength=(str, str_e))
    end
    mkpath(dirname(path))
    h5open(path, "w") do h
        h["Ls"] = SIZES
        h["ps"] = PS
        for (L, r) in zip(SIZES, results), obs in (:spanning, :strength)
            y, e = r[obs]
            h["$obs/L$L/x"] = PS
            h["$obs/L$L/y"] = y
            h["$obs/L$L/e"] = e
            attrs(h["$obs/L$L"])["nsamples"] = NSAMPLES[L]
        end
    end
end

write_percolation(isempty(ARGS) ? joinpath(@__DIR__, "..", "testdata", "percolation", "bond_square.h5") : ARGS[1])
