# Exact finite-N Curie–Weiss model (H = −M²/(2N), T_c = 1). Writes the Binder ratio
# U = 1 − ⟨M⁴⟩/(3⟨M²⟩²) and m_rms = sqrt(⟨M²⟩)/N for N = 64, 128, 256, 512, 1024 on 61
# temperatures T ∈ [0.8, 1.2] to test/testdata/exact/curie_weiss.h5. Not run in CI.
#
# Z_N(T) = Σ_k binomial(N, k) exp(M²/(2NT)) with M = 2k − N, summed in log space.
#
# Layout: `Ls` (the N values), `binder/L{N}/{x, y}` and `m_rms/L{N}/{x, y}`, with x = T.
#
# Usage (in an environment with HDF5):
#     julia curie_weiss.jl [out.h5]

using HDF5

# log binomial(N, k) for k = 0:N by the recursion C(N, k+1) = C(N, k) (N − k)/(k + 1)
function logbinomials(N)
    lb = zeros(N + 1)
    for k in 0:N-1
        lb[k+2] = lb[k+1] + log((N - k) / (k + 1))
    end
    return lb
end

# ⟨M²⟩ and ⟨M⁴⟩ at temperature T
function moments(N, T)
    Ms = [2k - N for k in 0:N]
    logw = logbinomials(N) .+ Ms .^ 2 ./ (2N * T)
    w = exp.(logw .- maximum(logw))
    Z = sum(w)
    return sum(w .* Ms .^ 2) / Z, sum(w .* Ms .^ 4) / Z
end

binder(N, T) = ((m2, m4) = moments(N, T); 1 - m4 / (3m2^2))
m_rms(N, T) = sqrt(moments(N, T)[1]) / N

# sanity checks: infinite temperature (independent spins) and the N → ∞ limits around T_c
@assert isapprox(moments(64, 1e12)[1], 64; rtol=1e-9)
@assert isapprox(binder(64, 1e12), 1 - (3 * 64^2 - 2 * 64) / (3 * 64^2); rtol=1e-6)
@assert binder(4096, 0.5) > 0.66 && binder(4096, 2.0) < 0.01
# universal U(T_c) = 1 − Γ(1/4)⁴/(24π²) for N → ∞ (Γ(1/4) = 3.625609908221908…)
@assert abs(binder(2^16, 1.0) - (1 - 3.625609908221908^4 / (24π^2))) < 1e-2

out = isempty(ARGS) ? joinpath(@__DIR__, "..", "testdata", "exact", "curie_weiss.h5") : ARGS[1]
Ns = [64, 128, 256, 512, 1024]
Ts = collect(range(0.8, 1.2; length=61))
h5open(out, "w") do h
    h["Ls"] = Ns
    for N in Ns
        h["binder/L$N/x"] = Ts
        h["binder/L$N/y"] = [binder(N, T) for T in Ts]
        h["m_rms/L$N/x"] = Ts
        h["m_rms/L$N/y"] = [m_rms(N, T) for T in Ts]
    end
end
