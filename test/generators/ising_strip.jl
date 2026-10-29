# Exact correlation length ξ_L of the 2D Ising model (J = 1) on a periodic strip of width L,
# from Kaufman's transfer-matrix spectrum. Writes y = ξ_L/L for L = 8, 16, 32, 64 on 61
# temperatures T ∈ [2.1, 2.45] to test/testdata/exact/ising_strip.h5. Not run in CI.
#
# With K = 1/T, tanh K* = e^{-2K} and
#     cosh γ_q = cosh 2K* cosh 2K − sinh 2K* sinh 2K cos(πq/L),
# the two largest eigenvalues (even and odd under global spin flip) give
#     ξ_L⁻¹ = log(λ₀/λ₁) = ½ Σ_{r=0}^{L−1} (γ_{2r+1} − γ_{2r}),
# where γ₀ = 2(K − K*) carries the sign of K − K*. This convention is asserted below against
# brute-force diagonalisation of the 2^L × 2^L transfer matrix.
#
# Usage (in an environment with HDF5):
#     julia ising_strip.jl [out.h5]

using HDF5, LinearAlgebra

Kstar(K) = atanh(exp(-2K))

# γ_q for q ≥ 1 in the cancellation-free form sinh²(γ/2) = sinh²(K*−K) + sinh 2K* sinh 2K sin²(πq/2L)
function γ(K, L, q)
    Ks = Kstar(K)
    q == 0 && return 2(K - Ks)
    return 2asinh(sqrt(sinh(Ks - K)^2 + sinh(2Ks) * sinh(2K) * sin(π * q / (2L))^2))
end

xi_kaufman(K, L) = 1 / (sum(γ(K, L, 2r + 1) - γ(K, L, 2r) for r in 0:L-1) / 2)

# 1/log(λ₀/λ₁) with λ₀ (λ₁) the largest eigenvalue of the symmetric transfer matrix
# T(s, s') = exp(K E(s)/2 + K s·s' + K E(s')/2) in the flip-even (flip-odd) sector
function xi_bruteforce(K, L)
    spins(i) = [isodd(i >> j) ? -1 : 1 for j in 0:L-1]
    E(s) = sum(s[i] * s[mod1(i + 1, L)] for i in 1:L)
    T(s, t) = exp(K * (E(s) + E(t)) / 2 + K * dot(s, t))
    S = [spins(i) for i in 0:2^(L-1)-1]   # representatives with s_L = +1
    even = [T(s, t) + T(s, -t) for s in S, t in S]
    odd = [T(s, t) - T(s, -t) for s in S, t in S]
    return 1 / log(maximum(eigvals(Symmetric(even))) / maximum(eigvals(Symmetric(odd))))
end

const Tc = 2 / log(1 + √2)

for L in (4, 6, 8), T in (2.0, 2.269, 2.6)
    exact, brute = xi_kaufman(1 / T, L), xi_bruteforce(1 / T, L)
    @assert isapprox(exact, brute; rtol=1e-10) "Kaufman ξ = $exact ≠ brute force $brute at L = $L, T = $T"
end
let r = xi_kaufman(1 / Tc, 256) / 256
    @assert abs(r - 4 / π) < 1e-3 "ξ_L/L at T_c and L = 256 is $r, expected → 4/π"
end

out = isempty(ARGS) ? joinpath(@__DIR__, "..", "testdata", "exact", "ising_strip.h5") : ARGS[1]
Ls = [8, 16, 32, 64]
Ts = collect(range(2.1, 2.45; length=61))
h5open(out, "w") do h
    h["Ls"] = Ls
    for L in Ls
        h["L$L/x"] = Ts
        h["L$L/y"] = [xi_kaufman(1 / T, L) / L for T in Ts]
    end
end
