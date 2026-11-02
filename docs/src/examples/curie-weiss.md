```@meta
CurrentModule = ScalingCollapse
```

# Mean-field model: exact data

The Curie–Weiss model couples all ``N`` Ising spins with each other, with the energy
``H = -M^2/(2N)`` of the total magnetization ``M``. It is the mean-field limit of the Ising
model, with ``T_c = 1`` and the mean-field exponents ``\beta = 1/2`` and ``\gamma = 1``. The
role of the system size is played by ``N``, and the scaling variable is
``(T - T_c) N^{1/2}``, so the exponent `nu` of the presets is ``\nu = 2``.

For this model, the partition function is a sum over the ``N + 1`` values of ``M``, so the
observables of a finite system can be computed exactly. This example computes them on the fly
and fits noise-free data.

## Exact finite-size data

```@example cw
using ScalingCollapse

# probabilities of the magnetization per spin m = M/N at temperature T
function distribution(N, T)
    logw = zeros(N + 1)            # log binomial(N, k), by recursion over k
    for k in 0:N-1
        logw[k+2] = logw[k+1] + log((N - k) / (k + 1))
    end
    Ms = [2k - N for k in 0:N]
    logw .+= Ms .^ 2 ./ (2N * T)   # Boltzmann weight exp(M²/(2NT))
    w = exp.(logw .- maximum(logw))
    return Ms ./ N, w ./ sum(w)
end

# Binder cumulant and susceptibility
function observables(N, T)
    m, w = distribution(N, T)
    m1, m2, m4 = sum(w .* abs.(m)), sum(w .* m .^ 2), sum(w .* m .^ 4)
    return 1 - m4 / (3m2^2), N * (m2 - m1^2) / T
end

Ns = [256, 1024, 4096, 16384]
Ts = collect(range(0.9, 1.1; length=41))
obs = [observables.(N, Ts) for N in Ns]
U = [first.(o) for o in obs]
chi = [last.(o) for o in obs]
nothing # hide
```

## Binder cumulant

For noise-free data, cubic splines interpolate more accurately than the default monotone
interpolation, so we use `MultipleSplines(interp=:cubic)` (see [Quality functions](@ref)):

```@example cw
sp_U = ScalingProblem(Ts, U, Ns;
    sf=ScalingFunction(:x; p_names=["T_c", "nu"]),
    p_space=[0.9:0.01:1.1, 1.0:0.1:3.0],
    dx=[-2.0, 2.0],
    quality=MultipleSplines(interp=:cubic),
)
```

## Susceptibility

The susceptibility ``\chi = N(\langle m^2 \rangle - \langle |m| \rangle^2)/T`` has a peak that
grows as ``N^{\gamma/\nu}``, the preset `:xny`:

```@example cw
sp_chi = ScalingProblem(Ts, chi, Ns;
    sf=ScalingFunction(:xny; p_names=["T_c", "nu", "gamma"]),
    p_space=[0.9:0.01:1.1, 1.0:0.1:3.0, 0.5:0.1:1.5],
    dx=[-1.0, 1.0],
    quality=MultipleSplines(interp=:cubic),
)
```

```@example cw
using CairoMakie

colors = cgrad(:Blues)[range(0.35, 1.0; length=length(Ns))]
fig = Figure(size=(780, 560))
for (row, sp, data, ylabel, ylabel_scaled) in (
    (1, sp_U, U, L"U", L"U"),
    (2, sp_chi, chi, L"\chi", L"\chi\, N^{-\gamma/\nu}"))
    ax1 = Axis(fig[row, 1]; xlabel=L"T", ylabel)
    ax2 = Axis(fig[row, 2]; xlabel=L"(T - T_c)\,N^{1/\nu}", ylabel=ylabel_scaled)
    xs_s, ys_s, _, _ = scaled_data(sp)
    for i in eachindex(Ns)
        scatterlines!(ax1, Ts, data[i]; color=colors[i], markersize=5, label="N = $(Ns[i])")
        inside = -4 .<= xs_s[i] .<= 4
        scatter!(ax2, xs_s[i][inside], ys_s[i][inside]; color=colors[i], markersize=7)
    end
    vlines!(ax2, sp.dx; color=:gray, linestyle=:dash)
    row == 1 && Legend(fig[1:2, 3], ax1)
end
fig
```

## Exact data are not free of errors

Both fits agree with ``T_c = 1``, ``\nu = 2`` and ``\gamma = 1`` to within two to three
errors, but not exactly, although the data have no noise. The errors here measure the corrections to
scaling, which decay as ``N^{-1/2}`` in this model: the scaling form is only the leading term
for large ``N``. As in the [percolation example](@ref "Corrections to scaling"), larger
systems or fewer small ones reduce this bias.
