# Generates the README figure, docs/src/assets/readme-collapse-{light,dark}.png, in the style
# of the docs figures (see figure_theme.jl).
# Run from the repository root: julia --project=docs docs/readme_figure.jl
using ScalingCollapse, HDF5, Printf
include("figure_theme.jl")

file = joinpath(pkgdir(ScalingCollapse), "test", "testdata", "IsingMonteCarlo.h5")
Ts, Ls, binder = h5open(h -> (read(h, "Ts"), read(h, "Ls"), read(h, "binder")), file)
sp = ScalingProblem(Ts, binder, Ls;
    sf=ScalingFunction(:x; p_names=["T_c", "nu"]),
    p_space=[2.0:0.02:2.6, 0.5:0.05:1.5],
    dx=[-2.0, 2.0],
)

colors = cgrad(:Blues)[range(0.35, 1.0; length=length(Ls))]
fig = Figure(size=(780, 300); backgroundcolor=:transparent)
ax1 = Axis(fig[1, 1]; xlabel=L"T", ylabel=L"U", title="Binder cumulant, 2D Ising model")
ax2 = Axis(fig[1, 2]; xlabel=L"(T - T_c)\,L^{1/\nu}", ylabel=L"U",
    title=Makie.LaTeXString("\$T_c = " *
        @sprintf("%.3f \\pm %.3f", sp.optimal_ps[1], sp.optimal_ps_error[1]) *
        ",\\quad \\nu = " * @sprintf("%.2f \\pm %.2f", sp.optimal_ps[2], sp.optimal_ps_error[2]) *
        "\$"))
xs, ys, _, _ = scaled_data(sp)
for i in eachindex(Ls)
    scatterlines!(ax1, Ts, binder[:, i]; color=colors[i], markersize=5, label="L = $(Ls[i])")
    inside = -6 .<= xs[i] .<= 6
    scatter!(ax2, xs[i][inside], ys[i][inside]; color=colors[i], markersize=7)
end
Legend(fig[1, 3], ax1)

assets = joinpath(@__DIR__, "src", "assets")
write(joinpath(assets, "readme-collapse-light.png"), _png(fig))
write(joinpath(assets, "readme-collapse-dark.png"), _dark_png(fig))
