# Generates the glow behind the logo on the home page,
# docs/src/assets/hero-glow-{light,dark}.png (used in src/.vitepress/theme/overrides.css).
# Run from the repository root: julia --project=docs docs/hero_glow.jl
#
# The glow is one Gaussian blob per dot of the logo. It is an image rather than a CSS
# gradient because a faint gradient spans only a few color levels on the dark background,
# and Safari draws CSS gradients without dithering, so every level shows as a ring. Here
# the alpha channel is dithered with triangular noise before it is rounded to 8 bits, with
# an amplitude of about one level of the composited color on the page background.
using CairoMakie, Random

const SIZE = 560  # pixels; the largest size of the element in CSS pixels

# Blobs as (center x, center y, radius) in units of the image width, with the color of
# the dot. The first one is drawn on top.
const BLOBS = [
    (0.64, 0.36, 0.30, :purple),
    (0.50, 0.50, 0.38, :green),
    (0.36, 0.64, 0.30, :red),
]

# Page background, and color and peak opacity of each blob, per theme.
const THEMES = Dict(
    "light" => Dict(:background => RGBf(1, 1, 1),
                    :purple => (RGBf(170, 121, 193) / 255, 0.29),
                    :green => (RGBf(56, 152, 38) / 255, 0.27),
                    :red => (RGBf(213, 99, 92) / 255, 0.29)),
    "dark" => Dict(:background => RGBf(27, 27, 31) / 255,
                   :purple => (RGBf(149, 88, 178) / 255, 0.08),
                   :green => (RGBf(56, 152, 38) / 255, 0.10),
                   :red => (RGBf(203, 60, 51) / 255, 0.08)),
)

# Gaussian exp(-4.5 t^2) in the distance t from the center in units of the radius,
# shifted and rescaled so that it reaches zero at t = 1.
profile(t) = t >= 1 ? 0.0 : (exp(-4.5t^2) - exp(-4.5)) / (1 - exp(-4.5))

function glow(theme; rng=Xoshiro(1))
    img = Matrix{Makie.RGBA{Makie.N0f8}}(undef, SIZE, SIZE)
    for j in 1:SIZE, i in 1:SIZE
        x, y = (j - 0.5) / SIZE, (i - 0.5) / SIZE
        # Premultiplied color and alpha, composited from the bottom blob up.
        r = g = b = a = 0.0
        for (cx, cy, radius, name) in reverse(BLOBS)
            color, peak = theme[name]
            α = peak * profile(hypot(x - cx, y - cy) / radius)
            r = α * color.r + (1 - α) * r
            g = α * color.g + (1 - α) * g
            b = α * color.b + (1 - α) * b
            a = α + (1 - α) * a
        end
        if a == 0
            img[i, j] = Makie.RGBA{Makie.N0f8}(0, 0, 0, 0)
            continue
        end
        # Triangular noise before rounding. One level of alpha changes the composited color
        # by `contrast` levels in its strongest channel; the noise spans one color level.
        bg = theme[:background]
        contrast = max(abs(r / a - bg.r), abs(g / a - bg.g), abs(b / a - bg.b))
        noise = (rand(rng) - rand(rng)) / max(contrast, 0.05)
        a8 = clamp(round(Int, 255a + noise), 0, 255)
        img[i, j] = a8 == 0 ? Makie.RGBA{Makie.N0f8}(0, 0, 0, 0) :
            Makie.RGBA{Makie.N0f8}(min(r / a, 1), min(g / a, 1), min(b / a, 1), a8 / 255)
    end
    return img
end

for (name, theme) in THEMES
    path = joinpath(@__DIR__, "src", "assets", "hero-glow-$name.png")
    save(path, glow(theme))
    println(path, ": ", filesize(path) ÷ 1024, " KiB")
end
