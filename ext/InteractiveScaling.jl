module InteractiveScaling

using Scaling
using GLMakie


@info "Loading InteractiveScaling.jl"

import GLMakie: plot
export plot

# this function returns the maximum and the minimum of the y data for x in the interval [x1, x2]
function get_limits(sp::ScalingProblem, sx, sy)

    xlims = (-2 * maximum(sp.dx .|> abs), 2 * maximum(sp.dx .|> abs))
    y_min = Inf
    y_max = -Inf
    for (l, d) in enumerate(sp.data)
        for (x, y) in zip(sx[l], sy[l])
            if x >= xlims[1] && x <= xlims[2]
                y_min = min(y_min, y)
                y_max = max(y_max, y)
            end
        end
    end
    ylims = (y_min - 0.1 * (y_max - y_min), y_max + 0.1 * (y_max - y_min))
    return xlims, ylims
end




function plot(sp::ScalingProblem;
    size=(800, 800),  # size of the figure
    errorbars=true,  # plot errorbars
    splines=false,  # plot splines
    dims=[2, 3],  # defined the dimensions to plot
    logcolorbar=true,
    n_levels=10,
    kwargs...
)

    # get the quality landscape and set up an observable for the residuals
    p_space, res = residuals(sp; dims=dims, kwargs...)
    reso = Observable(res)

    # set ps to optimal_ps and define px and py observable
    ps = copy(sp.optimal_ps)
    pxo = Observable(ps[dims[1]])
    pyo = Observable(ps[dims[2]])

    # set observable for dx
    dxo = Observable(sp.dx)


    fig = Figure(size=size)

    ##### raw data ############################
    ax_data = Axis(fig[1, 1:3])
    for (l, d) in enumerate(sp.data)
        if errorbars
            errorbars!(ax_data, d.xs, d.ys, d.es; whiskerwidth=7.0, color=:black)
        end
        scatter!(ax_data, d.xs, d.ys; label="L = $(d.L)")
    end
    axislegend(ax_data)



    ##### scaled data #########################
    if splines
        sx, sy, se, sL, xvals, yvals, evals = scaled_data(sp, ps; splines=splines)
    else
        sx, sy, se, sL = scaled_data(sp, ps; splines=splines)
    end

    ax_scaled = Axis(fig[2, 1:3],
        limits=get_limits(sp, sx, sy),
    )

    Makie.deactivate_interaction!(ax_scaled, :rectanglezoom)
    sline = select_line(ax_scaled.scene)

    on(sline) do s
        start = s[1]
        stop = s[2]

        absdx = abs(dxo[][2] - dxo[][1])
        index = 0  # index which border to move
        if abs(start[1] - dxo[][1]) < absdx / 4
            index = 1
        elseif abs(start[1] - dxo[][2]) < absdx / 4
            index = 2
        end

        if index != 0
            if index == 1
                dxo[] = [stop[1], dxo[][2]]
            elseif index == 2
                dxo[] = [dxo[][1], stop[1]]
            end
            sp.dx = dxo[]
            xlims, ylims = get_limits(sp, sx, sy)
            ax_scaled.limits = xlims, ylims

            if splines
                sx, sy, se, sL, xvals, yvals, evals = scaled_data(sp, ps; splines=splines)
            else
                sx, sy, se, sL = scaled_data(sp, ps; splines=splines)
            end
            for (l, d) in enumerate(sp.data)
                xo[l][] = sx[l]
                yo[l][] = sy[l]
                eo[l][] = se[l]
                if splines
                    xspo[l][] = xvals[l]
                    yspo[l][] = yvals[l]
                end
            end
            @async _, reso[] = residuals(sp; dims=dims, kwargs...)
        end
    end






    slider_x = Slider(fig[4, 2],
        # label=Scaling.scaled_p_names(sp.sf)[dims[1]],
        range=range(p_space[1][1], stop=p_space[1][end], length=200),
        startvalue=sp.optimal_ps[dims[1]],
    )

    slider_y = Slider(fig[3, 1],
        # label=Scaling.scaled_p_names(sp.sf)[dims[2]],
        range=range(p_space[2][1], stop=p_space[2][end], length=200),
        startvalue=sp.optimal_ps[dims[2]],
        horizontal=false,
    )

    xo = [
        Observable(sx[l]) for l in eachindex(sp.data)
    ]
    yo = [
        Observable(sy[l]) for l in eachindex(sp.data)
    ]
    eo = [
        Observable(se[l]) for l in eachindex(sp.data)
    ]
    if splines
        xspo = [
            Observable(xvals[l]) for l in eachindex(sp.data)
        ]
        yspo = [
            Observable(yvals[l]) for l in eachindex(sp.data)
        ]
    end

    band!(
        ax_scaled,
        dxo,
        -1000 .* ones(length(sp.dx)),
        1000 .* ones(length(sp.dx));
        color=(:black, 0.1),
    )

    for (l, d) in enumerate(sp.data)
        if errorbars
            errorbars!(ax_scaled, xo[l], yo[l], eo[l]; whiskerwidth=7.0, color=:black)
        end
        if splines
            lines!(ax_scaled, xspo[l], yspo[l])
        end
        scatter!(ax_scaled, xo[l], yo[l]; label="L = $(d.L)")
    end

    on(slider_x.value) do px
        pxo[] = px
        ps[dims[1]] = px
        if splines
            sx, sy, se, sL, xvals, yvals, evals = scaled_data(sp, ps; splines=splines)
        else
            sx, sy, se, sL = scaled_data(sp, ps; splines=splines)
        end
        for (l, d) in enumerate(sp.data)
            xo[l][] = sx[l]
            yo[l][] = sy[l]
            eo[l][] = se[l]
            if splines
                xspo[l][] = xvals[l]
                yspo[l][] = yvals[l]
            end
        end
    end

    on(slider_y.value) do py
        pyo[] = py
        ps[dims[2]] = py
        if splines
            sx, sy, se, sL, xvals, yvals, evals = scaled_data(sp, ps; splines=splines)
        else
            sx, sy, se, sL = scaled_data(sp, ps; splines=splines)
        end
        for (l, d) in enumerate(sp.data)
            xo[l][] = sx[l]
            yo[l][] = sy[l]
            eo[l][] = se[l]
            if splines
                xspo[l][] = xvals[l]
                yspo[l][] = yvals[l]
            end
        end
    end




    ##### residuals ############################
    ax_res = Axis(fig[3, 2];
        xlabel=Scaling.scaled_p_names(sp.sf)[dims[1]],
        ylabel=Scaling.scaled_p_names(sp.sf)[dims[2]],
    )

    Makie.deactivate_interaction!(ax_res, :rectanglezoom)




    levels = if logcolorbar
        [exp(x) for x in range(log(minimum(res)), stop=log(maximum(res)), length=n_levels)]
    else
        range(minimum(res), stop=maximum(res), length=n_levels)
    end
    co = contourf!(ax_res, p_space[1], p_space[2], reso;
        levels=levels,
    )
    scatter!(ax_res, pxo, pyo; color=:red)

    Colorbar(fig[3, 3], co, label="Residuals")

    spoint = select_point(ax_res.scene; marker=:circle)
    on(spoint) do s
        set_close_to!(slider_x, s[1])
        set_close_to!(slider_y, s[2])
    end
    return fig
end




end
