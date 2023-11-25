struct ScalingFunction

    f::Function
    p_names::Vector{String}
    x_scale::String
    y_scale::String
    fixed_ps::Vector{Float64}

    function ScalingFunction(preset::Symbol; kwargs...)

        #=
        Presets:
            :x    -> x scaling with positive exponent
            :xy   -> x and y scaling with positive exponents
            :xny  -> x and y scaling with (pos x), (neg y) exponent
        =#

        if preset == :x
            f = get(kwargs, :f, _power_scaling)
            p_names = get(kwargs, :p_names, ["p1", "p2"])
            x_scale = get(
                kwargs,
                :x_scale,
                "(x - $(p_names[1]))/$(p_names[1]) * L^(1/$(p_names[2]))"
            )
            y_scale = get(kwargs, :y_scale, "y")

            # check whether a parameter was provided as kwarg
            fixed_ps = [Inf for _ in 1:length(p_names)]
            for (i, p) in enumerate(p_names)
                if haskey(kwargs, Symbol(p))
                    fixed_ps[i] = kwargs[Symbol(p)]
                end
            end

            return new(
                (d, args...) -> f(
                    d,
                    (_i = 0; fixed_ps[1] == Inf) ? (_i += 1; args[_i]) : fixed_ps[1],
                    fixed_ps[2] == Inf ? (_i += 1; args[_i]) : fixed_ps[2]
                ),
                p_names,
                x_scale,
                y_scale,
                fixed_ps
            )

        elseif preset == :xy
            f = get(kwargs, :f, _power_scaling)
            p_names = get(kwargs, :p_names, ["p1", "p2", "p3"])
            x_scale = get(
                kwargs,
                :x_scale,
                "(x - $(p_names[1]))/$(p_names[1]) * L^(1/$(p_names[2]))")
            y_scale = get(
                kwargs,
                :y_scale,
                "y * L^($(p_names[3])/$(p_names[2]))"
            )

            # check whether a parameter was provided as kwarg
            fixed_ps = [Inf for _ in 1:length(p_names)]
            for (i, p) in enumerate(p_names)
                if haskey(kwargs, Symbol(p))
                    fixed_ps[i] = kwargs[Symbol(p)]
                end
            end

            return new(
                (d, args...) -> f(
                    d,
                    (_i = 0; fixed_ps[1] == Inf) ? (_i += 1; args[_i]) : fixed_ps[1],
                    fixed_ps[2] == Inf ? (_i += 1; args[_i]) : fixed_ps[2],
                    fixed_ps[3] == Inf ? (_i += 1; args[_i]) : fixed_ps[3]
                ),
                p_names,
                x_scale,
                y_scale,
                fixed_ps
            )

            #return new(f, p_names, x_scale, y_scale, zeros(3)) # TODO fixed_ps

        elseif preset == :xny
            # put negative exponent on y
            f = get(
                kwargs,
                :f,
                (d::Data, p1, p2, p3) -> _power_scaling(d::Data, p1, p2, -p3)
            )
            p_names = get(kwargs, :p_names, ["p1", "p2", "p3"])
            x_scale = get(
                kwargs,
                :x_scale,
                "(x - $(p_names[1]))/$(p_names[1]) * L^(1/$(p_names[2]))"
            )
            y_scale = get(
                kwargs,
                :y_scale,
                "y * L^(-$(p_names[3])/$(p_names[2]))"
            )

            # check whether a parameter was provided as kwarg
            fixed_ps = [Inf for _ in 1:length(p_names)]
            for (i, p) in enumerate(p_names)
                if haskey(kwargs, Symbol(p))
                    fixed_ps[i] = kwargs[Symbol(p)]
                end
            end

            return new(
                (d, args...) -> f(
                    d,
                    (_i = 0; fixed_ps[1] == Inf) ? (_i += 1; args[_i]) : fixed_ps[1],
                    fixed_ps[2] == Inf ? (_i += 1; args[_i]) : fixed_ps[2],
                    fixed_ps[3] == Inf ? (_i += 1; args[_i]) : fixed_ps[3]
                ),
                p_names,
                x_scale,
                y_scale,
                fixed_ps
            )
        else
            error("Unknown preset: $preset")
        end
    end

    function ScalingFunction(f::Function; kwargs...)
        p_names = if haskey(kwargs, :p_names)
            kwargs[:p_names]
        else
            if haskey(kwargs, :N_parameters)
                ["p$i" for i in 1:kwargs[:N_parameters]]
            else
                error("Please specify either :p_names or :N_parameters")
            end
        end
        x_scale = get(kwargs, :x_scale, "tmp")
        y_scale = get(kwargs, :y_scale, "tmp")
        return new(f, p_names, x_scale, y_scale, zeros(3)) # TODO fixed_ps
    end

    function ScalingFunction(p_names::Vector{String}; kwargs...)
        length(p_names) == 2 && return ScalingFunction(:x, p_names=p_names; kwargs...)
        length(p_names) == 3 && return ScalingFunction(:xy, p_names=p_names; kwargs...)
    end

    function ScalingFunction(; kwargs...)
        haskey(kwargs, :p_names) && return ScalingFunction(kwargs[:p_names]; kwargs...)
        return ScalingFunction(:xy; kwargs...)
    end
end

function n_parameters(sf::ScalingFunction)
    n = length(sf.fixed_ps)
    for p in sf.fixed_ps
        if p != Inf
            n -= 1
        end
    end
    return n
end

function scaled_p_names(sf::ScalingFunction)
    res = sf.p_names |> copy
    for i in length(sf.fixed_ps):-1:1
        if sf.fixed_ps[i] != Inf
            deleteat!(res, i)
        end
    end
    return res
end

function fixed_p_names(sf::ScalingFunction)
    res = sf.p_names |> copy
    for i in length(sf.fixed_ps):-1:1
        if sf.fixed_ps[i] == Inf
            deleteat!(res, i)
        end
    end
    return res
end

# standart scaling function for 2 parameters
function _power_scaling(d::Data, p1, p2)

    xs = zeros(length(d.xs))
    ys = zeros(length(d.ys))
    es = zeros(length(d.es))

    for (i, x) in enumerate(d.xs)
        xs[i] = (x - p1) / p1 * d.L^(1 / p2)
    end
    for (i, y) in enumerate(d.ys)
        ys[i] = y
    end
    for (i, e) in enumerate(d.es)
        es[i] = e
    end

    return Data(d.L, xs, ys, es)
end

# standart scaling function for 3 parameters
function _power_scaling(d::Data, p1, p2, p3)

    xs = zeros(length(d.xs))
    ys = zeros(length(d.ys))
    es = zeros(length(d.es))

    for (i, x) in enumerate(d.xs)
        xs[i] = (x - p1) / p1 * d.L^(1 / p2)
    end
    for (i, y) in enumerate(d.ys)
        ys[i] = y * d.L^(p3 / p2)
    end
    for (i, e) in enumerate(d.es)
        es[i] = e * d.L^(p3 / p2)
    end

    return Data(d.L, xs, ys, es)
end
