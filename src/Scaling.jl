module Scaling

# Write your package code here.

using Dierckx
using Optim
using Statistics



export ScalingProblem, ScalingFunction

include("data.jl")
include("scalingfunction.jl")
include("quality_houdayer.jl")
include("quality_spline.jl")
include("optimization.jl")
include("error_analysis.jl")
include("scalingproblem.jl")

#=
# TODO add tests
# TODO add functions to extract 'ready to plot' data from ScalingProblem
# TODO add Measurements.jl support - and parameter export (maybe)
# TODO package extension for plotting (maybe)
=#
end
