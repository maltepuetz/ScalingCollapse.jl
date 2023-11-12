module Scaling

# Write your package code here.

using Dierckx
using Optim
using Statistics

include("Data.jl")
include("ScalingFunction.jl")
include("optimization.jl")
include("error_analysis.jl")
include("ScalingProblem.jl")

#=
# TODO add functions to extract 'ready to plot' data from ScalingProblem
# TODO add Measurements.jl support - and parameter export
=#
end
