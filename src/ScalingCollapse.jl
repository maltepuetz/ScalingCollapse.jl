module ScalingCollapse

using Dierckx
using KissSmoothing
using Measurements
using Optim
using Statistics

export ScalingProblem
export ScalingFunction
export residuals
export scaled_data

# export quality functions
export MultipleSplines, SingleSpline, Linear, Houdayer

include("data.jl")
include("scalingfunction.jl")
include("quality.jl")
include("optimization.jl")
include("error_analysis.jl")
include("scalingproblem.jl")
include("residuals.jl")

end
