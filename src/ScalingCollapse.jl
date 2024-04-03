module ScalingCollapse

using Dierckx
using Measurements
using Optim
using Statistics

export ScalingProblem
export ScalingFunction
export SingleMasterCurve
export Spline
export Houdayer
export SingleSpline
export residuals
export scaled_data

include("data.jl")
include("scalingfunction.jl")
include("quality.jl")
include("optimization.jl")
include("error_analysis.jl")
include("scalingproblem.jl")
include("residuals.jl")

end
