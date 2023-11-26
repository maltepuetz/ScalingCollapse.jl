module Scaling

# Write your package code here.

using Dierckx
using Measurements
using Optim
using Statistics



export ScalingProblem, ScalingFunction, Spline, Houdayer

include("data.jl")
include("scalingfunction.jl")
include("quality.jl")
#include("quality_houdayer.jl")
#include("quality_spline.jl")
#include("quality_hybrid.jl")
include("optimization.jl")
include("error_analysis.jl")
include("scalingproblem.jl")
include("residual_landscape.jl")

end
