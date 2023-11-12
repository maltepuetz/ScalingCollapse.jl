var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = Scaling","category":"page"},{"location":"#Scaling","page":"Home","title":"Scaling","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for Scaling.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [Scaling]","category":"page"},{"location":"#Scaling.Data","page":"Home","title":"Scaling.Data","text":"Data(L::Int, xs::Vector{Float64}, ys::Vector{Float64}, es::Vector{Float64})\n\nA Data object stores the data for a single system size.\n\nFields\n\nL::Int: system size\nxs::Vector{Float64}: x values\nys::Vector{Float64}: y values\nes::Vector{Float64}: y error values\n\n\n\n\n\n","category":"type"},{"location":"#Scaling.ScalingProblem","page":"Home","title":"Scaling.ScalingProblem","text":"ScalingProblem(args...; kwargs...)\n\nCreate a scaling problem which solves on initialization.\n\nArguments\n\nxs: x values of the data\nys: y values of the data\nes: y error values of the data (optional)\nLs: system sizes of the data\n\nFor more information on the arguments, see methods(Scaling.unzip_data).\n\nKeyword Arguments\n\nsf::ScalingFunction=ScalingFunction(; kwargs...): scaling function\np_space::Vector{StepRangeLen}=[0.1:0.1:3.0 for _ in sf.p_names]: search parameter space\ndx::Vector{Float64}=[-Inf, Inf]: optimization interval\nverbose::Bool=false: print information during optimization\nstarting_ps::Vector{Float64}: If starting_ps are given, there will be no initial   parameter space scan. Instead, the optimization will start at the given points. This is   much faster, but might not find the global minimum.\nerror::Bool=false: If error=true, the error analysis will be performed to give   estimates of the uncertainties of the optimal parameters.\n\nFields\n\ndata::Vector{Data}: data to be scaled\nsf::ScalingFunction: scaling function\np_space::Vector{StepRangeLen}: search parameter space\ndx::Vector{Float64}: optimization interval\noptimal_ps::Vector{Float64}: optimal parameters\noptimal_ps_error::Vector{Float64}: uncertainties of the optimal parameters\nminimum::Float64: minimum of the objective function\n\nExamples\n\nusing Scaling\n\nrescale the x axis only\n\n(e.g. for ys == binder cumulant of Ising model)\n\nScalingProblem(xs, ys, Ls;\n    sc=ScalingFunction(:x; p_names=[\"T_c\", \"nu\"]),\n    p_space=[0.1:0.1:3.0, 0.1:0.1:3.0],\n    dx=[-1.0, 1.0],\n)\n\n# or the same but with errors\nScalingProblem(xs, ys, es, Ls;\n    sc=ScalingFunction(:x; p_names=[\"T_c\", \"nu\"]),\n    p_space=[0.1:0.1:3.0, 0.1:0.1:3.0],\n    dx=[-1.0, 1.0],\n)\n\nrescale the x and y axis\n\n(e.g. for ys == magnetization of Ising model)\n\nScalingProblem(xs, ys, Ls;\n    sc=ScalingFunction(:xy; p_names=[\"T_c\", \"nu\", \"beta\"]),\n    dx=[-1.0, 1.0],\n)\n\nrescale the x and y axis, but y axis with negative exponent\n\n(e.g. for ys == susceptibility of Ising model)\n\nScalingProblem(xs, ys, Ls;\n    sc=ScalingFunction(:xny; p_names=[\"T_c\", \"nu\", \"gamma\"]),\n    dx=[-1.0, 1.0],\n)\n\n\n\n\n\n","category":"type"},{"location":"#Scaling.squared_sum_residuals-Tuple{Vector{Scaling.Data}, Any, Any}","page":"Home","title":"Scaling.squared_sum_residuals","text":"squared_sum_residuals(data::Vector{Data}, scaling_function, parameters; kwargs...)\n\nCalculate the squared sum of residuals for a given scaling function and parameters.\n\nArguments\n\ndata::Vector{Data}: data to be scaled\nscaling_function::Function: scaling function\nparameters: parameters for the scaling function\n\n\n\n\n\n","category":"method"}]
}
