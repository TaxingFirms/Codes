
@everywhere using Grid:CoordInterpGrid, BCnan, BCnearest, InterpLinear
using NPZ
using Sobol

#using DataFrames

include("markov_approx.jl")
include("mc_tools.jl")

@everywhere include("Main.jl")
@everywhere include("Firms.jl")
@everywhere include("FreeEntry.jl")
@everywhere include("Distribution.jl")
@everywhere include("Aggregation.jl")
@everywhere include("SolveSteadyState.jl")
@everywhere include("calibrate.jl")
@everywhere include("Transitions.jl")



include("TaxMaximization.jl")
include("Temp.jl")

# 1. Create tax sequence from initial point
  #First component is tauc, second is taui
initialpoint = [0.0,0.28]
create_taxspace(initialpoint, 1500; jldfile=false)

npzwrite("taxspace.npy",taxspace)

#govexp = 0.046970754282188165;
#taul = 0.28 ; wguess0=0.532;
#maximize_welfare_tesla(taxspace::Array{Float64,2},govexp::Float64, taul::Float64, wguess0::Float64, pa::Param)
