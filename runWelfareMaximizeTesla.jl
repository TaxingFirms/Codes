
@everywhere using Grid:CoordInterpGrid, BCnan, BCnearest, InterpLinear
using JLD
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

ref,pa = load("reform2.jld","ref","pa");

include("TaxMaximization.jl")
include("Temp.jl")



taxspace=load("taxspace.jld","taxspace")


govexp = ref[end].eq.a.G;
taul = ref[end].tau.l ; wguess0=ref[end].eq.w;
maximize_welfare_tesla(taxspace::Array{Float64,2},govexp::Float64, taul::Float64, wguess0::Float64, pa::Param)
