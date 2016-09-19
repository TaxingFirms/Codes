
@everywhere using Grid:CoordInterpGrid, BCnan, BCnearest, InterpLinear
using JLD

include("markov_approx.jl")
include("mc_tools.jl")

@everywhere include("Main.jl")
@everywhere include("Firms.jl")
@everywhere include("FreeEntry.jl")
@everywhere include("Distribution.jl")
@everywhere include("Aggregation.jl")
@everywhere include("SolveSteadyState.jl")

@everywhere include("LafferCurve.jl")


tauivec = 0.02:0.02:0.4
@time laffer_taui_parallel(tauivec, "LafferTauI.jld","ModelResults.jld")

tauevec = 0.16:0.02:0.54
@time laffer_taue_parallel(tauevec, "LafferTauE.jld","ModelResults.jld")

taudvec = 0.02:0.02:0.4
@time laffer_taud_parallel(taudvec, "LafferTauD.jld","ModelResults.jld")

taucvec = 0.16:0.02:0.54
@time laffer_tauc_parallel(taucvec, "LafferTauC.jld","ModelResults.jld")
