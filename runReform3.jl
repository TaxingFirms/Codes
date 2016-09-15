
@everywhere using Grid:CoordInterpGrid, BCnan, BCnearest, InterpLinear
using JLD
using StatsFuns

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

pr,eq,tau,pa=load("ModelResults.jld", "pr","eq","tau","pa");

include("TaxReforms.jl")
include("Reforms.jl")


taucvec1 = [0.33 0.30 0.27 0.24 0.21 0.18 0.17 0.16 0.155 ];
reform1 = Reform1Vector("reform1.jld", taucvec1, pr, eq, tau, pa)

taucvec = [0.33 0.30 0.27 0.24 0.21 0.18 0.15 0.12 0.09 0.06 0.03 0.00];
reform3 = Reform3Vector("reform3.jld", taucvec, pr, eq, tau, pa)
