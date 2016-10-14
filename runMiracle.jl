
@everywhere using Grid:CoordInterpGrid, BCnan, BCnearest, InterpLinear
using JLD
using DataFrames

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


include("TaxReforms.jl")
include("Reforms.jl")

pr,eq,tau,pa = load("ModelResults.jld","pr","eq","tau","pa");
tau.i=0.0;
taucvec = [0.33 0.31 0.29 0.27 0.25 0.23 0.21 0.19 0.17 0.15 0.13 0.11 0.09 0.07 0.05 0.03 0.01 0.00];

ref1 = Reform1Vector("reform1ZeroTauI.jld", taucvec, pr, eq, tau, pa; bctol= 5*10.0^-4.0, update=0.75);
ref2 = Reform2Vector("reform2ZeroTauI.jld", taucvec, pr, eq, tau, pa; bctol= 5*10.0^-4.0, update=0.75);
ref3 = Reform3Vector("reform3ZeroTauI.jld", taucvec, pr, eq, tau, pa; bctol= 5*10.0^-4.0, update=0.75);

tau.i=0.28;
ref2 = Reform2Vector("reform2.jld", taucvec, pr, eq, tau, pa; bctol= 5*10.0^-4.0, update=0.75);
ref3 = Reform3Vector("reform3.jld", taucvec, pr, eq, tau, pa; bctol= 5*10.0^-4.0, update=0.75);
