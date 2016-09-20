
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

pr,eq,tau,pa=load("ModelResults.jld", "pr","eq","tau","pa");

taucvec = [0.34 0.33 0.32 0.31 0.30 0.29 0.28 0.27 0.26 0.25 0.24 0.23 0.22 0.21 0.20 0.19 0.18 0.17 0.16 0.15 0.14 0.13 0.12 0.11 0.10 0.09 0.08 0.07 0.06 0.05 0.04 0.03 0.02 0.01 0.00];
ref2 = Reform2Vector("reform2.jld", taucvec, pr, eq, tau, pa; bctol= 10.0^-3.0, update=0.75);
ref3 = Reform3Vector("reform3.jld", taucvec, pr, eq, tau, pa; bctol= 10.0^-3.0, update=0.75);

tauivec = [0.27 0.26 0.25 0.24 0.23 0.22 0.21 0.20 0.19 0.18 0.17 0.16 0.15 0.14 0.13 0.12 0.11 0.10 0.09 0.08 0.07 0.06 0.05 0.04 0.03 0.02 0.01 0.00 ];
re2taui = Reform2TauIVector("Reform2TauIVector.jld", tauivec, pr, eq, tau, pa; bctol= 10.0^-3.0, update=0.75);



prg,eqg,taug,pag=load("ModelResultsNoTaxG.jld", "pr","eq","tau","pa");

taucvec = [0.30 0.27 0.25 0.23 0.21 0.19 0.17 0.15 0.13 0.11 0.09 0.07 0.05 0.03 0.01 0.00];
ref4 = Reform4Vector("reform4.jld", taucvec, pr, eq, tau, pa; bctol= 5.0*10.0^-3.0, update=0.75);

taucvec =  [0.30 0.27 0.25 0.23 0.21 0.19 0.17 0.15 ];
ref1 = Reform1Vector("reform1.jld", taucvec, prg, eqg, taug, pag; bctol= 5.0*10.0^-3.0, update=0.75);
