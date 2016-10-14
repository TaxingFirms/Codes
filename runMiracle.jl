
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

pa =init_parameters(H = 3.4, psi=0.5, ddelta = 0.0768, ttheta = 0.2295 , rhoz = 0.75, ssigmaz = 0.0993, llambda0 = 0.0255, llambda1 = 0.23958, ff = 1.44459, e=0.041674);

tau = init_taxes(ttaud =0.15, ttauc= 0.35, ttaui= 0.0, ttaug= 0.15, ttaul=0.28, ttauexit=0.0);
@time pr,eq= SolveSteadyState(tau,pa;wguess=0.5, VFItol=10.0^-5.0, displayit0=true, displayw = true);
moments=computeMomentsCutoff(eq.E,pr,eq,tau,pa,cutoffCapital=0.0;toPrint=true);
save("MiracleResults.jld","pr",pr,"eq",eq,"tau",tau,"pa",pa);

include("TaxReforms.jl")
include("Reforms.jl")



#pr,eq,tau,pa = load("ModelResults.jld","pr","eq","tau","pa");
tau.i=0.0;
taucvec = [0.33 0.31 0.29 0.27 0.25 0.23 0.21 0.19 0.17 0.15 0.13 0.11 0.09 0.07 0.05 0.03 0.01 0.00];

ref1 = Reform1Vector("reform1ZeroTauI.jld", taucvec, pr, eq, tau, pa; bctol= 5*10.0^-4.0, update=0.75);
ref2 = Reform2Vector("reform2ZeroTauI.jld", taucvec, pr, eq, tau, pa; bctol= 5*10.0^-4.0, update=0.75);
ref3 = Reform3Vector("reform3ZeroTauI.jld", taucvec, pr, eq, tau, pa; bctol= 5*10.0^-4.0, update=0.75);

tau.i=0.28;
ref2 = Reform2Vector("reform2.jld", taucvec, pr, eq, tau, pa; bctol= 5*10.0^-4.0, update=0.75);
ref3 = Reform3Vector("reform3.jld", taucvec, pr, eq, tau, pa; bctol= 5*10.0^-4.0, update=0.75);
