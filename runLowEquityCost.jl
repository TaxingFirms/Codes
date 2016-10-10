
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

#Tesla
pa =init_parameters(H = 3.4, psi=0.5, ddelta = 0.0769, ttheta = 0.22956 , rhoz = 0.75, ssigmaz = 0.10608, llambda0 = 0.0456, llambda1 = 0.08897, ff = 1.3267, e=0.0260);

tau = init_taxes(ttaud = 0.15, ttauc = 0.35, ttaui = 0.0, ttaug = 0.0, ttaul = 0.28);
@time pr,eq= SolveSteadyState(tau,pa;wguess=0.5, VFItol=10.0^-5.0, displayit0=true, displayw = true);
moments=computeMomentsCutoff(eq.E,pr,eq,tau,pa,cutoffCapital=0.0;toPrint=true);
save("ResultsZeroTauI.jld","pr",pr,"eq",eq,"tau",tau,"pa",pa);

#pr,eq,tau,pa =load("ResultsZeroTauI.jld","pr","eq","tau","pa");

include("TaxReforms.jl")
include("Reforms.jl")


taucvec = [0.33 0.31 0.29 0.27 0.25 0.23 0.21 0.19 0.17 0.15 0.13 0.11 0.09 0.07 0.05 0.03 0.01 0.00];
ref5 = Reform5Vector("ref5LowEquityCost.jld", taucvec, pr, eq, tau, pa; bctol= 5*10.0^-3.0, update=0.0);
#ref5,pa=load("ref5CoarseZgrid.jld","ref","pa");
