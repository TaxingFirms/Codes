# All you need to run the code
# Modify ~/.juliarc.jl and add the following line:
# push!(LOAD_PATH, "/Path/To/My/Module/")

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
pa =init_parameters(H = 1.094, ddelta=0.07557, ttheta = 0.2290 , rhoz =0.7451, ssigmaz = 0.1067, llambda0 = 0.02605, llambda1 = 0.2467, ff = 1.3856, e=0.01820);

#pa =init_parameters(H = 1.094, ddelta=0.0768, ttheta = 0.14067 , rhoz =0.75, ssigmaz = 0.086, llambda0 = 0.01074, llambda1 = 0.1581, ff = 1.3267, e=0.073289);
tau = init_taxes(ttaud =0.15, ttauc= 0.35, ttaui= 0.28, ttaug= 0.15, ttaul=0.28);
@time pr,eq= SolveSteadyState(tau,pa;wguess=0.5, VFItol=10.0^-3.0, displayit0=false, displayw = false);
moments=computeMomentsCutoff(eq.E,pr,eq,tau,pa,cutoffCapital=0.0;toPrint=false);
pcterror_params(pr,eq,tau,pa)
save("ModelResults.jld","pr",pr,"eq",eq,"tau",tau,"pa",pa);
#pr,eq,tau,pa=load("ModelResults.jld", "pr","eq","tau","pa");


#####################################################################
include("TaxReforms.jl")
include("Reforms.jl")

taucvec = [0.30 0.27 0.24 0.21 0.18 0.15 0.12 0.09 0.06 0.03 0.00];
reform2 = Reform2Vector("reform2.jld", taucvec, pr, eq, tau, pa)
reform3 = Reform3Vector("reform3.jld", taucvec, pr, eq, tau, pa)

taucvec1 = [0.30 0.27 0.24 0.21 0.18 0.17 0.16 0.155 ];
reform1 = Reform1Vector("reform1.jld", taucvec1, pr, eq, tau, pa)
