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
pa =init_parameters(H = 3.4, psi=0.5, ddelta = 0.0768, ttheta = 0.2295 , rhoz = 0.75, ssigmaz = 0.0993, llambda0 = 0.0255, llambda1 = 0.23958, ff = 1.44459, e=0.041674);

tau = init_taxes(ttaud =0.15, ttauc= 0.35, ttaui= 0.28, ttaug= 0.15, ttaul=0.28);
@time pr,eq= SolveSteadyState(tau,pa;wguess=0.5, VFItol=10.0^-5.0, displayit0=true, displayw = true);
moments=computeMomentsCutoff(eq.E,pr,eq,tau,pa,cutoffCapital=0.0;toPrint=true);
save("ModelResults.jld","pr",pr,"eq",eq,"tau",tau,"pa",pa);
#pr,eq,tau,pa=load("ModelResults.jld", "pr","eq","tau","pa");


# NO TAU G
pa =init_parameters(H = 3.4, psi=0.5, ddelta = 0.0769, ttheta = 0.22956 , rhoz = 0.75, ssigmaz = 0.10608, llambda0 = 0.0456, llambda1 = 0.08897, ff = 1.3267, e=0.0260);

tau = init_taxes(ttaud =0.15, ttauc= 0.35, ttaui= 0.28, ttaug= 0.0, ttaul=0.28);
@time pr,eq= SolveSteadyState(tau,pa;wguess=0.5, VFItol=10.0^-3.0, displayit0=false, displayw = false);
moments=computeMomentsCutoff(eq.E,pr,eq,tau,pa,cutoffCapital=0.0;toPrint=true);
save("ModelResultsNoTaxG.jld","pr",pr,"eq",eq,"tau",tau,"pa",pa);
