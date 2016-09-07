# All you need to run the code
# Modify ~/.juliarc.jl and add the following line:
# push!(LOAD_PATH, "/Path/To/My/Module/")

@everywhere using Grid:CoordInterpGrid, BCnan, BCnearest, InterpLinear
using JLD
using DataFrames
using StatsFuns

include("markov_approx.jl")
include("mc_tools.jl")

@everywhere include("Main.jl")
@everywhere include("Firms.jl")
@everywhere include("FreeEntry.jl")
@everywhere include("Distribution.jl")
@everywhere include("Aggregation.jl")
@everywhere include("SolveSteadyState.jl")
@everywhere include("TaxReforms.jl")
@everywhere include("calibrate.jl")
@everywhere include("Transitions.jl")


include("Simulations.jl")
include("Magnitudes.jl")
include("PlotFunctions.jl")
include("Reforms.jl")



#pa =init_parameters(H=1.07, ddelta=0.125,ttheta=0.2, llambda1 = 0.2, rhoz = 0.75, ssigmaz = 0.072  ,ff= 0.5, e= 0.16);

pa =init_parameters(H = 1.2135, ddelta=0.0768, ttheta = 0.274 , rhoz =0.75, ssigmaz = 0.086, llambda0 = 0.03, llambda1 = 0.1507, ff = 0.796, e=0.0388);

tau = init_taxes(ttaud =0.15, ttauc= 0.35, ttaui= 0.28, ttaug= 0.15, ttaul=0.28);
@time pr,eq= SolveSteadyState(tau,pa;wguess=0.556, VFItol=10.0^-3.0,displayit0=false, displayw=false);
moments=computeMomentsCutoff(eq.E,pr,eq,tau,pa,cutoffCapital=0.0;toPrint=false);
pcterror_params(pr,eq,tau,pa)
save("ModelResults.jld","pr",pr,"eq",eq,"tau",tau,"pa",pa);
#pr,eq,tau,pa=load("ModelResults.jld", "pr","eq","tau","pa");
