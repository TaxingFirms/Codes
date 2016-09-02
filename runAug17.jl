
@everywhere using Grid:CoordInterpGrid, BCnan, BCnearest, InterpLinear
@everywhere using Roots:fzero
using QuantEcon:tauchen
using JLD
using DataFrames
using StatsFuns

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
include("Reforms.jl")

#See August13jl for parameter selection
pa =init_parameters( H=1.176, bbeta=0.972, ff= 0.5, aalphak=0.23, aalphal=0.64, llambda0=0.004, llambda1= 0.04, ddelta = 0.13,
                      allowance=0.86, ttheta = 0.25,rhoz= 0.75, ssigmaz= 0.08, e=0.0, A=1.0);
tau = init_taxes(ttaud =0.15, ttauc= 0.35, ttaui= 0.29, ttaug= 0.15, ttaul=0.28);
@time pr,eq= SolveSteadyState(tau,pa;wguess=0.54, VFItol=10.0^-3.0, verbose=false);
moments=computeMomentsCutoff(eq.E,pr,eq,tau,pa,cutoffCapital=0.0;toPrint=false);
pcterror_params(pr,eq,tau,pa)
save("Model1.jld","pr",pr,"eq",eq,"tau",tau,"pa",pa);
#pr,eq,tau,pa=load("ModelResults.jld", "pr","eq","tau","pa");

#taxesb=[0.33 0.31 0.29 0.27 0.25 0.23 0.21 0.19 0.17 0.15 0.13 0.11 0.09 0.07 0.05 0.03 0.01];
taxesb=[0.3 0.25 0.2 0.15 0.1 0.05 0.0];
#taxesa=[0.4 0.45 0.5];

ref=Reform2Vector("Benchmark1Reforms2.jld", taxesb, pr, eq, tau, pa)


pa =init_parameters( H=1.176, bbeta=0.972, ff= 0.5, aalphak=0.23, aalphal=0.64, llambda0=0.006, llambda1= 0.04, ddelta = 0.1,
                      allowance=1.00, ttheta = 0.25,rhoz= 0.75, ssigmaz= 0.09, e=0.0, A=1.0);
tau = init_taxes(ttaud =0.15, ttauc= 0.35, ttaui= 0.29, ttaug= 0.15, ttaul=0.28);
@time pr,eq= SolveSteadyState(tau,pa;wguess=0.54, VFItol=10.0^-3.0, verbose=false);
moments=computeMomentsCutoff(eq.E,pr,eq,tau,pa,cutoffCapital=0.0;toPrint=false);
pcterror_params(pr,eq,tau,pa)
ref=Reform2Vector("Benchmark2Reforms2.jld", taxesb, pr, eq, tau, pa)
