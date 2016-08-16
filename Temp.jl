cd("C:\\Users\\IF User\\Codes")


using Grid:CoordInterpGrid, BCnan, BCnearest, InterpLinear
using Roots:fzero
using QuantEcon:tauchen
using JLD
using DataFrames
using StatsFuns

include("Main.jl")
include("Firms.jl")
include("FreeEntry.jl")
include("Distribution.jl")
include("Aggregation.jl")
include("SolveSteadyState.jl")
include("TaxReforms.jl")
include("calibrate.jl")
include("Transitions.jl")
include("Simulations.jl")
include("Magnitudes.jl")
include("Reforms.jl")



pa =init_parameters( H=1.176, bbeta=0.972, ff= 0.5, aalphak=0.23, aalphal=0.64, llambda0=0.004, llambda1= 0.04, ddelta = 0.13,
                      allowance=0.86, ttheta = 0.25,rhoz= 0.75, ssigmaz= 0.08, e=0.0, A=1.0);
tau = init_taxes(ttaud =0.15, ttauc= 0.35, ttaui= 0.29, ttaug= 0.15, ttaul=0.28);
@time pr,eq= SolveSteadyState(tau,pa;wguess=0.54, VFItol=10.0^-3.0, verbose=false);
moments=computeMomentsCutoff(eq.E,pr,eq,tau,pa,cutoffCapital=0.0);
save("ModelResults.jld","pr",pr,"eq",eq,"tau",tau,"pa",pa);


initialguess=copy(pr.firmvaluegrid);
pr1, eq1, taunew =taxreform2(0.3, eq, tau, pa; momentsprint=true, verbose=true, firmvalueguess=initialguess);
moments1=computeMomentsCutoff(eq1.E,pr1,eq1,taunew,pa,cutoffCapital=0.0);


taxesb=[0.3 0.25 0.2 0.15 0.1 0.05 0.0];

Reform2Vector("BenchmarkReforms2Below.jld", taxesb, pr, eq, tau, pa)

taxesa=[0.4 0.45 0.5];

Reform2Vector("BenchmarkReforms2Above.jld", taxesa, pr, eq, tau, pa)