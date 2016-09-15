
@everywhere using Grid:CoordInterpGrid, BCnan, BCnearest, InterpLinear
using NPZ
using Sobol

#using DataFrames

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



include("TaxMaximization.jl")
include("Temp.jl")

# 1. Create tax sequence from initial point
  #First component is tauc, second is taui
initialpoint = [0.0,0.28]
taxspace = create_taxspace(initialpoint, 10; jldfile=false)

#npzwrite("taxspace.npy",taxspace)

govexp = 0.046970754282188165;
taul = 0.28 ; wguess0=0.532;
pa =init_parameters(H = 1.094, ddelta=0.07557, ttheta = 0.2290 , rhoz =0.7451, ssigmaz = 0.1067, llambda0 = 0.02605, llambda1 = 0.2467, ff = 1.3856, e=0.01820);
tau0=Taxes(0.0,0.0,0.28,0.0,taul)
close_gov_taue!(govexp,tau0, pa; update=0.75, verbose = true, wguess= 0.53, updateVFIguess = true, outsideparallel= true)

Nexp=3
args = Array(Argument,(Nexp,))
for j=1:Nexp
  tau0 = Taxes(0.0,taxspace[j,1],taxspace[j,2],0.0,taul)
  args[j] = Argument(govexp,tau0,pa)
end

welf_and_taxes = map(taxequilibrium,args);
