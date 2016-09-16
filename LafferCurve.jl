

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

pr,eq,tau,pa=load("ModelResults.jld", "pr","eq","tau","pa");


taucvec = 0.0:0.05:0.5

~,Ntau=size(taucvec)
collections = Array(Float64,(Ntau,));
collections2gdp = Array(Float64,(Ntau,));

type AllInfo
  pr::FirmProblem
  eq::Equilibrium
end
allinfo = Array(AllInfo,(Ntau,));

benchmarkG=eq.a.G;

for j=1:Ntau
  println(j,j,j,j,j,j,j,j,j,j,j,j,j,j,j)
  pr,eq= SolveSteadyState(Taxes(tau.d,taucvec[j],tau.i,tau.g,tau.l),pa;wguess=eq.w, VFItol=10.0^-3.0, displayit0=false, displayw = false);
  collections[j]= eq.a.G/benchmarkG
  collections2gdp[j]= eq.a.G/eq.a.output
  allinfo[j]= AllInfo(pr,eq)
  println(eq.a.G/benchmarkG)
end

save("LafferTauC.jld","collections",collections,"collections2gdp", collections2gdp,"allinfo",allinfo)
