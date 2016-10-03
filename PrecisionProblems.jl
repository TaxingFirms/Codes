
@everywhere using Grid:CoordInterpGrid, BCnan, BCnearest, InterpLinear
using JLD
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


ref2,pa2=load("reform2.jld", "ref","pa");

Nref = size(ref2)
cev2 = Array(Float64,Nref);
taucvec = Array(Float64,Nref);
tauevec2 = Array(Float64,Nref);
output =  Array(Float64,Nref);
consumption = Array(Float64,Nref);
wage =  Array(Float64,Nref);
labor = Array(Float64,Nref);
Evec =  Array(Float64,Nref);
sumdistr = Array(Float64,Nref);
welfare =  Array(Float64,Nref);
welfare2 =  Array(Float64,Nref);


for j=1:Nref[1]
  cev2[j] = ref2[j].cev
  taucvec[j] = ref2[j].tau.c
  tauevec2[j] = ref2[j].tau.d
  output[j] = ref2[j].eq.a.output
  consumption[j] = ref2[j].eq.a.consumption
  wage[j] = ref2[j].eq.w
  labor[j] = ref2[j].eq.a.laborsupply
  Evec[j] = ref2[j].eq.E
  sumdistr[j] = sum(ref2[j].eq.distr)
  welfare[j] = ref2[j].eq.a.welfare
  welfare2[j]=log(consumption[j] - (pa2.H/(1+1/pa2.psi))* labor[j]^(1+1/pa2.psi));
end
turnover = Evec./sumdistr

using PyPlot
close("all")

figure()
plot(taucvec,cev2)

figure()
plot(taucvec,welfare2)

figure()
plot(taucvec,tauevec2)

figure()
plot(taucvec,output)

figure()
plot(taucvec,consumption)

figure()
plot(taucvec,labor)

figure()
plot(taucvec,wage)

figure()
plot(taucvec,Evec)

figure()
plot(taucvec,sumdistr)

figure()
plot(taucvec,turnover)
