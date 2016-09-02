using Grid:CoordInterpGrid, BCnan, BCnearest, InterpLinear
using Roots:fzero
using QuantEcon:tauchen
using JLD
using DataFrames
using StatsFuns
using PyPlot

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


#ref,pa=load("BenchmarkReforms3b.jld", "ref","pa");
#pr0,eq0,tau0,cev0 = ref[1].pr,ref[1].eq,ref[1].tau, ref[1].cev;
#pr1,eq1,tau1,cev1 = ref[end].pr,ref[end].eq,ref[end].tau, ref[end].cev;
cd("/Users/danielwillsr/Dropbox/TeslaResults")
rpr,req,rtau,pa=load("Reforms5.jld", "rpr","req","rtau","pa");
cd("/Users/danielwillsr/Dropbox/Codes")
pr,eq,tau,pa=load("ModelResults.jld", "pr","eq","tau","pa");
pr0,eq0,tau0= pr,eq,tau;
cev0=0.0;
pr1,eq1,tau1 = rpr,req,rtau;
cev1=(req.a.consumption - eq.a.consumption)/eq.a.consumption - pa.H/(eq.a.consumption*(1+pa.psi))*( (req.w*(1-rtau.l)/pa.H)^(1+pa.psi) - (eq.w*(1-tau.l)/pa.H)^(1+pa.psi) );
labels = ["Output","Labor Demand","Consumption","Turnover", "TFP", "CEV"];
################################################################################
#Counterfactual results
initialVec=[eq0.a.output, eq0.a.laborsupply, eq0.a.consumption, eq0.E/sum(eq0.distr), eq0.a.output/(eq0.a.capital^pa.alphak*eq0.a.laborsupply^pa.alphal), cev0 ];
finalVec=[eq1.a.output, eq1.a.laborsupply, eq1.a.consumption, eq1.E/sum(eq1.distr), eq1.a.output/(eq1.a.capital^pa.alphak*eq1.a.laborsupply^pa.alphal), cev1 ];
changeVec = finalVec./initialVec-1;

println(DataFrame(Var=labels, initial=initialVec ,final=finalVec, change=changeVec ))
################################################################################
#productivity histograms
hist0=sum(eq0.distr/sum(eq0.distr),1)
hist0p=reshape(hist0,(9,))

hist1=sum(eq1.distr/sum(eq1.distr),1)
hist1p=reshape(hist1,(9,))

ind=1:pa.Nz;
width=0.35;
figure()
rects1 = bar(ind -width, hist0p, width, color="r", label =L"$\tau_c = 0.35$")
recst2 =bar(ind, hist1p, width, color="y", label=L"$\tau_c = 0.0$")
ticklabels = round(pa.zgrid*100)/100
xticks(ind, ticklabels )
xlabel("Productivity")
ylabel("Fraction of capital")
title("Distribution of Productivity")
legend(loc="best")
savefig("../1_Firmtaxation/1FirmTaxation/Figures/ProductivityHist.pdf")


################################################################################
include("Decomposition.jl")

decompostion(pr0::FirmProblem, eq0::Equilibrium, tau0::Taxes, pr1::FirmProblem, eq1::Equilibrium, tau1::Taxes, pa::Param)
  #Decompose the effects of a policy change

################################################################################
