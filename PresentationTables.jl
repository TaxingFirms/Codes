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


ref3,pa=load("reform3.jld", "ref","pa");
ref1g,pa=load("reform1g_a.jld", "ref","pa");
pr,eq,tau,pa=load("ModelResults.jld", "pr","eq","tau","pa");
pr0,eq0,tau0= pr,eq,tau;
cev0=0.0;
pr1,eq1,tau1 = ref3[end].pr,ref3[end].eq,ref3[end].tau;
pr2,eq2,tau2 = ref1g[5].pr,ref1g[5].eq,ref1g[5].tau;

cev1=(eq1.a.consumption - eq0.a.consumption)/eq0.a.consumption - pa.H/(eq0.a.consumption*(1+pa.psi))*( (eq1.w*(1-tau1.l)/pa.H)^(1+pa.psi) - (eq0.w*(1-tau0.l)/pa.H)^(1+pa.psi) );
cev2=(eq2.a.consumption - eq0.a.consumption)/eq0.a.consumption - pa.H/(eq0.a.consumption*(1+pa.psi))*( (eq2.w*(1-tau2.l)/pa.H)^(1+pa.psi) - (eq0.w*(1-tau0.l)/pa.H)^(1+pa.psi) );

labels = ["Output","Labor Demand","Consumption","Turnover", "TFP","wage" , "CEV"];
################################################################################
#Counterfactual results
# 1. Capital gains
initialVec=[eq0.a.output, eq0.a.laborsupply, eq0.a.consumption, eq0.E/sum(eq0.distr), eq0.a.output/(eq0.a.capital^pa.alphak*eq0.a.laborsupply^pa.alphal), eq0.w, cev0 ];
finalVec=[eq1.a.output, eq1.a.laborsupply, eq1.a.consumption, eq1.E/sum(eq1.distr), eq1.a.output/(eq1.a.capital^pa.alphak*eq1.a.laborsupply^pa.alphal), eq1.w, cev1 ];
changeVec = finalVec./initialVec-1;

println(DataFrame(Var=labels, initial=initialVec ,final=finalVec, change=changeVec ))

# 1. NO Capital gains
initialVec=[eq0.a.output, eq0.a.laborsupply, eq0.a.consumption, eq0.E/sum(eq0.distr), eq0.a.output/(eq0.a.capital^pa.alphak*eq0.a.laborsupply^pa.alphal), eq0.w, cev0 ];
finalVec2=[eq2.a.output, eq2.a.laborsupply, eq2.a.consumption, eq2.E/sum(eq2.distr), eq2.a.output/(eq2.a.capital^pa.alphak*eq2.a.laborsupply^pa.alphal), eq2.w, cev2 ];
changeVec2 = finalVec./initialVec-1;

println(DataFrame(Var=labels, initial=initialVec ,final=finalVec2, change=changeVec2 ))

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
