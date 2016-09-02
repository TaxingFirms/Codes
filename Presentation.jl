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

pr,pr1,pr2,pr3,pr4,pr5,pa=load("ShiftPolicies.jld", "pr", "pr1", "pr2", "pr3", "pr4","pr5","pa");
~,eq,tau,~=load("ModelResults.jld", "pr","eq","tau","pa");
getpolicies!(pr,eq,tau,pa);  #Update distributions on the finer grid

using PyPlot
 rc("text",usetex=true)
close("all")


 figure()
 d= plot(pa.omega.grid, pr.distributions[:,7], label="distributions")
 k= plot(pa.omega.grid, pr.kpolicy[:,7], label="""k'""")
 q= plot(pa.omega.grid, pr.qpolicy[:,7], label="""b'""")
 xlabel("Net worth")
 title("Policy functions")
 legend(loc="best")
savefig("../1_Firmtaxation/1FirmTaxation/Figures/Policies.pdf")
close()


figure()
k= plot(pa.omega.grid, pr.kpolicy[:,7] , label=L"$\tau_c = 0.35$")
k1= plot(pa.omega.grid, pr1.kpolicy[:,7] ,label=L"$\tau_c = 0.3$")
  xlabel("Net worth")
  ylabel("""k' """)
  title("Capital policy")
  legend(loc="best")
savefig("../1_Firmtaxation/1FirmTaxation/Figures/ShiftTauC.pdf")
close()


figure()
k= plot(pa.omega.grid, pr.kpolicy[:,7] , label=L"$\tau_d = 0.15$")
k1= plot(pa.omega.grid, pr2.kpolicy[:,7] ,label=L"$\tau_d = 0.20$")
  xlabel("Net worth")
  ylabel("""k' """)
  title("Capital policy")
  legend(loc="best")
  savefig("../1_Firmtaxation/1FirmTaxation/Figures/ShiftTauD.pdf")
  close()


figure()
k= plot(pa.omega.grid, pr.kpolicy[:,7] , label=L"$\tau_c = \tau_g = 0.15$")
k1= plot(pa.omega.grid, pr3.kpolicy[:,7] ,label=L"$\tau_c = \tau_g = 0.1$")
xlabel("Net worth")
ylabel("""k' """)
title("Capital policy")
legend(loc="best")
savefig("../1_Firmtaxation/1FirmTaxation/Figures/ShiftTauG.pdf")
close()


figure()
k= plot(pa.omega.grid, pr.kpolicy[:,7] , label=L"$\tau_i = 0.29$")
k1= plot(pa.omega.grid, pr4.kpolicy[:,7] ,label=L"$\tau_i = 0.24$")
xlabel("Net worth")
ylabel("""k' """)
title("Capital policy")
legend(loc="best")
savefig("../1_Firmtaxation/1FirmTaxation/Figures/ShiftTauI.pdf")
close()


figure()
k= plot(pa.omega.grid, pr.kpolicy[:,7] , label=L"$\tau_c = 0.35, \tau_d=\tau_g= 0.15$")
k1= plot(pa.omega.grid, pr5.kpolicy[:,7] ,label=L"$\tau_c = 0.30, \tau_d=\tau_g= 0.20$")
xlabel("Net worth")
ylabel("""k' """)
title("Capital policy")
legend(loc="best")
savefig("../1_Firmtaxation/1FirmTaxation/Figures/ShiftTauI.pdf")
close()


##################################################################
T=19;
z_history_ind = 7*ones(Int,(T,1));

capital, debt, networth, dividends, investment=timeseries(z_history_ind,T,pr,pa);
capital1, debt1, networth1, dividends1, investment1=timeseries(z_history_ind,T,pr1,pa);
capital2, debt2, networth2, dividends2, investment2=timeseries(z_history_ind,T,pr2,pa);
capital3, debt3, networth3, dividends3, investment3=timeseries(z_history_ind,T,pr3,pa);
capital4, debt4, networth4, dividends4, investment4=timeseries(z_history_ind,T,pr4,pa);
capital5, debt5, networth5, dividends5, investment5=timeseries(z_history_ind,T,pr5,pa);


figure()
plot(1:T-1, capital[2:T],"+-", label="Benchmark")
plot(1:T-1, capital1[2:T],"+-", label=L"$\tau_c = 0.3$")
plot(1:T-1, capital3[2:T],"+-", label=L"$\tau_d = \tau_g = 0.1")
xlabel("Period")
ylabel("""k' """)
legend(loc="best")
savefig("../1_Firmtaxation/1FirmTaxation/Figures/TimeSeries7.pdf")
close()


##################################################################
T=500;
S=5000;
simcap, simdebt, simnetworth, simdividends, siminvestment, z_history_ind, survivor= simulation(S, T,pr,pa; seed=1234);
simcap1, ~, ~, ~, ~, ~, survivor1= simulation(S, T,pr1,pa; seed=1234);
simcap3, ~, ~, ~, ~, ~, survivor3= simulation(S, T,pr3,pa; seed=1234);
simcap5, ~, ~, ~, ~, ~, survivor5= simulation(S, T,pr5,pa; seed=1234);
simcap_bar= sum(simcap,2)./sum(survivor,2);
simcap1_bar= sum(simcap1,2)./sum(survivor1,2);
simcap3_bar= sum(simcap3,2)./sum(survivor3,2);
simcap5_bar= sum(simcap5,2)./sum(survivor5,2);

lifetime=sum(survivor,1);
ave_lifetime=sum(lifetime)/S;
Tmax=round(Int64,ave_lifetime);

figure()
plot(1:Tmax, simcap_bar[2:Tmax+1],"+-", label=L"Benchmark")
plot(1:Tmax, simcap1_bar[2:Tmax+1],"+-", label=L"$\tau_c = 0.30 $")
plot(1:Tmax, simcap3_bar[2:Tmax+1],"+-", label=L"$\tau_d = \tau_g = 0.10 $")
xlabel("Period")
ylabel("""k""")
#title("Time Series of Capital (5000 simulations average)")
legend(loc="lower right")
xlim(1,Tmax)
savefig("../1_Firmtaxation/1FirmTaxation/Figures/TimeSeriesSimul.pdf")


##################################################################
