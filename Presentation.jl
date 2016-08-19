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

pr,eq,tau,pa=load("ModelResults.jld", "pr","eq","tau","pa");

using PyPlot
 rc("text",usetex=true)

 figure()
 d= plot(pa.omega.grid, pr.distributions[:,7], label="distributions")
 k= plot(pa.omega.grid, pr.kpolicy[:,7], label="""k'""")
 q= plot(pa.omega.grid, pr.qpolicy[:,7], label="""b'""")
 xlabel("Net worth")
 title("Policy functions")
 legend(loc="best")
savefig("../1_Firmtaxation/1FirmTaxation/Figures/Policies.pdf")
close()

pa =init_parameters( H=1.176, bbeta=0.972, ff= 0.5, aalphak=0.23, aalphal=0.64, llambda0=0.004, llambda1= 0.04, ddelta = 0.13,
                      allowance=0.86, ttheta = 0.25,rhoz= 0.75, ssigmaz= 0.08, e=0.0, A=1.0, Nk=300, Nq=50);
#Increase precision: Nk=300
firmVFIParallelOmega!(pr,eq,tau,pa;tol=10.0^-3.0, verbose=true);


tau = init_taxes(ttaud =0.15, ttauc= 0.30, ttaui= 0.29, ttaug= 0.15, ttaul=0.28);
pr1=deepcopy(pr);
firmVFIParallelOmega!(pr1,eq,tau,pa;tol=10.0^-3.0, verbose=true);
figure()
k= plot(pa.omega.grid, pr.kpolicy[:,7] , label=L"$\tau_c = 0.35$")
k1= plot(pa.omega.grid, pr1.kpolicy[:,7] ,label=L"$\tau_c = 0.3$")
  xlabel("Net worth")
  ylabel("""k' """)
  title("Capital policy")
  legend()
savefig("../1_Firmtaxation/1FirmTaxation/Figures/ShiftTauC.pdf")
close()


tau = init_taxes(ttaud =0.14, ttauc= 0.35, ttaui= 0.29, ttaug= 0.15, ttaul=0.28);
pr2=deepcopy(pr);
firmVFIParallelOmega!(pr2,eq,tau,pa;tol=10.0^-3.0, verbose=true);
figure()
k= plot(pa.omega.grid, pr.kpolicy[:,7] , label=L"$\tau_c = 0.15$")
k1= plot(pa.omega.grid, pr2.kpolicy[:,7] ,label=L"$\tau_c = 0.1$")
  xlabel("Net worth")
  ylabel("""k' """)
  title("Capital policy")
  legend()
  savefig("../1_Firmtaxation/1FirmTaxation/Figures/ShiftTauD.pdf")
  close()


tau = init_taxes(ttaud =0.10, ttauc= 0.35, ttaui= 0.29, ttaug= 0.10, ttaul=0.28);
pr3=deepcopy(pr);
firmVFIParallelOmega!(pr3,eq,tau,pa;tol=10.0^-3.0, verbose=true);
figure()
k= plot(pa.omega.grid, pr.kpolicy[:,7] , label=L"$\tau_c = \tau_g = 0.15$")
k1= plot(pa.omega.grid, pr3.kpolicy[:,7] ,label=L"$\tau_c = \tau_g = 0.1$")
xlabel("Net worth")
ylabel("""k' """)
title("Capital policy")
legend()
savefig("../1_Firmtaxation/1FirmTaxation/Figures/ShiftTauG.pdf")
close()


tau = init_taxes(ttaud =0.15, ttauc= 0.35, ttaui= 0.24, ttaug= 0.15, ttaul=0.28);
pr4=deepcopy(pr);
firmVFIParallelOmega!(pr4,eq,tau,pa;tol=10.0^-3.0, verbose=true);
figure()
k= plot(pa.omega.grid, pr.kpolicy[:,7] , label=L"$\tau_i = 0.29$")
k1= plot(pa.omega.grid, pr4.kpolicy[:,7] ,label=L"$\tau_i = 0.24$")
xlabel("Net worth")
ylabel("""k' """)
title("Capital policy")
legend()
savefig("../1_Firmtaxation/1FirmTaxation/Figures/ShiftTauI.pdf")
close()


save("ShiftPolicies", "pr", pr, "pr1", pr1, "pr2", pr2, "pr3", pr3, "pr4", pr4)
