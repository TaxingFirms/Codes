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


figure()
k= plot(pa.omega.grid, pr.kpolicy[:,7] , label=L"$\tau_c = 0.35$")
k1= plot(pa.omega.grid, pr1.kpolicy[:,7] ,label=L"$\tau_c = 0.3$")
  xlabel("Net worth")
  ylabel("""k' """)
  title("Capital policy")
  legend()
savefig("../1_Firmtaxation/1FirmTaxation/Figures/ShiftTauC.pdf")
close()


figure()
k= plot(pa.omega.grid, pr.kpolicy[:,7] , label=L"$\tau_c = 0.15$")
k1= plot(pa.omega.grid, pr2.kpolicy[:,7] ,label=L"$\tau_c = 0.1$")
  xlabel("Net worth")
  ylabel("""k' """)
  title("Capital policy")
  legend()
  savefig("../1_Firmtaxation/1FirmTaxation/Figures/ShiftTauD.pdf")
  close()


figure()
k= plot(pa.omega.grid, pr.kpolicy[:,7] , label=L"$\tau_c = \tau_g = 0.15$")
k1= plot(pa.omega.grid, pr3.kpolicy[:,7] ,label=L"$\tau_c = \tau_g = 0.1$")
xlabel("Net worth")
ylabel("""k' """)
title("Capital policy")
legend()
savefig("../1_Firmtaxation/1FirmTaxation/Figures/ShiftTauG.pdf")
close()


figure()
k= plot(pa.omega.grid, pr.kpolicy[:,7] , label=L"$\tau_i = 0.29$")
k1= plot(pa.omega.grid, pr4.kpolicy[:,7] ,label=L"$\tau_i = 0.24$")
xlabel("Net worth")
ylabel("""k' """)
title("Capital policy")
legend()
savefig("../1_Firmtaxation/1FirmTaxation/Figures/ShiftTauI.pdf")
close()
