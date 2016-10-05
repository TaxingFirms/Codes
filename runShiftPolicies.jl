
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

# 1. Run benchmark model in a finer grid for graphs
pa =init_parameters(Nk =500, Nomega =500, H = 3.4, psi=0.5, ddelta = 0.0768, ttheta = 0.2295 , rhoz = 0.75, ssigmaz = 0.0993, llambda0 = 0.0255, llambda1 = 0.23958, ff = 1.44459, e=0.041674 );
tau = init_taxes(ttaud =0.15, ttauc= 0.35, ttaui= 0.28, ttaug= 0.15, ttaul=0.28);
@time pr,eq= SolveSteadyState(tau,pa;wguess=0.5, VFItol=10.0^-5.0, displayit0=true, displayw = true);
moments=computeMomentsCutoff(eq.E,pr,eq,tau,pa,cutoffCapital=0.0;toPrint=true);
save("ModelResults500.jld","pr",pr,"eq",eq,"tau",tau,"pa",pa);


# 2. Initialize equilibrium object: just to have prices
# Important to set r at the initial level
eqaux = init_equilibirium(eq.w,tau,pa);

# 3. Shift policies
## 3.1 Shift tauc to 0.3
pr1  = init_firmproblem(pa, firmvalueguess = pr.firmvaluegrid);
tau1=Taxes(tau.d, 0.3, tau.i, tau.g, tau.l)
firmVFIParallelOmega!(pr1, eqaux, tau1, pa; tol = 10^-5.0 );
getpolicies!(pr1,eqaux,tau1,pa);
#Fix entry and change distributions
dist1 = stationarydist(eq.E, pr1, eqaux, tau1, pa)
#Fix wage consistent with free entry
eq1 = init_equilibirium(eq.w,tau1,pa);
pr11 = deepcopy(pr1);
w1 = free_entry!(eq1, pr11, tau1, pa, firmVFIParallelOmega, maximizationconstraint, 10.0^-5.0);
getpolicies!(pr11,eq1,tau1,pa);
mass_of_entrantsGHH!( pr11, eq1, tau1, pa, stationarydist ; verbose = false);
aggregates!(pr11, eq1, tau1, pa);
##OUTPUT: pr1, dist1, pr1.firmvaluegrid[0], pr11, w1, eq1.E, eq1.distr

save("ShiftTauC.jld","pr1", pr1," dist1", dist1, "valentry" ,pr1.firmvaluegrid[0], "pr11",pr11, "w1", w1, "E", eq1.E, "distr11" ,eq1.distr,"pa",pa);


#  ## 3.2 Shift taud to 0.2
#  pr2  = init_firmproblem(pa, firmvalueguess = pr.firmvaluegrid);
#  tau2=Taxes(0.2, tau.c, tau.i, tau.g, tau.l)
#  firmVFIParallelOmega!(pr2, eqaux, tau2, pa; tol = 10^-5.0 );
#  getpolicies!(pr1,eqaux,tau2,pa);
#  ## 3.3 Shift taui to 0.23
#  pr2  = init_firmproblem(pa, firmvalueguess = pr.firmvaluegrid);
#  tau2=Taxes(tau.d, tau.c, 0.23, tau.g, tau.l)
#  firmVFIParallelOmega!(pr2, eqaux, tau2, pa; tol = 10^-5.0 );
#  ## 3.4 Shift taug to 0.2
#  pr2  = init_firmproblem(pa, firmvalueguess = pr.firmvaluegrid);
#  tau2=Taxes(tau.d, tau.c, tau.i, 0.2, tau.l)
#  firmVFIParallelOmega!(pr2, eqaux, tau2, pa; tol = 10^-5.0 );
