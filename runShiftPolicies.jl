
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
# pr,eq,tau,pa =load("ModelResults500.jld","pr","eq","tau","pa");

# 2. Initialize equilibrium object: just to have prices
# Important to set r at the initial level
eqaux = init_equilibirium(eq.w,tau,pa);

# 3. Shift policies
## 3.1 Shift tauc to 0.3
pr1  = init_firmproblem(pa, firmvalueguess = pr.firmvaluegrid);
tau1=Taxes(tau.d, 0.3, tau.i, tau.g, tau.l)
firmVFIParallelOmega!(pr1, eqaux, tau1, pa; tol = 10^-5.0 );
getpolicies!(pr1,eqaux,tau1,pa);
expvalentry= compute_expvalentry(pr1,pa,eq1,tau1);
#Fix entry and change distributions
dist1 = stationarydist(eq.E, pr1, eqaux, tau1, pa);
#Fix wage consistent with free entry
eq1 = init_equilibirium(eq.w,tau1,pa);
pr11 = deepcopy(pr1);
w1 = free_entry!(eq1, pr11, tau1, pa, firmVFIParallelOmega!, maximizationconstraint, 10.0^-5.0);
getpolicies!(pr11,eq1,tau1,pa);
mass_of_entrantsGHH!( pr11, eq1, tau1, pa, stationarydist ; verbose = false);
aggregates!(pr11, eq1, tau1, pa);
#
dist111 = stationarydist(eq1.E, pr1, eqaux, tau1, pa);

##OUTPUT: pr1, dist1, pr1.firmvaluegrid[0], pr11, w1, eq1.E, eq1.distr

save("ShiftTauC.jld","pr1", pr1,"dist1", dist1, "valentry" ,expvalentry, "pr11",pr11, "w1", w1, "E", eq1.E, "distr11" ,eq1.distr,"pa",pa);

# Generate histogram for the distribution
Nbins =20;
hstep = pa.Nomega/Nbins;
hist = Array(Float64,(Nbins,));
hist1 = Array(Float64,(Nbins,));
hist11 = Array(Float64,(Nbins,));
hdist0= eq.distr[:,7];
hdist1= dist1[:,7];
hdist11 = distr11[:,7]

for j=1:Nbins
  ind0 = convert(Int64,(j-1)*hstep +1);
  indend = convert(Int64,(j)*hstep);
  hist[j] = sum(hdist0[ind0:indend]);
  hist1[j] = sum(hdist1[ind0:indend]);
  hist11[j] = sum(hdist11[ind0:indend]);
end
hist = hist/sum(hist);
hist1 = hist1/sum(hist1);
hist11 = hist11/sum(hist11);


hstep2 = convert(Int64,round(pa.omega.ub/Nbins));
hend = convert(Int64,round(pa.omega.ub))
ind=1:hstep2:hend+1;
width= hstep2*0.4;


figure()
title("Change in Corporate Income Tax \n Partial equilibrium effect", fontsize=16)
xlabel("Net worth", fontsize=14)
k= plot(pa.omega.grid, pr.kpolicy[:,7] , color="r", linewidth = 2.0, label=L"$\tau_c = 0.35$")
k1= plot(pa.omega.grid, pr1.kpolicy[:,7], color="y", linewidth = 2.0, label=L"$\tau_c = 0.3$")
ylabel("""k' """, fontsize=14)
ylim(0,90)
tick_params(labelsize=13)
legend(loc="best")
twinx()
rects = bar(ind -width, hist, width, color="r", label =L"$\tau_c = 0.35$")
rects1 = bar(ind, hist1, width, color="y", label=L"$\tau_c = 0.30$")
ylabel("Mass of firms", fontsize=14)
ylim(0,1)
tick_params(labelsize=13)
savefig("../1_Firmtaxation/1FirmTaxation/Figures/ShiftTauC.pdf")
close()

width= hstep2*0.25;
figure()
title("Change in Corporate Income Tax \n Partial equilibrium effect", fontsize=16)
xlabel("Net worth", fontsize=14)
k= plot(pa.omega.grid, pr.kpolicy[:,7] , color="r", linewidth = 2.0, label=L"$\tau_c = 0.35$")
k1= plot(pa.omega.grid, pr1.kpolicy[:,7], color="y", linewidth = 2.0, label=L"$\tau_c = 0.3$")
k11= plot(pa.omega.grid, pr11.kpolicy[:,7], color="g", linewidth = 2.0, label=L"$\tau_c = 0.3$, GE")
ylabel("""k' """, fontsize=14)
ylim(0,90)
tick_params(labelsize=13)
legend(loc="best")
twinx()
rects = bar(ind -width, hist, width, color="r", label =L"$\tau_c = 0.35$")
rects1 = bar(ind, hist1, width, color="y", label=L"$\tau_c = 0.30$")
rects11 = bar(ind + width, hist11, width, color="g", label=L"$\tau_c = 0.30, GE$")
ylabel("Mass of firms", fontsize=14)
ylim(0,1)
tick_params(labelsize=13)
savefig("../1_Firmtaxation/1FirmTaxation/Figures/ShiftTauC_GE.pdf")
close()



## 3.2 Shift taud to 0.2
pr2  = init_firmproblem(pa, firmvalueguess = pr.firmvaluegrid);
tau2=Taxes(0.2, tau.c, tau.i, tau.g, tau.l)
firmVFIParallelOmega!(pr2, eqaux, tau2, pa; tol = 10^-5.0 );
getpolicies!(pr2,eqaux,tau2,pa);

getpolicies!(pr2,eqaux,tau2,pa);
expvalentry= compute_expvalentry(pr2,pa,eq2,tau2);
#Fix entry and change distributions
dist2 = stationarydist(eq.E, pr2, eqaux, tau2, pa);
#Fix wage consistent with free entry
eq2 = init_equilibirium(eq.w,tau2,pa);
pr22 = deepcopy(pr2);
w2 = free_entry!(eq2, pr22, tau2, pa, firmVFIParallelOmega!, maximizationconstraint, 10.0^-5.0);
getpolicies!(pr22,eq2,tau2,pa);
mass_of_entrantsGHH!( pr22, eq2, tau2, pa, stationarydist ; verbose = false);
aggregates!(pr22, eq2, tau2, pa);
dist222 = stationarydist(eq2.E, pr2, eqaux, tau2, pa);

save("ShiftTauD.jld","pr2", pr2,"dist2", dist2, "eq2" ,eq2, "pr22",pr22, "w2", w2, "dist222" , dist222, "pa",pa);











## 3.3 Shift taui to 0.23
pr3  = init_firmproblem(pa, firmvalueguess = pr.firmvaluegrid);
tau3=Taxes(tau.d, tau.c, 0.23, tau.g, tau.l)
firmVFIParallelOmega!(pr2, eqaux, tau2, pa; tol = 10^-5.0 );
## 3.4 Shift taug to 0.2
pr4  = init_firmproblem(pa, firmvalueguess = pr.firmvaluegrid);
tau4=Taxes(tau.d, tau.c, tau.i, 0.2, tau.l)
firmVFIParallelOmega!(pr2, eqaux, tau2, pa; tol = 10^-5.0 );
