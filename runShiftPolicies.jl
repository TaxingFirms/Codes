
@everywhere using Grid:CoordInterpGrid, BCnan, BCnearest, InterpLinear
using JLD
using DataFrames
using PyPlot

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
#  pa =init_parameters(Nk =500, Nomega =500, H = 3.4, psi=0.5, ddelta = 0.0768, ttheta = 0.2295 , rhoz = 0.75, ssigmaz = 0.0993, llambda0 = 0.0255, llambda1 = 0.23958, ff = 1.44459, e=0.041674 );
#  tau = init_taxes(ttaud =0.15, ttauc= 0.35, ttaui= 0.28, ttaug= 0.15, ttaul=0.28);
#  @time pr,eq= SolveSteadyState(tau,pa;wguess=0.5, VFItol=10.0^-5.0, displayit0=true, displayw = true);
#  moments=computeMomentsCutoff(eq.E,pr,eq,tau,pa,cutoffCapital=0.0;toPrint=true);
#  save("ModelResults500.jld","pr",pr,"eq",eq,"tau",tau,"pa",pa);
pr,eq,tau,pa =load("ModelResults500.jld","pr","eq","tau","pa");

# Generate histogram for the initial distribution
indomegamax = 360;
Nbins =20;
hstep = indomegamax/Nbins;
hist = Array(Float64,(Nbins,));
hdist0= sum(eq.distr,2);#eq.distr[:,7];

proddisplay = 5 ;

for j=1:Nbins
  ind0 = convert(Int64,(j-1)*hstep +1);
  indend = convert(Int64,(j)*hstep);
  hist[j] = sum(hdist0[ind0:indend]);
end
hist = hist/sum(hist);


# 2. Initialize equilibrium object: just to have prices
# Important to set r at the initial level
eqaux = init_equilibirium(eq.w,tau,pa);

# 3. Shift policies
## 3.1 Shift tauc to 0.3
pr1  = init_firmproblem(pa, firmvalueguess = pr.firmvaluegrid);
tau1=Taxes(tau.d, 0.3, tau.i, tau.g, tau.l)
firmVFIParallelOmega!(pr1, eqaux, tau1, pa; tol = 10^-5.0 );
getpolicies!(pr1,eqaux,tau1,pa); #save("Temp.jld","pr1",pr1)
expvalentry1= compute_expvalentry(pr1,pa,eq1,tau1);
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

save("ShiftTauC.jld","pr1", pr1,"dist1", dist1, "eq1" ,eq1, "pr11",pr1, "w1", w1, "dist111" , dist111, "pa",pa);
# pr1, dist1,eq1 ,pr11, w1, dist111,pa=load("ShiftTauC.jld","pr1","dist1", "eq1", "pr11", "w1", "dist111", "pa");

# Generate histograms fr counterfactuals
hist1 = Array(Float64,(Nbins,));
hist11 = Array(Float64,(Nbins,));
hdist1= sum(dist1,2); #dist1[:,7];
hdist11 = sum(eq1.distr,2); #distr11[:,7];

for j=1:Nbins
  ind0 = convert(Int64,(j-1)*hstep +1);
  indend = convert(Int64,(j)*hstep);
  hist1[j] = sum(hdist1[ind0:indend]);
  hist11[j] = sum(hdist11[ind0:indend]);
end
hist1 = hist1/sum(hist1);
hist11 = hist11/sum(hist11);


hstep2 = pa.omega.grid[indomegamax]/Nbins;
hend = pa.omega.grid[indomegamax] - 0.5*hstep2;
ind=0.5*hstep2:hstep2:hend;
width= hstep2*0.4;


figure()
title("Change in Corporate Income Tax \n Partial equilibrium effect", fontsize=16)
xlabel("Net worth", fontsize=14)
k= plot(pa.omega.grid[1:indomegamax], pr.kpolicy[1:indomegamax,proddisplay] , color="r", linewidth = 2.0, label=L"$\tau_c = 0.35$")
k1= plot(pa.omega.grid[1:indomegamax], pr1.kpolicy[1:indomegamax,proddisplay], color="y", linewidth = 2.0, label=L"$\tau_c = 0.30$")
ylabel("""k' """, fontsize=14)
ylim(0,35)
tick_params(labelsize=13)
legend(loc="center right")
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
title("Change in Corporate Income Tax \n General equilibrium effect", fontsize=16)
xlabel("Net worth", fontsize=14)
k= plot(pa.omega.grid[1:indomegamax], pr.kpolicy[1:indomegamax,proddisplay] , color="r", linewidth = 2.0, label=L"$\tau_c = 0.35$")
k1= plot(pa.omega.grid[1:indomegamax], pr1.kpolicy[1:indomegamax,proddisplay], color="y", linewidth = 2.0, label=L"$\tau_c = 0.3$")
k11= plot(pa.omega.grid[1:indomegamax], pr11.kpolicy[1:indomegamax,proddisplay], color="g", linewidth = 2.0, label=L"$\tau_c = 0.3$, GE")
ylabel("""k' """, fontsize=14)
ylim(0,35)
tick_params(labelsize=13)
legend(loc="center right")
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
expvalentry2= compute_expvalentry(pr2,pa,eq2,tau2);
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

#pr2, dist2, eq2, pr22, w2, dist222, pa = load("ShiftTauD.jld","pr2", "dist2", "eq2", "pr22","w2", "dist222" , "pa");

hist2 = Array(Float64,(Nbins,));
hist22 = Array(Float64,(Nbins,));
hist222 = Array(Float64,(Nbins,));
hdist2 = sum(dist2,2); #dist2[:,7];
hdist22 = sum(eq2.distr,2); #eq2.distr[:,7];
hdist222 = sum(dist222,2); #dist222[:,7];


for j=1:Nbins
  ind0 = convert(Int64,(j-1)*hstep +1);
  indend = convert(Int64,(j)*hstep);
  hist2[j] = sum(hdist2[ind0:indend]);
  hist22[j] = sum(hdist22[ind0:indend]);
  hist222[j] = sum(hdist222[ind0:indend]);
end
hist2 = hist2/sum(hist2);
hist22 = hist22/sum(hist22);
hist222 = hist222/sum(hist222);

width= hstep2*0.4;


figure()
title("Change in Dividend Tax \n Partial equilibrium effect", fontsize=16)
xlabel("Net worth", fontsize=14)
k= plot(pa.omega.grid[1:indomegamax], pr.kpolicy[1:indomegamax,proddisplay] , color="r", linewidth = 2.0, label=L"$\tau_d = 0.15$")
k1= plot(pa.omega.grid[1:indomegamax], pr2.kpolicy[1:indomegamax,proddisplay], color="b", linewidth = 2.0, label=L"$\tau_d = 0.20$")
ylabel("""k' """, fontsize=14)
#ylim(-20,90)
tick_params(labelsize=13)
legend(loc="best")
twinx()
rects = bar(ind -width, hist, width, color="r", label =L"$\tau_d = 0.15$")
rects2 = bar(ind, hist2, width, color="b", label=L"$\tau_d = 0.20$")
ylabel("Mass of firms", fontsize=14)
ylim(0,1)
tick_params(labelsize=13)
savefig("../1_Firmtaxation/1FirmTaxation/Figures/ShiftTauD.pdf")
close()

expvalentry2= compute_expvalentry(pr2,pa,eq2,tau2);


width= hstep2*0.25;
figure()
title("Change in Dividend Tax \n General equilibrium effect", fontsize=16)
xlabel("Net worth", fontsize=14)
k=plot(pa.omega.grid[1:indomegamax], pr.kpolicy[1:indomegamax,proddisplay] , color="r", linewidth = 2.0, label=L"$\tau_d = 0.15$")
k2= plot(pa.omega.grid[1:indomegamax], pr2.kpolicy[1:indomegamax,proddisplay], color="b", linewidth = 2.0, label=L"$\tau_d = 0.20$")
k2= plot(pa.omega.grid[1:indomegamax], pr22.kpolicy[1:indomegamax,proddisplay], color="c", linewidth = 2.0, label=L"$\tau_d = 0.20, GE$")
ylabel("""k' """, fontsize=14)
#ylim(-20,90)
tick_params(labelsize=13)
legend(loc="best")
twinx()
ects = bar(ind -width, hist, width, color="r", label =L"$\tau_d = 0.15$")
rects2 = bar(ind, hist2, width, color="b", label=L"$\tau_d = 0.20$")
rects22 = bar(ind + width, hist22, width, color="c", label=L"$\tau_c = 0.30, GE$")
ylabel("Mass of firms", fontsize=14)
ylim(0,1)
tick_params(labelsize=13)
savefig("../1_Firmtaxation/1FirmTaxation/Figures/ShiftTauD_GE.pdf")
close()




## 3.3 Shift taud to 0.2, taug to 0.2
pr3  = init_firmproblem(pa, firmvalueguess = pr.firmvaluegrid);
tau3=Taxes(0.2, tau.c, tau.i, 0.2, tau.l)
firmVFIParallelOmega!(pr3, eqaux, tau3, pa; tol = 10^-5.0 );
getpolicies!(pr3,eqaux,tau3,pa);

getpolicies!(pr3,eqaux,tau3,pa);
expvalentry= compute_expvalentry(pr3,pa,eq3,tau3);
#Fix entry and change distributions
dist3 = stationarydist(eq.E, pr3, eqaux, tau3, pa);
#Fix wage consistent with free entry
eq3 = init_equilibirium(eq.w,tau3,pa);
pr33 = deepcopy(pr3);
w3 = free_entry!(eq3, pr33, tau3, pa, firmVFIParallelOmega!, maximizationconstraint, 10.0^-5.0);
getpolicies!(pr33,eq3,tau3,pa);
mass_of_entrantsGHH!( pr33, eq3, tau3, pa, stationarydist ; verbose = false);
aggregates!(pr33, eq3, tau3, pa);

save("ShiftTauDG.jld","pr3", pr3,"dist3", dist3, "eq3" ,eq3, "pr33", pr33, "w3", w3,  "pa",pa);
pr3, dist3, eq3, pr33, w3, pa = load("ShiftTauDG.jld","pr3","dist3","eq3", "pr33", "w3", "pa");

hist3 = Array(Float64,(Nbins,));
hist33 = Array(Float64,(Nbins,));
hist333 = Array(Float64,(Nbins,));
hdist3 = sum(dist2,2); #dist2[:,7];
hdist33 = sum(eq2.distr,2); #eq2.distr[:,7];

for j=1:Nbins
  ind0 = convert(Int64,(j-1)*hstep +1);
  indend = convert(Int64,(j)*hstep);
  hist3[j] = sum(hdist3[ind0:indend]);
  hist33[j] = sum(hdist33[ind0:indend]);
end
hist3 = hist3/sum(hist3);
hist33 = hist33/sum(hist33);

width= hstep2*0.4;


figure()
title("Change in Dividend Tax and Capital Gains\n Partial equilibrium effect", fontsize=16)
xlabel("Net worth", fontsize=14)
k= plot(pa.omega.grid[1:indomegamax], pr.kpolicy[1:indomegamax,proddisplay] , color="r", linewidth = 2.0, label=L"$\tau_d = \tau_g = 0.15$")
k1= plot(pa.omega.grid[1:indomegamax], pr3.kpolicy[1:indomegamax,proddisplay], color="m", linewidth = 2.0, label=L"$\tau_d = \tau_g = 0.20$")
ylabel("""k' """, fontsize=14)
#ylim(-20,90)
tick_params(labelsize=13)
legend(loc="best")
twinx()
rects = bar(ind -width, hist, width, color="r", label =L"$\tau_d = \tau_g = 0.15$")
rects2 = bar(ind, hist3, width, color="m", label=L"$\tau_d = \tau_g = 0.15$")
ylabel("Mass of firms", fontsize=14)
ylim(0,1)
tick_params(labelsize=13)
savefig("../1_Firmtaxation/1FirmTaxation/Figures/ShiftTauDG.pdf")
close()

expvalentry3= compute_expvalentry(pr3,pa,eq3,tau3);



width= hstep2*0.25;
figure()
title("Change in Dividend Tax \n General equilibrium effect", fontsize=16)
xlabel("Net worth", fontsize=14)
k=plot(pa.omega.grid[1:indomegamax], pr.kpolicy[1:indomegamax,proddisplay] , color="r", linewidth = 2.0, label=L"$\tau_d = \tau_g = 0.15$")
k3= plot(pa.omega.grid[1:indomegamax], pr3.kpolicy[1:indomegamax,proddisplay], color="m", linewidth = 2.0, label=L"$\tau_d = 0.20$")
k33= plot(pa.omega.grid[1:indomegamax], pr33.kpolicy[1:indomegamax,proddisplay], color="y", linewidth = 2.0, label=L"$\tau_d = 0.20, GE$")
ylabel("""k' """, fontsize=14)
#ylim(-20,90)
tick_params(labelsize=13)
legend(loc="best")
twinx()
rects = bar(ind -width, hist, width, color="r", label =L"$\tau_d = 0.15$")
rects3 = bar(ind, hist3, width, color="m", label=L"$\tau_d = 0.20$")
rects33 = bar(ind + width, hist33, width, color="y", label=L"$\tau_c = 0.30, GE$")
ylabel("Mass of firms", fontsize=14)
ylim(0,1)
tick_params(labelsize=13)
savefig("../1_Firmtaxation/1FirmTaxation/Figures/ShiftTauD_GE.pdf")
close()
