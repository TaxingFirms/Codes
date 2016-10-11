#0. Load stuff
@everywhere using Grid:CoordInterpGrid, BCnan, BCnearest, InterpLinear
using JLD

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

pr,eq,tau,pa =load("ModelResultsNoTaxG600.jld","pr","eq","tau","pa");

#1. Initialize vector of type pr, tau,
taudvec= [0.15 0.20 0.25 0.30 0.35 0.40 0.45]
Nvec=size(taudvec)[2];
prvec=Array(FirmProblem,size(taudvec))
for j=1:Nvec
  prvec[1,j]  = init_firmproblem(pa, firmvalueguess = pr.firmvaluegrid);
end
eqaux = init_equilibirium(eq.w,tau,pa);

#2.  compute value fcn for different levels of taud

for j=1:Nvec
  taux=Taxes(taudvec[j], tau.c, tau.i, tau.g, tau.l)
  firmVFIParallelOmega!(prvec[j], eqaux, taux, pa; tol = 10^-5.0 );
  getpolicies!(prvec[j],eqaux,taux,pa); #save("Temp.jld","pr1",pr1)
end

save("ShiftTauDPolicy.jld","prvec",prvec)

#3. same for tauc

#clean prvec
for j=1:Nvec
  prvec[1,j]  = init_firmproblem(pa, firmvalueguess = pr.firmvaluegrid);
end

taucvec= [0.20 0.25 0.30 0.35 0.40 0.45 0.50]
Nvec=size(taucvec)[2];

for j=1:Nvec
  taux=Taxes(tau.d,taucvec[j], tau.i, tau.g, tau.l)
  firmVFIParallelOmega!(prvec[j], eqaux, taux, pa; tol = 10^-5.0 );
  getpolicies!(prvec[j],eqaux,taux,pa); #save("Temp.jld","pr1",pr1)
end

save("ShiftTauCPolicy.jld","prvec",prvec)

#3. same for taud and taug at the same time

#clean prvec
for j=1:Nvec
  prvec[1,j]  = init_firmproblem(pa, firmvalueguess = pr.firmvaluegrid);
end

taudvec= [0.15 0.20 0.25 0.30 0.35 0.40 0.45]
Nvec=size(taudvec)[2];

for j=1:Nvec
  taux=Taxes(taudvec[j], tau.c, tau.i, taudvec[j], tau.l)
  firmVFIParallelOmega!(prvec[j], eqaux, taux, pa; tol = 10^-5.0 );
  getpolicies!(prvec[j],eqaux,taux,pa); #save("Temp.jld","pr1",pr1)
end

save("ShiftTauDGPolicy.jld","prvec",prvec)

#4. Now taui

for j=1:Nvec
  prvec[1,j]  = init_firmproblem(pa, firmvalueguess = pr.firmvaluegrid);
end

tauivec= [0.15 0.20 0.25 0.30 0.35 0.40 0.45]
Nvec=size(tauivec)[2];

for j=1:Nvec
  taux=Taxes(tau.d, tau.c, tauivec[j], tau.g, tau.l)
  firmVFIParallelOmega!(prvec[j], eqaux, taux, pa; tol = 10^-5.0 );
  getpolicies!(prvec[j],eqaux,taux,pa); #save("Temp.jld","pr1",pr1)
end

save("ShiftTauIPolicy.jld","prvec",prvec)


################################################################################
################################ Get Results ###################################
################################################################################

#1. TauD should decrease the value of entrants
prvecd=load("ShiftTauDPolicy.jld","prvec");
for j=1:size(prvecd)[2]
  expval=pa.invariant_distr'*prvecd[j].firmvaluegrid[1,:]';
  println(@sprintf("%9.3f \t %9.3f \t %9.3f \t %9.3f \t %9.3f  \t %9.3f \t %9.3f \t %9.3f \t %9.3f \t ExpVal = %9.3f", prvecd[j].firmvaluegrid[1,1],prvecd[j].firmvaluegrid[1,2], prvecd[j].firmvaluegrid[1,3], prvecd[j].firmvaluegrid[1,4], prvecd[j].firmvaluegrid[1,5], prvecd[j].firmvaluegrid[1,6], prvecd[j].firmvaluegrid[1,7], prvecd[j].firmvaluegrid[1,8], prvecd[j].firmvaluegrid[1,9] , expval[1]))
end

#2. TauC should decrease the value of entrants by more
prvecc=load("ShiftTauCPolicy.jld","prvec");
for j=1:size(prvecc)[2]
  println(@sprintf("%9.3f \t %9.3f \t %9.3f \t %9.3f \t %9.3f  \t %9.3f \t %9.3f \t %9.3f \t %9.3f ", prvecc[j].firmvaluegrid[1,1],prvecc[j].firmvaluegrid[1,2], prvecc[j].firmvaluegrid[1,3], prvecc[j].firmvaluegrid[1,4], prvecc[j].firmvaluegrid[1,5], prvecc[j].firmvaluegrid[1,6], prvecc[j].firmvaluegrid[1,7], prvecc[j].firmvaluegrid[1,8], prvecc[j].firmvaluegrid[1,9] ))
end

#Increasing taud decreases value by less than increasing tauc by 5 points
##The expression below is always positive
prvecd[2].firmvaluegrid[1,:]-prvecc[5].firmvaluegrid[1,:]

#Increasing taud decreases value by less than decreasing tauc by 5 points increases value
(prvecd[2].firmvaluegrid[1,:] - prvecd[1].firmvaluegrid[1,:]) + (prvecc[3].firmvaluegrid[1,:] - prvecc[4].firmvaluegrid[1,:])
compute_expvalentry(prvecd[2], pa, eq, tau) + compute_expvalentry(prvecc[3], pa, eq, tau)

# Now increase taud by 10 pp and decrease tauc by only 5pp
(prvecd[3].firmvaluegrid[1,:] - prvecd[1].firmvaluegrid[1,:]) + (prvecc[3].firmvaluegrid[1,:] - prvecc[4].firmvaluegrid[1,:])
compute_expvalentry(prvecd[3], pa, eq, tau) + compute_expvalentry(prvecc[3], pa, eq, tau)


#Convexity
(prvecc[7].firmvaluegrid[1,:] - prvecc[6].firmvaluegrid[1,:]) - (prvecc[6].firmvaluegrid[1,:] - prvecc[5].firmvaluegrid[1,:])
(prvecc[6].firmvaluegrid[1,:] - prvecc[5].firmvaluegrid[1,:]) - (prvecc[5].firmvaluegrid[1,:] - prvecc[4].firmvaluegrid[1,:])
(prvecc[5].firmvaluegrid[1,:] - prvecc[4].firmvaluegrid[1,:]) - (prvecc[4].firmvaluegrid[1,:] - prvecc[3].firmvaluegrid[1,:])
(prvecc[4].firmvaluegrid[1,:] - prvecc[3].firmvaluegrid[1,:]) - (prvecc[3].firmvaluegrid[1,:] - prvecc[2].firmvaluegrid[1,:])
(prvecc[3].firmvaluegrid[1,:] - prvecc[2].firmvaluegrid[1,:]) - (prvecc[2].firmvaluegrid[1,:] - prvecc[1].firmvaluegrid[1,:])


(prvecd[7].firmvaluegrid[1,:] - prvecd[6].firmvaluegrid[1,:]) - (prvecd[6].firmvaluegrid[1,:] - prvecd[5].firmvaluegrid[1,:])
(prvecd[6].firmvaluegrid[1,:] - prvecd[5].firmvaluegrid[1,:]) - (prvecd[5].firmvaluegrid[1,:] - prvecd[4].firmvaluegrid[1,:])
(prvecd[5].firmvaluegrid[1,:] - prvecd[4].firmvaluegrid[1,:]) - (prvecd[4].firmvaluegrid[1,:] - prvecd[3].firmvaluegrid[1,:])
(prvecd[4].firmvaluegrid[1,:] - prvecd[3].firmvaluegrid[1,:]) - (prvecd[3].firmvaluegrid[1,:] - prvecd[2].firmvaluegrid[1,:])
(prvecd[3].firmvaluegrid[1,:] - prvecd[2].firmvaluegrid[1,:]) - (prvecd[2].firmvaluegrid[1,:] - prvecd[1].firmvaluegrid[1,:])


## Distribution
getpolicies!(prvecc[4],eq,tau,pa);
dist = stationarydist(eq.E, prvecc[4], eq, tau, pa);

eqaux = init_equilibirium(eq.w,tau,pa);
eqaux.E=eq.E;
eqaux.distr=dist;

taux=Taxes(tau.d,0.35, tau.i, tau.g, tau.l)

aggregates!(prvecc[4], eqaux, taux, pa);


x=sum(dist,2)/sum(dist)


getpolicies!(prvecc[3],eq,tau,pa);
dist30 = stationarydist(eq.E, prvecc[3], eq, tau, pa);
y=sum(dist30,2)/sum(dist30)


eqaux30 = init_equilibirium(eq.w,tau,pa);
eqaux30.E=eq.E;
eqaux30.distr=dist30;
taux30=Taxes(tau.d,0.3, tau.i, tau.g, tau.l)

aggregates!(prvecc[3], eqaux30, taux30, pa; consistencychecks=false);


#2. Check if by changing tauD and tauG at the same time, firm value is decreased. It is, but by less than it is when taud only is increased
prvecdg=load("ShiftTauDGPolicy.jld","prvec");
for j=1:size(prvecdg)[2]
  expval=pa.invariant_distr'*prvecdg[j].firmvaluegrid[1,:]';
  println(@sprintf("%9.3f \t %9.3f \t %9.3f \t %9.3f \t %9.3f  \t %9.3f \t %9.3f \t %9.3f \t %9.3f \t ExpVal = %9.3f", prvecdg[j].firmvaluegrid[1,1],prvecdg[j].firmvaluegrid[1,2], prvecdg[j].firmvaluegrid[1,3], prvecdg[j].firmvaluegrid[1,4], prvecdg[j].firmvaluegrid[1,5], prvecdg[j].firmvaluegrid[1,6], prvecdg[j].firmvaluegrid[1,7], prvecdg[j].firmvaluegrid[1,8], prvecdg[j].firmvaluegrid[1,9], expval[1] ))
end

#3. Increasing Tau i should increase firm value
prveci=load("ShiftTauIPolicy.jld","prvec");
for j=1:size(prveci)[2]
  expval=pa.invariant_distr'*prveci[j].firmvaluegrid[1,:]';
  println(@sprintf("%9.3f \t %9.3f \t %9.3f \t %9.3f \t %9.3f  \t %9.3f \t %9.3f \t %9.3f \t %9.3f \t ExpVal = %9.3f", prveci[j].firmvaluegrid[1,1],prveci[j].firmvaluegrid[1,2], prveci[j].firmvaluegrid[1,3], prveci[j].firmvaluegrid[1,4], prveci[j].firmvaluegrid[1,5], prveci[j].firmvaluegrid[1,6], prveci[j].firmvaluegrid[1,7], prveci[j].firmvaluegrid[1,8], prveci[j].firmvaluegrid[1,9], expval[1] ))
end
