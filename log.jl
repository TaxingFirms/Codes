@everywhere using QuantEcon:tauchen
@everywhere using Grid:CoordInterpGrid, BCnan, BCnearest, InterpLinear
@everywhere using Roots:fzero
@everywhere include("Main.jl")
@everywhere include("Firms.jl")
@everywhere include("Policies.jl")
@everywhere include("Distribution.jl")
@everywhere include("FreeEntry.jl")
@everywhere include("Aggregation.jl")
@everywhere include("TaxReforms.jl")
@everywhere include("FirmsConstantExit.jl")
@everywhere include("labormarketNE.jl")


hp = init_hhparameters();
fp  = init_firmparameters(hp);

using JLD
#####################
# TAX REFORMS
 pr,tau,fp,res,p = load("/home/dwills/firms/ModelResults.jld", "pr","tau","fp","res","p");
# pr,tau,fp,res,p = load("/home/gcam/firms/ModelResults.jld", "pr","tau","fp","res","p");
##############################################

E = p.E
totalfirms=sum(p.distr);
exitprob = sum(res.exitprobability.*p.distr)/sum(p.distr);


#############################
######## BENCHMARK ##########
#############################

tau = init_taxes();
p   = guess_prices(tau,fp,hp);
pr  = init_firmproblem(p,tau,fp,hp);

#Compute the model on first time
@time firmVFIParallelOmegaCE!(exitprob,pr,p,tau,fp); #pr is updated, computes Value Function
res=copy_opt_policies(pr);
getpolicies!(res,pr,p,tau,fp);  #r is updated exctracts policies

p.E = E
distr = stationarydistCE(p.E, pr, p, tau, fp)
#distr = distributionStupid(p.E,res, pr, p, tau, fp)

clear_labormarketCE!(p.E,exitprob,res,pr,p,tau,fp,tol=.01)

save("/home/dwills/firms/BenchmarkNOEXIT.jld", "pr", pr, "tau", tau, "fp", fp, "res",res, "p",p);


# for tauc in [0.25 0.3 0.4 0.45]
  C = p.collections.c / tau.c #corporate base
  D = p.collections.d / tau.d #corporate base
  I = p.collections.i / tau.i #corporate base
  G = p.collections.g / tau.g #corporate base


##############################################
tau = init_taxes() #0.15, 0.3, 0.3, 0.15);
tauc=0.3
  x= (tauc - tau.c)*C/(D+I+G);
  tau = Taxes(tau.d-x,tauc,tau.i-x,tau.g-x);

p = guess_prices(tau,fp,hp);
pr  = init_firmproblem(p,tau,fp,hp);
#Compute the model
@time firmVFIParallelOmegaCE!(exitprob,pr,p,tau,fp); #pr is updated, computes Value Function
res=copy_opt_policies(pr);
getpolicies!(res,pr,p,tau,fp);  #r is updated exctracts policies

p.E = E
distr = stationarydistCE(p.E, pr, p, tau, fp)
#distr = distributionStupid(p.E,res, pr, p, tau, fp)

clear_labormarketCE!(p.E,exitprob,res,pr,p,tau,fp,tol=.01)

save("/home/dwills/firms/NoExit30.jld", "pr", pr, "tau", tau, "fp", fp, "res",res, "p",p);

##############################################
tau = init_taxes() #0.15, 0.3, 0.3, 0.15);
tauc=0.4
  x= (tauc - tau.c)*C/(D+I+G);
  tau = Taxes(tau.d-x,tauc,tau.i-x,tau.g-x);

p = guess_prices(tau,fp,hp);
pr  = init_firmproblem(p,tau,fp,hp);
#Compute the model
@time firmVFIParallelOmegaCE!(exitprob,pr,p,tau,fp); #pr is updated, computes Value Function
res=copy_opt_policies(pr);
getpolicies!(res,pr,p,tau,fp);  #r is updated exctracts policies

p.E = E
distr = stationarydistCE(p.E, pr, p, tau, fp)
#distr = distributionStupid(p.E,res, pr, p, tau, fp)

clear_labormarketCE!(p.E,exitprob,res,pr,p,tau,fp,tol=.01)

save("/home/dwills/firms/NoExit40.jld", "pr", pr, "tau", tau, "fp", fp, "res",res, "p",p);

##############################################
tau = init_taxes() #0.15, 0.3, 0.3, 0.15);
tauc=0.45
  x= (tauc - tau.c)*C/(D+I+G);
  tau = Taxes(tau.d-x,tauc,tau.i-x,tau.g-x);

p = guess_prices(tau,fp,hp);
pr  = init_firmproblem(p,tau,fp,hp);
#Compute the model
@time firmVFIParallelOmegaCE!(exitprob,pr,p,tau,fp); #pr is updated, computes Value Function
res=copy_opt_policies(pr);
getpolicies!(res,pr,p,tau,fp);  #r is updated exctracts policies

p.E = E
distr = stationarydistCE(p.E, pr, p, tau, fp)
#distr = distributionStupid(p.E,res, pr, p, tau, fp)

clear_labormarketCE!(p.E,exitprob,res,pr,p,tau,fp,tol=.01)

save("/home/dwills/firms/NoExit45.jld", "pr", pr, "tau", tau, "fp", fp, "res",res, "p",p);

##############################################

tau = init_taxes() #0.15, 0.3, 0.3, 0.15);
tauc=0.25
  x= (tauc - tau.c)*C/(D+I+G);
  tau = Taxes(tau.d-x,tauc,tau.i-x,tau.g-x);

p = guess_prices(tau,fp,hp);
pr  = init_firmproblem(p,tau,fp,hp);
#Compute the model
@time firmVFIParallelOmegaCE!(exitprob,pr,p,tau,fp); #pr is updated, computes Value Function
res=copy_opt_policies(pr);
getpolicies!(res,pr,p,tau,fp);  #r is updated exctracts policies

p.E = E
distr = stationarydistCE(p.E, pr, p, tau, fp)
#distr = distributionStupid(p.E,res, pr, p, tau, fp)

clear_labormarketCE!(p.E,exitprob,res,pr,p,tau,fp,tol=.01)

save("/home/dwills/firms/NoExit25.jld", "pr", pr, "tau", tau, "fp", fp, "res",res, "p",p);


