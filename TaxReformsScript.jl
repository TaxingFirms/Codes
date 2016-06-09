
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

hp = init_hhparameters();
fp  = init_firmparameters(hp);

using JLD
#####################
# TAX REFORMS
pr,tau,fp,res,p=load("/home/dwills/firms/ModelResultsTxExit.jld", "pr","tau","fp","res","p")

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

#Compute the model on first time
@time firmVFIParallel!(pr,p,tau,fp); #pr is updated, computes Value Function
@time w=free_entry!(pr, p, tau, fp,hp)
res=copy_opt_policies(pr);
getpolicies!(res,pr,p,tau,fp);  #r is updated exctracts policies
equilibrium_aggregates!( res, pr, p, tau, fp);

save("/home/dwills/firms/Reform03TxExit.jld", "pr", pr, "tau", tau, "fp", fp, "res",res, "p",p);

##############################################
tau = init_taxes() #0.15, 0.3, 0.3, 0.15);
tauc=0.4
  x= (tauc - tau.c)*C/(D+I+G);
  tau = Taxes(tau.d-x,tauc,tau.i-x,tau.g-x);

p = guess_prices(tau,fp,hp);
pr  = init_firmproblem(p,tau,fp,hp);

#Compute the model on first time
@time firmVFIParallel!(pr,p,tau,fp); #pr is updated, computes Value Function
@time w=free_entry!(pr, p, tau, fp,hp)
res=copy_opt_policies(pr);
getpolicies!(res,pr,p,tau,fp);  #r is updated exctracts policies
equilibrium_aggregates!( res, pr, p, tau, fp);

save("/home/dwills/firms/Reform04TxExit.jld", "pr", pr, "tau", tau, "fp", fp, "res",res, "p",p);

##############################################
tau = init_taxes() #0.15, 0.3, 0.3, 0.15);
tauc=0.45
  x= (tauc - tau.c)*C/(D+I+G);
  tau = Taxes(tau.d-x,tauc,tau.i-x,tau.g-x);

p = guess_prices(tau,fp,hp);
pr  = init_firmproblem(p,tau,fp,hp);

#Compute the model on first time
@time firmVFIParallel!(pr,p,tau,fp); #pr is updated, computes Value Function
@time w=free_entry!(pr, p, tau, fp,hp)
res=copy_opt_policies(pr);
getpolicies!(res,pr,p,tau,fp);  #r is updated exctracts policies
equilibrium_aggregates!( res, pr, p, tau, fp);

save("/home/dwills/firms/Reform045TxExit.jld", "pr", pr, "tau", tau, "fp", fp, "res",res, "p",p);

##############################################
tau = init_taxes() #0.15, 0.3, 0.3, 0.15);
tauc=0.25
  x= (tauc - tau.c)*C/(D+I+G);
  tau = Taxes(tau.d-x,tauc,tau.i-x,tau.g-x);

p = guess_prices(tau,fp,hp);
pr  = init_firmproblem(p,tau,fp,hp);

#Compute the model on first time
@time firmVFIParallel!(pr,p,tau,fp); #pr is updated, computes Value Function
@time w=free_entry!(pr, p, tau, fp,hp)
res=copy_opt_policies(pr);
getpolicies!(res,pr,p,tau,fp);  #r is updated exctracts policies
equilibrium_aggregates!( res, pr, p, tau, fp);

save("/home/dwills/firms/Reform025TxExit.jld", "pr", pr, "tau", tau, "fp", fp, "res",res, "p",p);

