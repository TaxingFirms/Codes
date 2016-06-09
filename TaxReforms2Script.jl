
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
pr,tau,fp,res,p=load("/home/dwills/firms/ModelResults.jld", "pr","tau","fp","res","p")

##############################################

tau = init_taxes(ttaud=0.2,ttaui=0.25);
p = guess_prices(tau,fp,hp);
pr  = init_firmproblem(p,tau,fp,hp);

#Compute the model on first time
@time firmVFIParallel!(pr,p,tau,fp); #pr is updated, computes Value Function
@time w=free_entry!(pr, p, tau, fp,hp) #0.6895019531249997
res=copy_opt_policies(pr);
getpolicies!(res,pr,p,tau,fp);  #r is updated exctracts policies
equilibrium_aggregates!( res, pr, p, tau, fp);

save("/home/dwills/firms/DividendTx20.jld", "pr", pr, "tau", tau, "fp", fp, "res",res, "p",p);

##############################################
