@everywhere using QuantEcon:tauchen
@everywhere using Grid:CoordInterpGrid, BCnan, BCnearest, InterpLinear
@everywhere using Roots:fzero
@everywhere include("Main.jl")
@everywhere include("FirmsNOEXIT.jl")
@everywhere include("Policies.jl")
@everywhere include("Distribution.jl")
@everywhere include("Aggregation.jl")
@everywhere include("labormarketNE.jl")



hp = init_hhparameters();
fp  = init_firmparameters(hp);

using JLD
 pr,tau,fp,res,p=load("/home/dwills/firms/ModelResults1.jld", "pr","tau","fp","res","p");
#pr,tau,fp,res,p=load("/home/gcam/firms/ModelResults1.jld", "pr","tau","fp","res","p");
E = copy(p.E);
exitrule=copy(res.exitrule);


#Make at the initial equilibirium, markets are clearing
#clear_labormarket!(E, exitrule, res, pr, p, tau, fp)


###################################################
################## TAX REFORM #####################
###################################################

tau = init_taxes(ttaud=0.2,ttaui=0.25);
p = guess_prices(tau,fp,hp);
pr  = init_firmproblem(p,tau,fp,hp);


@time firmVFIParallelOmega!(pr,p,tau,fp); #pr is updated, computes Value Function
res=copy_opt_policies(pr);
getpolicies!(res,pr,p,tau,fp);  #r is updated exctracts policies

clear_labormarket!(E, exitrule, res, pr, p, tau, fp)
res=copy_opt_policies(pr);
getpolicies!(res,pr,p,tau,fp);  #r is updated exctracts policies


save("/home/dwills/firms/DividendTx20NOEXIT.jld", "pr", pr, "tau", tau, "fp", fp, "res",res, "p",p);

##############################################

