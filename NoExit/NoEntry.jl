#Model witout entry or exit
using QuantEcon:tauchen
using Grid:CoordInterpGrid, BCnan, BCnearest, InterpLinear
using Roots:fzero
using JLD
using PyPlot

include("Main.jl")
include("Firms.jl")
include("Policies.jl")
include("Distribution.jl")
include("FreeEntry.jl")
include("Aggregation.jl")
include("TaxReforms.jl")


function equilibrium_aggregatesNOENTRY!( res::ResultsFP, pr::FirmProblem, p::Equilibrium, tau::Taxes, fp::FirmParam)
  #Computes the mass of entrants such that the labor market clears,
  distr1=stationarydist(0,res, pr, p, tau,fp);
  capital1, bonds1, labor_d1, gdp1, corptax1, inctax1, netdistributions1, liquidations1, liquidationcosts1 = aggregates(0, distr1, res, pr, p, tau, fp);
  investment1 = sum(distr1.*res.kprime) - (1-fp.delta)*capital1;
  grossdividends1=sum(distr1.*res.grossdividends);
  divtax1= tau.d*grossdividends1;
  financialcosts1= - sum(distr1.*res.financialcosts);
  G1 = divtax1 + corptax1 + inctax1;
end

tau = init_taxes() #0.15, 0.3, 0.3, 0.15);
hp = init_hhparameters();
fp  = init_firmparameters(hp,ff=0);
p = guess_prices(tau,fp,hp);
pr  = init_firmproblem(p,tau,fp,hp);


@time firmVFIParallel!(pr,p,tau,fp;mp=true); #pr is updated, computes Value Function

res=copy_opt_policies(pr);
getpolicies!(res,pr,p,tau,fp);  #r is updated exctracts policies

equilibrium_aggregatesNOENTRY!( res, pr, p, tau, fp);






