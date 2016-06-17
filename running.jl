# All you need to run the code
# Modify ~/.juliarc.jl and add the following line:
# push!(LOAD_PATH, "/Path/To/My/Module/")
# push!(LOAD_PATH, "/home/dwills/firms/Codes/")
# push!(LOAD_PATH,"/Users/danielwillsr/Dropbox/1FirmTaxation/SimpleDiscreteTime/Codes")
# /data/global/bin/julia

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
@everywhere include("Temp.jl")

hp = init_hhparameters();
fp  = init_firmparameters(hp);

tau = init_taxes() #0.15, 0.3, 0.3, 0.15);
p = guess_prices(tau,fp,hp);
pr  = init_firmproblem(p,tau,fp,hp);



#Compute the model on first time
@time firmVFIParallelOmega!(pr,p,tau,fp); #pr is updated, computes Value Function
#597.841274 seconds on tesla, tol = 10^-3.

#Compute wage such that free entry condition holds
@time w=free_entry!(pr, p, tau, fp,hp,tol=.001)
#1339.480311 seconds on tesla, tol = 10^-2, w = 0.719

#Extract policies and other idiosyncratic results of interest
res=copy_opt_policies(pr);
getpolicies!(res,pr,p,tau,fp);  #r is updated exctracts policies

#Compute mass of entrants and stationary distribution
# both are updated in p.
mass_of_entrants!( res, pr, p, tau, fp);

# Compute aggregate results of interest and moments
aggregates!(res, pr, p, tau, hp, fp);
p.E/sum(p.distr)
using JLD
#save("/home/gcam/firms/Codes/FreeEntryResults.jld", "pr", pr, "tau", tau, "fp", fp, "res",res,"p",p);
save("/home/dwills/firms/ModelResults.jld", "pr", pr, "tau", tau, "fp", fp, "res",res, "p",p);



#####################
# PLOTS

using UnicodePlots
a = lineplot(collect(pr.omega.grid),squeeze(mean(res.distributions,2),2),title="Average Policies", name = "dividends");
lineplot!(a,collect(pr.omega.grid),squeeze(mean(res.kprime,2),2),name="K Policy");
lineplot!(a,collect(pr.omega.grid),squeeze(mean(res.qprime,2),2),name="Q Policy");
a

b = lineplot(collect(pr.omega.grid),pr.kpolicygrid[:,1],title="K Policy", name = "z=1");
map(x->lineplot!(b,collect(pr.omega.grid),pr.kpolicygrid[:,x],name=string("z=",x)),2:pr.Nz);
b

c = lineplot(collect(pr.omega.grid),pr.qpolicygrid[:,1],title="Q Policy", name = "z=1");
map(x->lineplot!(c,collect(pr.omega.grid),pr.qpolicygrid[:,x],name=string("z=",x)),2:pr.Nz);
c

initPlot = 5 # Sometimes 0s in policy function
d = lineplot(collect(pr.omega.grid)[initPlot:pr.Nomega],pr.qpolicygrid[initPlot:pr.Nomega,1]./pr.kpolicygrid[initPlot:pr.Nomega,1],title="Ratio Q/K Policy", name = "z=1");
map(x->lineplot!(d,collect(pr.omega.grid)[initPlot:pr.Nomega],pr.qpolicygrid[initPlot:pr.Nomega,x]./pr.kpolicygrid[initPlot:pr.Nomega,x],name=string("z=",x)),2:pr.Nz);
d

e = lineplot(collect(pr.omega.grid),res.exitprobability[:,1],title="Exit Probability", name = "z=1");
map(x->lineplot!(e,collect(pr.omega.grid),res.exitprobability[:,x],name=string("z=",x)),2:pr.Nz);
e


f = lineplot(collect(pr.omega.grid),distr[:,1],title="Invariant Distribution", name = "z=1");
map(x->lineplot!(f,collect(pr.omega.grid),distr[:,x],name=string("z=",x)),2:pr.Nz);
f

g = lineplot(collect(pr.omega.grid),res.distributions[:,1],title="Net Distributions", name = "z=1");
map(x->lineplot!(g,collect(pr.omega.grid),dist[:,x],name=string("z=",x)),2:pr.Nz);
g
