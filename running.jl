# All you need to run the code
# Modify ~/.juliarc.jl and add the following line:
# push!(LOAD_PATH, "/Path/To/My/Module/")
# push!(LOAD_PATH, "/home/dwills/firms/Codes/")
# push!(LOAD_PATH,"/Users/danielwillsr/Dropbox/1FirmTaxation/SimpleDiscreteTime/Codes")
# /data/global/bin/julia

@everywhere using QuantEcon:tauchen
@everywhere using Grid:CoordInterpGrid, BCnan, BCnearest, InterpLinear
@everywhere using Roots:fzero
@everywhere using JLD
@everywhere include("Main.jl")
@everywhere include("Firms.jl")
@everywhere include("Policies.jl")
@everywhere include("Distribution.jl")
@everywhere include("FreeEntry.jl")
@everywhere include("Aggregation.jl")
@everywhere include("TaxReforms.jl")
@everywhere include("SolveModel.jl")

hp = init_hhparameters();
fp  = init_firmparameters(hp; llambda0=0.02, llambda1=0.04, ttheta=0.42, ff= 0.015, ddelta=0.12);
tau = init_taxes() #0.15, 0.3, 0.3, 0.15);


#p,res,pr= SolveModel!(tau,fp,hp)

#save("/home/dwills/firms/ModelResults.jld", "pr", pr, "tau", tau, "fp", fp, "res",res, "p",p);

pr,tau,fp,res,p=load("/home/dwills/firms/ModelResults.jld", "pr","tau","fp","res","p");


taxreform2(0.3, p, tau, fp, hp)

taxreform1(0.3, p, tau, fp, hp)
