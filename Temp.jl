
pr,eq,tau,pa=load("ModelResults.jld", "pr","eq","tau","pa");

include("TaxMaximization.jl")
include("TaxMaximizationAlpha.jl")

#maximize_welfare()




tauhat = deepcopy(tau)
tauhat.c=0.3;
govexp=eq.a.G;

#update =0.5; verbose = true; tol =10.0^-3.0; returnall = false; updateVFIguess= true;
####################
include("ComparativeStatics.jl")
pr,eq,tau,pa=load("ModelResults.jld", "pr","eq","tau","pa");
@time comparativestatics( eq, pr, tau, pa; tauvec=[0.15, 0.20, 0.25, 0.30, 0.35, 0.40])
