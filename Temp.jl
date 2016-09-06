
pr,eq,tau,pa=load("ModelResults.jld", "pr","eq","tau","pa");

include("TaxMaximization.jl")

tauhat = deepcopy(tau)
tauhat.c=0.3;

govexp=eq.a.G;

@time close_gov_bc(govexp,tauhat,tau, eq, pr, pa; update =0.9, tol = 10^-4.0)






@time close_gov_bc!(govexp,tauhat,tau, eq, pr, pa; update =0.9, tol = 10^-4.0)



govexp =eq.a.G
update =0.5; verbose = true; tol =10.0^-3.0; returnall = false; updateVFIguess= true;
