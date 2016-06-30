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
@everywhere include("FreeEntry.jl")
@everywhere include("Distribution.jl")
@everywhere include("Aggregation.jl")
@everywhere include("SolveModel.jl")
@everywhere include("TaxReforms.jl")


pa  = init_parameters();
tau = init_taxes();

@time pr,eq= SolveModel!(tau,pa);
save("ModelResults.jld","pr",pr,"eq",eq,"tau",tau,"pa",pa);

# Speed Benchmark: 41 seconds, 41 M, 3.7GB
## @time pr,eq= SolveModel!(tau,pa; maxroutine=maximizationbf);
# 650 s, 42 M aaloc 3.7 gb
#
## pa  = init_parameters( Nz=9, Nk=150, Nq=75, Nomega=150 );
## @time pr,eq= SolveModel!(tau,pa);
# 117 seconds
#
## pa  = init_parameters( Nz=15 );
## @time pr,eq= SolveModel!(tau,pa);
# 113 second



save("ModelResults.jld","pr",pr,"eq",eq,"tau",tau,"pa",pa);


pr1,eq1,tau1 = taxreform1(0.3, eq, tau, pa);
save("Counterfactual1.jld","pr",pr1,"eq",eq1,"tau",tau1,"pa",pa);

pr2,eq2,tau2 = taxreform2(0.3, eq, tau, pa);
save("Counterfactual2.jld","pr",pr2,"eq",eq2,"tau",tau2,"pa",pa);
