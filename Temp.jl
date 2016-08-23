#cd("C:\\Users\\IF User\\Codes")


using Grid:CoordInterpGrid, BCnan, BCnearest, InterpLinear
using Roots:fzero
using QuantEcon:tauchen
using JLD
using DataFrames
using StatsFuns

include("Main.jl")
include("Firms.jl")
include("FreeEntry.jl")
include("Distribution.jl")
include("Aggregation.jl")
include("SolveSteadyState.jl")
include("TaxReforms.jl")
include("calibrate.jl")
include("Transitions.jl")
include("Simulations.jl")
include("Magnitudes.jl")
include("Reforms.jl")


pr,eq,tau,pa=load("ModelResults.jld", "pr","eq","tau","pa");

#taxesb=[0.33 0.31 0.29 0.27 0.25 0.23 0.21 0.19 0.17 0.15 0.13 0.11 0.09 0.07 0.05 0.03 0.01];
taxesb=[0.3 0.25 0.2 0.15 0.1 0.05 0.0];
taxesa=[0.4 0.45 0.5];


rpr,req,rtau = taxreform6(0.0, eq, tau, pa; tol=10.0^-3.0,update=0.5,momentsprint=true);
  save("Reforms6","rpr",rpr,"req",req,"rtau",rtau,"pa",pa);

#On TEsla
rpr,req,rtau = taxreform5(0.0, eq, tau, pa; tol=10.0^-3.0,update=0.7,momentsprint=true);
save("Reforms5","rpr",rpr,"req",req,"rtau",rtau,"pa",pa);


#ref=Reform5Vector("BenchmarkReforms5b.jld", taxesb, pr, eq, tau, pa)

#ref=Reform2Vector("BenchmarkReforms2b.jld", taxesb, pr, eq, tau, pa)
#ref=Reform3Vector("BenchmarkReforms3b.jld", taxesb, pr, eq, tau, pa)
#ref=Reform1Vector("BenchmarkReforms1b.jld", taxesb, pr, eq, tau, pa)

#ref=Reform2Vector("BenchmarkReforms2a.jld", taxesa, pr, eq, tau, pa)
#ref=Reform3Vector("Benchmark1Reforms3a.jld", taxesa, pr, eq, tau, pa)
#ref=Reform1Vector("Benchmark1Reforms1a.jld", taxesa, pr, eq, tau, pa)
