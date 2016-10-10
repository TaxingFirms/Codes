@everywhere using Grid:CoordInterpGrid, BCnan, BCnearest, InterpLinear
using JLD
using DataFrames
using PyPlot

include("markov_approx.jl")
include("mc_tools.jl")

@everywhere include("Main.jl")
@everywhere include("Firms.jl")
@everywhere include("FreeEntry.jl")
@everywhere include("Distribution.jl")
@everywhere include("Aggregation.jl")
@everywhere include("SolveSteadyState.jl")
@everywhere include("calibrate.jl")
@everywhere include("Transitions.jl")

# ModelResults500.jld is built and saved in runShiftPolicies.jl
pr,eq,tau,pa =load("ModelResults500.jld","pr","eq","tau","pa");

pctdividendpayer = Array(Float64,(pa.Nz,));
outputdividendpayer = Array(Float64,(pa.Nz,));

pctequityissuer = Array(Float64,(pa.Nz,));
outputequityissuer = Array(Float64,(pa.Nz,));
