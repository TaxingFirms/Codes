

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

p35=load("/Users/danielwillsr/Dropbox/1FirmTaxation/SimpleDiscreteTime/ModelResults.jld","p")
p25=load("/Users/danielwillsr/Dropbox/1FirmTaxation/SimpleDiscreteTime/Reform025.jld","p")
p30=load("/Users/danielwillsr/Dropbox/1FirmTaxation/SimpleDiscreteTime/Reform03.jld","p")
p40=load("/Users/danielwillsr/Dropbox/1FirmTaxation/SimpleDiscreteTime/Reform04.jld","p")
p45=load("/Users/danielwillsr/Dropbox/1FirmTaxation/SimpleDiscreteTime/Reform045.jld","p")


prodis25 = sum(p25.distr,1)
prodis30 = sum(p30.distr,1)
prodis35 = sum(p35.distr,1)
prodis40 = sum(p40.distr,1)
prodis45 = sum(p45.distr,1)

normdis25= sum(p25.distr,1)/sum(p25.distr)
normdis30 = sum(p30.distr,1)
normdis35 = sum(p35.distr,1)
normdis40 = sum(p40.distr,1)
normdis45 = sum(p45.distr,1)

save("/Users/danielwillsr/Dropbox/1FirmTaxation/SimpleDiscreteTime/Codes/ProductivityDistributions.jld", "prodis25", prodis25, "prodis30", prodis30, "prodis35", prodis35, "prodis40", prodis40,"prodis45", prodis45,"normdis25", normdis25, "normdis30", normdis30, "normdis35", normdis35, "normdis40", normdis40,"normdis45", normdis45);
