# Check Reform

@everywhere using Grid:CoordInterpGrid, BCnan, BCnearest, InterpLinear
using JLD
using DataFrames

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


reformFiles = readdir()[map(x-> contains(x,"ZeroTauI.jld"),readdir())]

for file in reformFiles

    ref,pa = load(file,"ref","pa");

    nReforms = length(ref)

    saveReform = Array(Float64,nReforms,16)

    for index = 1:nReforms

        value = ref[index]
        # Save the aggregate variables    

        saveReform[index,1]  = value.eq.a.consumption
        saveReform[index,2]  = value.eq.a.output
        saveReform[index,3]  = value.eq.a.netdistributions
        saveReform[index,4]  = value.eq.a.agginterests
        saveReform[index,5]  = value.eq.a.laborsupply
        saveReform[index,6]  = value.eq.a.welfare
        saveReform[index,7]  = value.eq.a.bonds
        saveReform[index,8]  = value.eq.a.capital
        saveReform[index,9]  = value.eq.a.investment
        saveReform[index,10] = value.eq.a.grossdividends
        saveReform[index,11] = value.eq.a.G
        saveReform[index,12] = value.tau.c
        saveReform[index,13] = value.tau.d
        saveReform[index,14] = value.tau.i
        saveReform[index,15] = value.tau.g
        saveReform[index,16] = value.tau.l

    end

    varNames = ["consumption","output","netdistributions","agginterests","laborsupply","welfare","bonds","capital",
    "investment","grossdividends","govExp","corpProfitTax","dividendTax","interestTax","capGainsTax","laborTax"]

    reforms = DataFrame(saveReform)
    names!(reforms,map(symbol,varNames))
    writetable(string(file[1:(length(file)-3)],"csv"),reforms)

end
