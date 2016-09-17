
using ClusterManagers
ClusterManagers.addprocs_sge(59,queue="all.q",qsub_env="LD_LIBRARY_PATH");

@everywhere using Grid:CoordInterpGrid, BCnan, BCnearest, InterpLinear
@everywhere using Roots:fzero
using DataFrames

include("markov_approx.jl")
include("mc_tools.jl")

@everywhere include("Main.jl")
@everywhere include("Firms.jl")
@everywhere include("FreeEntry.jl")
@everywhere include("Distribution.jl")
@everywhere include("Aggregation.jl")
@everywhere include("SolveSteadyState.jl")
@everywhere include("TaxReforms.jl")
@everywhere include("calibrate.jl")
@everywhere include("Transitions.jl")


#################
# Calibration
#################

# Optimization
#         delta    theta    rhoz    sigmaz    lambda0    lambda1     f       e
LB  = [    0.0268,    0.074,   0.55,   0.056,    0.01,    0.0507,    0.0,     0.0 ]
#         delta    theta    rhoz    sigmaz    lambda0    lambda1     f      e
UB  = [    0.1268,    0.474,   0.95,   0.116,   0.05,     0.2507,      1.592,     0.0776]

initialGuess = [0.0768,0.274,0.75, 0.086,0.03, 0.1507, 0.796,0.038]
count        = 0

using NLopt
using Calculus

function f(x::Vector,grad::Vector)
    g(y) =  try
                computeDistance(y)
            catch eexception
                if isa(eexception,ErrorException)
                    println(eexception)
                    100000000000.0
                else
                    println(eexception)
                    throw(eexception)
                end
            end
	if length(grad) > 0
		grad[:] = Calculus.gradient(g,x)
	end
	answer = g(x)
	global count
    count::Int += 1
    println(calout, "it ",count,"  [",x, "] = ",answer)
    mod1(count,50)==1 && flush(calout)
    answer[1]
end


# Simulated Annealing First

opt = Opt(:GN_DIRECT_L,length(LB))

lower_bounds!(opt,LB)
upper_bounds!(opt,UB)
min_objective!(opt,f)
xtol_rel!(opt,.1)

calout=open("sobolcalout.txt","a")
println(calout, "-------------------------------------------------------------------------------------------------------")
flush(calout)

(minf,minx,ret) = optimize(opt,initialGuess)


println(calout, "=======================================================================================================")
close(calout)
println(calout,"got $minf at $minx after $count iterations (returned $ret)")


# Nelder-Mead Locally to improve

# opt = Opt(:LN_SBPLX,length(LB))
# lower_bounds!(opt,LB)
# upper_bounds!(opt,UB)
# min_objective!(opt,f)
# xtol_rel!(opt,.1)

# (minf,minx,ret) = optimize(opt,initialGuess)
# println("got $minf at $minx after $count iterations (returned $ret)")
