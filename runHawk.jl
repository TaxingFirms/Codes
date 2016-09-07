using ClusterManagers
ClusterManagers.addprocs_sge(36,queue="openmp.q",qsub_env="LD_LIBRARY_PATH");

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
#         delta    theta    rhoz    sigmaz    lambda0    lambda1     f     H    e
LB  = [    0.075,    0.1,   0.55,     0.052,   0.001,      0.05,   0.0,   0.57,  0.0 ]
#         delta    theta    rhoz    sigmaz    lambda0    lambda1     f     H    e
UB  = [    0.175,    0.3,   0.95,     0.092,   0.007,      0.45,   1.0,   1.57,  0.3]

initialGuess = [0.125,0.2,0.75,0.072,0.004, 0.25, 0.5,1.07,0.15]
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
    mod1(count,50)==50 && flush(calout)
    answer
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
println("got $minf at $minx after $count iterations (returned $ret)")


# Nelder-Mead Locally to improve

# opt = Opt(:LN_SBPLX,length(LB))
# lower_bounds!(opt,LB)
# upper_bounds!(opt,UB)
# min_objective!(opt,f)
# xtol_rel!(opt,.1)

# (minf,minx,ret) = optimize(opt,initialGuess)
# println("got $minf at $minx after $count iterations (returned $ret)")
