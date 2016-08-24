
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
#         delta    theta    rhoz    sigmaz    lambda0    lambda1     f     H
LB  = [    0.03,    0.05,   0.55,     0.01,       0.0,      0.0,   0.0,   0.0]
#         delta    theta    rhoz    sigmaz    lambda0    lambda1     f     H
UB  = [    0.23,    0.45,   0.95,     0.15,       0.5,     0.08,   1.0,   2.4]


initialGuess = [0.13,0.25,0.75,0.08,0.004, 0.04, 0.5,1.176]
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
    println("f",count,"  [",x, "] = ",answer)
    answer
end


# Simulated Annealing First

opt = Opt(:GN_DIRECT_L,length(LB))

lower_bounds!(opt,LB)
upper_bounds!(opt,UB)
min_objective!(opt,f)
xtol_rel!(opt,.1)

(minf,minx,ret) = optimize(opt,initialGuess)
println("got $minf at $minx after $count iterations (returned $ret)")


# Nelder-Mead Locally to improve

# opt = Opt(:LN_SBPLX,length(LB))
# lower_bounds!(opt,LB)
# upper_bounds!(opt,UB)
# min_objective!(opt,f)
# xtol_rel!(opt,.1)

# (minf,minx,ret) = optimize(opt,initialGuess)
# println("got $minf at $minx after $count iterations (returned $ret)")
