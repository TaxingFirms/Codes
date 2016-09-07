
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
LB  = [    0.0268,    0.074,   0.55,   0.056,   0.0507,      0.01,   0.0,     0.0 ]
#         delta    theta    rhoz    sigmaz    lambda0    lambda1     f      e
UB  = [    0.1268,    0.474,   0.95,   0.116,   0.2507,      0.05,   1.592,     0.0776]

initialGuess = [0.0768,0.274,0.75, 0.086, 0.1507,0.03, 0.796,0.038]
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

@everywhere using Grid:CoordInterpGrid, BCnan, BCnearest, InterpLinear
@everywhere using Roots:fzero
using QuantEcon:tauchen
using JLD
using DataFrames

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
LB  = [    0.0,      0.1,    .01,     0.0,       0.0,     .0001,   0.0,  0.01]
#         delta    theta    rhoz    sigmaz    lambda0    lambda1     f     H
UB  = [    0.2,     0.95,    .95,     0.5,       0.5,       .09,   1.0,  10.0]


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
