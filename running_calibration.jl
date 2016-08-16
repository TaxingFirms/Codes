
@everywhere using Grid:CoordInterpGrid, BCnan, BCnearest, InterpLinear
@everywhere using Roots:fzero
using QuantEcon:tauchen
using JLD
using DataFrames
using StatsFuns

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
#         delta     rhoz    sigmaz   theta   lambda0   lambda1
LB  = [     .01,      .5,    .01,    .01 ,   .01,       .0001, .005]
#         delta     rhoz    sigmaz   theta
UB  = [     .15,     .95,   .50 ,     .8,    .15,         .03, .06]

initialGuess = [0.14,0.76,0.0352,.45,.08,.028,.0145]
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
