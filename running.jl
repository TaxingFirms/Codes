# All you need to run the code
# Modify ~/.juliarc.jl and add the following line:
# push!(LOAD_PATH, "/Path/To/My/Module/")


@everywhere cd("C:/Users/m1dsw00.BOARD/Codes")

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


pa  = init_parameters();
tau = init_taxes();
@time pr,eq= SolveSteadyState(tau,pa);
save("ModelResults.jld","pr",pr,"eq",eq,"tau",tau,"pa",pa);
#pr,eq,tau,pa=load("ModelResults.jld", "pr","eq","tau","pa");

@time pr2,eq2,tau2 = taxreform2(0.3, eq, tau, pa;update=0.0);
save("Counterfactual2.jld","pr",pr2,"eq",eq2,"tau",tau2,"pa",pa);
#pr2,eq2,tau2,pa=load("Counterfactual2.jld", "pr","eq","tau","pa");

@time pr3,eq3,tau3 = taxreform3(0.3, eq, tau, pa; tol=10.0^-4.0, update = 0.9);
save("Counterfactual3.jld","pr",pr3,"eq",eq3,"tau",tau3,"pa",pa);

@time pr4,eq4,tau4 = taxreform3(0.0, eq, tau, pa);
save("Counterfactual4.jld","pr",pr4,"eq",eq4,"tau",tau4,"pa",pa);

#NOT WORKING
##@time pr1,eq1,tau1 = taxreform1(0.3, eq, tau, pa);
##save("Counterfactual1.jld","pr",pr1,"eq",eq1,"tau",tau1,"pa",pa);
##


#################
# Transitions
#################




#################
# Calibration
#################

# Optimization
#         delta     rhoz    sigmaz   theta   lambda0   lambda1
LB  = [     .01,      .5,    .01,    .01 ,   .01,       .0001 ]
#         delta     rhoz    sigmaz   theta
UB  = [     .15,     .95,   .50 ,     .8,    .15,         .03  ]

initialGuess = [0.14,0.76,0.0352,.45,.08,.028]
count        = 0



using NLopt
using Calculus

function f(x::Vector,grad::Vector)
	g(y) = computeDistance(x)

	if length(grad) > 0
		grad[:] = Calculus.gradient(g,x)
	end
	answer = g(x)
	global count
    count::Int += 1
    println("f_$count($x)=$answer")
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
