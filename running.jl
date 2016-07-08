# All you need to run the code
# Modify ~/.juliarc.jl and add the following line:
# push!(LOAD_PATH, "/Path/To/My/Module/")
# @everywhere cd("C:/Users/m1dsw00.BOARD/Codes")

@everywhere using QuantEcon:tauchen
@everywhere using Grid:CoordInterpGrid, BCnan, BCnearest, InterpLinear
@everywhere using Roots:fzero
@everywhere using JLD
@everywhere using Dierckx:Spline1D
@everywhere include("Main.jl")
@everywhere include("Firms.jl")
@everywhere include("FreeEntry.jl")
@everywhere include("Distribution.jl")
@everywhere include("Aggregation.jl")
@everywhere include("SolveSteadyState.jl")
@everywhere include("TaxReforms.jl")
@everywhere include("calibrate.jl")


pa  = init_parameters();
tau = init_taxes(ttaud=0.15, ttauc=0.3, ttaui = 0.35, ttaug = 0.0);

@time pr,eq= SolveSteadyState!(tau,pa;maxroutine=maximizationbf);
save("ModelResults.jld","pr",pr,"eq",eq,"tau",tau,"pa",pa);

# Speed Benchmark: 41 seconds, 41 M, 3.7GB
## @time pr,eq= SolveModel!(tau,pa; maxroutine=maximizationbf);
# 650 s, 42 M aaloc 3.7 gb
#
## pa  = init_parameters( Nz=9, Nk=150, Nq=75, Nomega=150 );
## @time pr,eq= SolveModel!(tau,pa);
# 117 seconds
#
## pa  = init_parameters( Nz=15 );
## @time pr,eq= SolveModel!(tau,pa);
# 113 second



save("ModelResults.jld","pr",pr,"eq",eq,"tau",tau,"pa",pa);


pr1,eq1,tau1 = taxreform1(0.3, eq, tau, pa);
save("Counterfactual1.jld","pr",pr1,"eq",eq1,"tau",tau1,"pa",pa);

#tax reform 2 is not converging after the
pr2,eq2,tau2 = taxreform2(0.3, eq, tau, pa);
save("Counterfactual2.jld","pr",pr2,"eq",eq2,"tau",tau2,"pa",pa);


#################
# Calibration
#################

# Optimization
#         delta     rhoz    sigmaz   theta   lambda0   lambda1
LB  = [     .01,      .5,    .01,    .01 ,   .01,       .0001 ]
#         delta     rhoz    sigmaz   theta
UB  = [     .15,     .95,   .50 ,     .8,    .15,         .03  ]

initialGuess = [0.14,0.76,0.0352,.45,.08,.028]
count = 0



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
