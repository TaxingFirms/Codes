# All you need to run the code
# Modify ~/.juliarc.jl and add the following line:
# push!(LOAD_PATH, "/Path/To/My/Module/")

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
include("Simulations.jl")
include("Magnitudes.jl")
include("Reforms.jl")

#See August13jl for parameter selection
pa =init_parameters( H=1.176, bbeta=0.972, ff= 0.5, aalphak=0.23, aalphal=0.64, llambda0=0.004, llambda1= 0.04, ddelta = 0.13,
                      allowance=0.86, ttheta = 0.25,rhoz= 0.75, ssigmaz= 0.08, e=0.0, A=1.0);
tau = init_taxes(ttaud =0.15, ttauc= 0.35, ttaui= 0.29, ttaug= 0.15, ttaul=0.28);
@time pr,eq= SolveSteadyState(tau,pa;wguess=0.54, VFItol=10.0^-3.0, verbose=false);
moments=computeMomentsCutoff(eq.E,pr,eq,tau,pa,cutoffCapital=0.0);
save("ModelResults.jld","pr",pr,"eq",eq,"tau",tau,"pa",pa);
#pr,eq,tau,pa=load("ModelResults.jld", "pr","eq","tau","pa");





#taxesb=[0.33 0.31 0.29 0.27 0.25 0.23 0.21 0.19 0.17 0.15 0.13 0.11 0.09 0.07 0.05 0.03 0.01];
taxesb=[0.3 0.25 0.2 0.15 0.1 0.05 0.0];
taxesa=[0.4 0.45 0.5];

Reform2Vector("BenchmarkReforms2.jld", taxesb, pr, eq, tau, pa)









("Parameters: delta ",initialParams[1], " ttheta ",initialParams[2], " rhoz ", initialParams[3], " ssigmaz ", initialParams[4],
  " llambda0 ",initialParams[5], " llambda1 ",initialParams[6] , " f ",initialParams[7], "H", initialParams[8]


#################
# Calibration
#################

# Optimization
#         delta   theta   rhoz    sigmaz   lambda0   lambda1    f       H
LB  = [    .09,      .2,   .55,      .03 ,   .01,       .001     0.001   0.001]
#         delta   theta   rhoz    sigmaz   lambda0   lambda1    f       H
UB  = [    .11,     .4,   .75 ,     .15,    1.0,         .06     5.0    5.0 ]

initialGuess = [0.12,0.3,0.7,0.041,0.02,0.04,0.15,1.3]
count        = 0



using NLopt
using Calculus

function f(x::Vector,grad::Vector)

    g(y) =  try
                computeDistance(y)
            catch eexception
                println(eexception)
                100000000000.0
            end


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
