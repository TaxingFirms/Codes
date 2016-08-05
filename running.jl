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

#pa  = init_parameters( bbeta=0.98, ssigma=1.0,psi=0.55, H=3.47, aalphak=0.3,
 #aalphal = 0.65, ff=0.014, llambda0= 0.02, llambda1= 0.04, ddelta= 0.012, ttheta=0.42,
 #kappa=1.0, e=0.0, k0=0.00, rhoz= 0.76, ssigmaz= 0.1, Nz=11, Nk=80, Nq=40, Nomega=100, A = 0.65);
#tau = init_taxes(ttaud =0.12, ttauc= 0.35, ttaui= 0.29, ttaug= 0.12, ttaul=0.28);

# Calibration Board presentation
pa  = init_parameters( H=1.3, ff= 0.15, llambda0=0.02, llambda1= 0.04, ddelta = 0.12, ttheta = 0.45);
tau = init_taxes();

capital_unc, capital_equity, profits_unc, profits_equity =magnitudes(tau, pa);

@time pr,eq= SolveSteadyState(tau,pa; wguess=0.55);
moments=computeMomentsCutoff(eq.E,pr,eq,tau,pa,cutoffCapital=0.0);
capital, debt, networth, dividends, investment, z_history_ind = simulation(10, 50,pr,pa; seed=1234);
figure()
plot(debt)
capital_unc, capital_equity, profits_unc, profits_equity = magnitudes(tau, pa; r= eq.r, w=eq.w);

for i=1:pa.Nz
figure()
d= plot(pa.omega.grid, pr.distributions[:,i] )
k= plot(pa.omega.grid, pr.kpolicy[:,i] )
q= plot(pa.omega.grid, pr.qpolicy[:,i] )
  xlabel("Net worth")
  title("Policy functions")
  legend("dkq", loc="best")
end

save("ModelResults.jld","pr",pr,"eq",eq,"tau",tau,"pa",pa);
#pr,eq,tau,pa=load("ModelResults.jld", "pr","eq","tau","pa");

@time pr2,eq2,tau2 = taxreform2(0.3, eq, tau, pa;update=0.0);
save("Counterfactual2.jld","pr",pr2,"eq",eq2,"tau",tau2,"pa",pa);
#pr2,eq2,tau2,pa=load("Counterfactual2.jld", "pr","eq","tau","pa");

@time pr3,eq3,tau3 = taxreform3(0.3, eq, tau, pa; tol=10.0^-4.0, update = 0.98);
save("Counterfactual3.jld","pr",pr3,"eq",eq3,"tau",tau3,"pa",pa);

@time pr4,eq4,tau4 = taxreform3(0.0, eq, tau, pa);
save("Counterfactual4.jld","pr",pr4,"eq",eq4,"tau",tau4,"pa",pa);

#NOT WORKING
##@time pr1,eq1,tau1 = taxreform1(0.3, eq, tau, pa);
##save("Counterfactual1.jld","pr",pr1,"eq",eq1,"tau",tau1,"pa",pa);
##

("Parameters: delta ",initialParams[1], " ttheta ",initialParams[2], " rhoz ", initialParams[3], " ssigmaz ", initialParams[4],
  " llambda0 ",initialParams[5], " llambda1 ",initialParams[6] , " f ",initialParams[7], "H", initialParams[8]


#################
# Calibration
#################

# Optimization
#         delta   theta   rhoz    sigmaz   lambda0   lambda1    f       H
LB  = [    .05,      .1,   .4,      .01 ,   .01,       .001     0.001   0.001]
#         delta   theta   rhoz    sigmaz   lambda0   lambda1    f       H
UB  = [    .15,     .6,   .95 ,     .2,    1.0,         .06     5.0    5.0 ]

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
