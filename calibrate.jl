




function computeDistance(initialParams)
	# Moments are Avg Inv,SD Inv, Avg Lev, SD Lev, Avg Divs, SD Divs, Avg Profits, SD Profits, Avg Eq Iss,
	# Freq Iss, AutoCovProfits
	dataMoments = [.06,.08,.27,.38,.02,.10,.05,.30,-.01,.65,.66]

	# Parameters to be calibrated are, in order
	# Delta  - depreciation
	# rhoz   - autocorrelation of profits
	# sigmaz - sd of profits
	# theta  -  leverage
	# llambda0 - fixed cost of issuance
	# llambda1 - variable cost of issuance

	pa  = init_parameters(ddelta=initialParams[1],rhoz=initialParams[2],ssigmaz=initialParams[3],ttheta=initialParams[4],llambda0=initialParams[5],llambda1=initialParams[6]);
	tau = init_taxes();
	pr,eq= SolveModel!(tau,pa);

	moments = computeMomentsCutoff(eq.E,pr, eq, tau, pa, cutoffCapital=0.0, toPrint=false)

	currentMomentsMatch = [moments.mean_inv_rate,moments.sd_inv_rate,moments.mean_leverage,
	moments.sd_leverage,moments.mean_dividends2k,moments.sd_dividends2k,moments.mean_profits2k,
	moments.sd_profits2k,moments.mean_eqis2k,moments.freq_equis2k,moments.autocov_profits2k]

	sum((currentMomentsMatch-dataMoments).^2)
end

# Optimization
#         delta     rhoz    sigmaz   theta
LB = [      .01,      .5,    .01,    .01 ,   .01, .0001 ]
	   #   theta        rho   sigma   k0       highOrLow
UB  = [     .15,     .95,   .50 ,     .4,    .15, .03  ]

initialGuess = [0.14,0.76,0.0352,0.0034917457985334703,.45,.08,.028]
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


# opt = Opt(:GN_DIRECT_L,length(LB))
#
# lower_bounds!(opt,LB)
# upper_bounds!(opt,UB)
# min_objective!(opt,f)
# xtol_rel!(opt,.1)
#
# (minf,minx,ret) = optimize(opt,initialGuess)
# println("got $minf at $minx after $count iterations (returned $ret)")




opt = Opt(:LN_SBPLX,length(LB))
lower_bounds!(opt,LB)
upper_bounds!(opt,UB)
min_objective!(opt,f)
xtol_rel!(opt,.1)

(minf,minx,ret) = optimize(opt,initialGuess)
println("got $minf at $minx after $count iterations (returned $ret)")