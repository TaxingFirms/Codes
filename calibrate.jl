




function computeDistance(initialParams)

	# Freq Iss, AutoCovProfits
	dataMoments = [.063,0.105,0.211,0.091,0.256,0.509,0.13,0.33]

	# Parameters to be calibrated are, in order
	# Delta  - depreciation
	# rhoz   - autocorrelation of profits
	# sigmaz - sd of profits
	# theta  -  leverage
	# llambda0 - fixed cost of issuance
	# llambda1 - variable cost of issuance

	println("Parameters: delta ",initialParams[1], " ttheta ",initialParams[2], " rhoz ", initialParams[3], " ssigmaz ", initialParams[4],
		" llambda0 ",initialParams[5], " llambda1 ",initialParams[6] , " f ",initialParams[7], " H ", initialParams[8] )

	factorToBenchmark = .5
	Nnz = round(Int64,factorToBenchmark*9); Nnk = round(Int64,factorToBenchmark*80); Nnq = round(Int64,factorToBenchmark*40); Nnomega = round(Int64,factorToBenchmark*100);
	pa  = init_parameters(ddelta=initialParams[1],ttheta=initialParams[2],rhoz=initialParams[3],
	ssigmaz=initialParams[4],llambda0=initialParams[5],llambda1=initialParams[6],  ff =initialParams[7], H= initialParams[8],  Nz=Nnz,Nk=Nnk,Nomega=Nnomega, Nq=Nnq);
	tau = init_taxes(ttaud =0.15, ttauc= 0.35, ttaui= 0.29, ttaug= 0.15, ttaul=0.28);

	pr,eq= SolveSteadyState(tau,pa, VFItol=10.0^-3.0, verbose=false);
	moments = computeMomentsCutoff(eq.E,pr, eq, tau, pa, cutoffCapital=0.0, toPrint=false)

	currentMomentsMatch = [moments.mean_inv_rate,moments.sd_profits2k, moments.mean_leverage,
	moments.mean_eqis2k, moments.freq_equis2k, moments.autocov_profits2k,
	moments.turnover, moments.labor]

	sum((currentMomentsMatch-dataMoments).^2.0)
end


function pcterror_params(pr::FirmProblem,eq::Equilibrium, tau::Taxes, pa::Param)
  dataMoments = [.063,0.105,0.211,0.091,0.256,0.509,0.13,0.33];
  namesMoments = ["Mean Investment","SD Profits","Mean Leverage",
  "Mean Equity Issuance","Frequency of Equity Issuance","Autocovariance Profits",
  "Turnover","Time At Work"];
  ################################################################################
  #Copy after changing params
  moments=computeMomentsCutoff(eq.E,pr,eq,tau,pa,cutoffCapital=0.0);
  currentMomentsMatch = [moments.mean_inv_rate,moments.sd_profits2k, moments.mean_leverage,
   moments.mean_eqis2k, moments.freq_equis2k, moments.autocov_profits2k,
  moments.turnover, moments.labor,(eq.a.collections.c/tau.c)/(eq.a.collections.d/tau.d)];
  relErrors=(currentMomentsMatch - dataMoments)./dataMoments;
  println(DataFrame(names=namesMoments,moments= currentMomentsMatch, errors=relErrors))
end
