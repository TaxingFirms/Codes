




function computeDistance(initialParams)
	# Moments are Avg Inv,SD Inv, Avg Lev, SD Lev, Avg Divs, SD Divs, Avg Profits, SD Profits, Avg Eq Iss,
	# Freq Iss, AutoCovProfits
	dataMoments = [.06,.08,.27,.38,.02,.10,.05,.30,.081,.65,.66]

	# Parameters to be calibrated are, in order
	# Delta  - depreciation
	# rhoz   - autocorrelation of profits
	# sigmaz - sd of profits
	# theta  -  leverage
	# llambda0 - fixed cost of issuance
	# llambda1 - variable cost of issuance

	println("Parameters: delta ",initialParams[1]," rhoz ", initialParams[2], " ssigmaz ", initialParams[3], " ttheta ",initialParams[4],
		" llambda0 ",initialParams[5], " llambda1 ",initialParams[6])

	pa  = init_parameters(ddelta=initialParams[1],rhoz=initialParams[2],ssigmaz=initialParams[3],ttheta=initialParams[4],llambda0=initialParams[5],llambda1=initialParams[6]);
	tau = init_taxes();
	pr,eq= SolveModel!(tau,pa);

	moments = computeMomentsCutoff(eq.E,pr, eq, tau, pa, cutoffCapital=0.0, toPrint=false)

	currentMomentsMatch = [moments.mean_inv_rate,moments.sd_inv_rate,moments.mean_leverage,
	moments.sd_leverage,moments.mean_dividends2k,moments.sd_dividends2k,moments.mean_profits2k,
	moments.sd_profits2k,moments.mean_eqis2k,moments.freq_equis2k,moments.autocov_profits2k]

	sum((currentMomentsMatch-dataMoments).^2)
end
