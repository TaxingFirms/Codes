




function computeDistance(initialParams)

	# Freq Iss, AutoCovProfits
	dataMoments = [.063,0.097,0.068,0.105,0.077,0.256,0.509,0.13]

	# Parameters to be calibrated are, in order
	# Delta  - depreciation
	# rhoz   - autocorrelation of profits
	# sigmaz - sd of profits
	# theta  -  leverage
	# llambda0 - fixed cost of issuance
	# llambda1 - variable cost of issuance

	println("Parameters: delta ",initialParams[1], " ttheta ",initialParams[2], " rhoz ", initialParams[3], " ssigmaz ", initialParams[4],
		" llambda0 ",initialParams[5], " llambda1 ",initialParams[6] , " f ",initialParams[7]  )

	pa  = init_parameters(ddelta=initialParams[1], ttheta=initialParams[2],  rhoz=initialParams[3], ssigmaz=initialParams[4],llambda0=initialParams[5],llambda1=initialParams[6], ff=initialParams[7]  );
	tau = init_taxes(ttaud =0.12, ttauc= 0.35, ttaui= 0.29, ttaug= 0.12, ttaul=0.28);
	pr,eq= SolveSteadyState!(tau,pa);

	moments = computeMomentsCutoff(eq.E,pr, eq, tau, pa, cutoffCapital=0.0, toPrint=false)

	currentMomentsMatch = [moments.mean_inv_rate,moments.sd_inv_rate,moments.mean_leverage,
	moments.sd_leverage,moments.mean_dividends2k,moments.sd_dividends2k,moments.mean_profits2k,
	moments.sd_profits2k,moments.mean_eqis2k,moments.freq_equis2k,moments.autocov_profits2k]

	sum((currentMomentsMatch-dataMoments).^2)
end
