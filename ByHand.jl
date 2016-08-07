dataMoments = [.063,0.105,0.211,0.091,0.256,0.509,0.13,0.33];
namesMoments = ["Mean Investment","SD Profits","Mean Leverage",
"Mean Equity Issuance","Frequency of Equity Issuance","Autocovariance Profits",
"Turnover","Time At Work"];
################################################################################
#Copy after changing params
@time pr,eq= SolveSteadyState(tau,pa;wguess=0.68);
moments=computeMomentsCutoff(eq.E,pr,eq,tau,pa,cutoffCapital=0.0);
currentMomentsMatch = [moments.mean_inv_rate,moments.sd_profits2k, moments.mean_leverage,
 moments.mean_eqis2k, moments.freq_equis2k, moments.autocov_profits2k,
moments.turnover, moments.labor];
relErrors=(currentMomentsMatch - dataMoments)./dataMoments;
println(DataFrame(names=namesMoments,moments= currentMomentsMatch, errors=relErrors))

pr.exitprobability #

capital, debt, networth, dividends, investment, z_history_ind = simulation(50, 50,pr,pa; seed=1234);
figure()
plot(capital)

figure()
plot(dividends)

figure()
plot(debt)


# Calibration Board presentation
pa  = init_parameters( H=1.3, ff= 0.015, llambda0=0.02, llambda1= 0.04, ddelta = 0.12, ttheta = 0.42,rhoz= 0.76, ssigmaz= 0.0352);
tau = init_taxes(ttaud =0.12, ttauc= 0.35, ttaui= 0.29, ttaug= 0.12, ttaul=0.28);

# Very nice
pa  = init_parameters( H=1.3, ff= 0.15, llambda0=0.02, llambda1= 0.04, ddelta = 0.14, ttheta = 0.3,rhoz= 0.76, ssigmaz= 0.0385);
comment:"Simulations are decent for large shocks. The selection mechanism is not doing much: productive firms stay and unproductive firms exit,
  independent of size. Investment is a little bit low."

#Decrease support of shocks
comment:"Turnover decreases, decreasing mean and freq equity shares. And increasing investment. We get some equity issuances by incumbents"
# Recalibrate
pa  = init_parameters( H=1.3, ff= 0.15, llambda0=0.02, llambda1= 0.04, ddelta = 0.14, ttheta = 0.3,rhoz= 0.76, ssigmaz= 0.045);
# The wage decreases substancially and profits increase
comment:"Once recalibrated, the moments look good (except the autocorr of profits, which increases but is still too low)"

# Decrease fix cost but increase e
pa  = init_parameters( H=1.3, ff= 0.05, llambda0=0.02, llambda1= 0.04, ddelta = 0.14, ttheta = 0.3,rhoz= 0.76, ssigmaz= 0.045, e=0.08);
comment:"Turnover decreases, decreasing mean and freq equity shares. And increasing investment. We get some equity issuances by incumbents
wage is 0.55"
pa  = init_parameters( H=1.3, ff= 0.05, llambda0=0.02, llambda1= 0.04, ddelta = 0.095, ttheta = 0.3,rhoz= 0.76, ssigmaz= 0.09, e=0.08);
comment:"Moments are very nice: but growing out of constraints mechanism is not present. Wage =0.67, Which helps with ave profits "

pa  = init_parameters( H=1.5, ff= 0.01, llambda0=0.06, llambda1= 0.04, ddelta = 0.095, ttheta = 0.25,rhoz= 0.76, ssigmaz= 0.09, e=0.08);
