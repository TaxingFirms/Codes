pa =init_parameters(H=1.07, ddelta=0.125, ttheta=0.2, llambda1 = 0.2, rhoz = 0.75, ssigmaz = 0.072  ,ff= 0.5, e= 0.16);

tau = init_taxes(ttaud =0.15, ttauc= 0.35, ttaui= 0.28, ttaug= 0.15, ttaul=0.28);
@time pr,eq= SolveSteadyState(tau,pa;wguess=0.534, VFItol=10.0^-3.0, displayit0=false, displayw = false);
moments=computeMomentsCutoff(eq.E,pr,eq,tau,pa,cutoffCapital=0.0;toPrint=false);
pcterror_params(pr,eq,tau,pa)


"Focus on investment, sd profits, ave, freq equity issuances, autocov profits, time at work"


# Investment
pa =init_parameters(H=1.07, ddelta=0.09, ttheta=0.2, llambda1 = 0.2, rhoz = 0.75, ssigmaz = 0.072  ,ff= 0.5, e= 0.16);
@time pr,eq= SolveSteadyState(tau,pa;wguess=0.534, VFItol=10.0^-3.0, displayit0=false, displayw = false);
moments=computeMomentsCutoff(eq.E,pr,eq,tau,pa,cutoffCapital=0.0;toPrint=false);
pcterror_params(pr,eq,tau,pa)

#Autocorrelation of shocks
pa =init_parameters(H=1.07, ddelta=0.09, ttheta=0.2, llambda1 = 0.2, rhoz = 0.65, ssigmaz = 0.072  ,ff= 0.5, e= 0.16);
@time pr,eq= SolveSteadyState(tau,pa;wguess=0.534, VFItol=10.0^-3.0, displayit0=false, displayw = false);
moments=computeMomentsCutoff(eq.E,pr,eq,tau,pa,cutoffCapital=0.0;toPrint=false);
pcterror_params(pr,eq,tau,pa)
