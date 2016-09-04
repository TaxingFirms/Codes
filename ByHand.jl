pa =init_parameters( H=1.176, bbeta=0.972, ff= 0.5, aalphak=0.23, aalphal=0.64, llambda0=0.004, llambda1= 0.04, ddelta = 0.13, allowance=0.86, ttheta = 0.25,rhoz= 0.75, ssigmaz= 0.08, e=0.0, A=1.0);

@time pr,eq= SolveSteadyState(tau,pa;wguess=0.54, VFItol=10.0^-3.0, verbose=false);
moments=computeMomentsCutoff(eq.E,pr,eq,tau,pa,cutoffCapital=0.0;toPrint=false);
pcterror_params(pr,eq,tau,pa)


pa =init_parameters( H=1.176, bbeta=0.972, ff= 0.5, aalphak=0.23, aalphal=0.64, llambda0=0.004, llambda1= 0.5, ddelta = 0.13, allowance=0.86, ttheta = 0.25,rhoz= 0.75, ssigmaz= 0.08, e=0.0, A=1.0);

@time pr,eq= SolveSteadyState(tau,pa;wguess=0.54, VFItol=10.0^-3.0, verbose=false);
moments=computeMomentsCutoff(eq.E,pr,eq,tau,pa,cutoffCapital=0.0;toPrint=false);
pcterror_params(pr,eq,tau,pa)


pa =init_parameters( H=1.176, bbeta=0.972, ff= 0.5, aalphak=0.23, aalphal=0.64, llambda0=0.004, llambda1= 0.2, ddelta = 0.13, allowance=0.86, ttheta = 0.25,rhoz= 0.75, ssigmaz= 0.08, e=0.0, A=1.0);

@time pr,eq= SolveSteadyState(tau,pa;wguess=eq.w, VFItol=10.0^-3.0, verbose=false);
moments=computeMomentsCutoff(eq.E,pr,eq,tau,pa,cutoffCapital=0.0;toPrint=false);
pcterror_params(pr,eq,tau,pa)

pa =init_parameters( H=1.176, bbeta=0.972, ff= 0.5, aalphak=0.23, aalphal=0.64, llambda0=0.004, llambda1= 0.2, ddelta = 0.13, allowance=0.86, ttheta = 0.25,rhoz= 0.75, ssigmaz= 0.08, e=0.0, A=1.0);
