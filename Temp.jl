

pa =init_parameters( ddelta=0.1, ttheta=0.525, rhoz=0.48, ssigmaz=0.25, llambda0=0.25, llambda1=0.04505, ff=0.5,H=5.005);
tau = init_taxes(ttaud =0.15, ttauc= 0.35, ttaui= 0.29, ttaug= 0.15, ttaul=0.28);
@time pr,eq= SolveSteadyState(tau,pa;wguess=0.54, VFItol=10.0^-3.0, verbose=true);
