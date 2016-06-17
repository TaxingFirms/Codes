## e=0.0

#Moments(0.1438869857760125,0.12854818418595268,0.3148153122062765,0.08400206845333515,0.017495474248765696,0.0012550009949749603,0.16224651411245722,0.0008443905483383638,-0.007900291424466016,0.0,0.7856350105355269,-2.1897194680585912e6)
#julia> p.E/sum(p.distr)
#0.01640858824055912





fp  = init_firmparameters(hp; e=0.001);

tau = init_taxes() #0.15, 0.3, 0.3, 0.15);
p = guess_prices(tau,fp,hp);
pr  = init_firmproblem(p,tau,fp,hp);



#Compute the model on first time
@time firmVFIParallelOmega!(pr,p,tau,fp); #pr is updated, computes Value Function
#597.841274 seconds on tesla, tol = 10^-3.

#Compute wage such that free entry condition holds
@time w=free_entry!(pr, p, tau, fp,hp,tol=.001)
#1339.480311 seconds on tesla, tol = 10^-2, w = 0.719

#Extract policies and other idiosyncratic results of interest
res=copy_opt_policies(pr);
getpolicies!(res,pr,p,tau,fp);  #r is updated exctracts policies

#Compute mass of entrants and stationary distribution
# both are updated in p.
mass_of_entrants!( res, pr, p, tau, fp);

# Compute aggregate results of interest and moments
aggregates!(res, pr, p, tau, hp, fp);
p.E/sum(p.distr)

#Moments(0.15397590985273868,0.14772653525744459,0.2603363068873757,0.14718752218538547,0.02589398603275505,0.0016620743781697925,0.1784559627915958,0.0009803529756961546,-0.006497613521945203,0.0,0.8474036399284285,-3.2174409005844877e6)
#julia> p.E/sum(p.distr)
#0.01589500284857363
