function SolveModel!(tau::Taxes,fp::FirmParam,hp::HouseholdParam;wguess::Float64=0.69)

  p = init_equilibirium(wguess,tau,fp,hp);
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

  # Compute mass of entrants and stationary distribution
  # both are updated in p.
  mass_of_entrants!( res, pr, p, tau, fp);

  # Compute aggregate results of interest and moments
  aggregates!(res, pr, p, tau, hp, fp);

  return p,res,pr
end
