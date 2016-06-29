function SolveModel!(tau::Taxes,fp::FirmParam,hp::HouseholdParam;wguess::Float64=0.69, VFIfunction::Function = firmVFIParallelOmega! )

  eq = init_equilibirium(wguess,tau,pa);
  pr  = init_firmproblem(eq,tau,pa);

  #Compute the model on first time
  @time VFIfunction(pr,eq,tau,pa); #pr is updated, computes Value Function
  #Compute wage such that free entry condition holds
  @time w=free_entry!(eq, pr, tau, pa; xtol=.001)
  #1339.480311 seconds on tesla, tol = 10^-2, w = 0.719

  #Extract policies and other idiosyncratic results of interest
  getpolicies!(pr,eq,tau,pa);  #r is updated exctracts policies

  # Compute mass of entrants and stationary distribution
  # both are updated in p.
  mass_of_entrants!( res, pr, p, tau, fp);

  # Compute aggregate results of interest and moments
  aggregates!(res, pr, p, tau, hp, fp);

  return p,res,pr
end
