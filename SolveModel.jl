function SolveModel!(tau::Taxes,pa::Param;wguess::Float64=0.69, VFIfunction::Function = firmVFIParallelOmega! )

  eq = init_equilibirium(wguess,tau,pa);
  pr  = init_firmproblem(eq,tau,pa);

  #Compute the model on first time
  @time VFIfunction(pr,eq,tau,pa); #pr is updated, computes Value Function

  #Compute wage such that free entry condition holds
  @time w=free_entry!(eq, pr, tau, pa, VFIfunction)

  #Extract policies and other idiosyncratic results of interest
  getpolicies!(pr,eq,tau,pa);  #r is updated exctracts policies

  # Compute mass of entrants and stationary distribution
  # both are updated in p.
  mass_of_entrants!( pr, eq, tau, pa);

  # Compute aggregate results of interest and moments
  aggregates!(pr, eq, tau, pa);

  return pr,eq
end
