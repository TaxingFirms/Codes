function SolveSteadyState(tau::Taxes,pa::Param;wguess::Float64=0.65,
  VFIfunction::Function = firmVFIParallelOmega!, distr_routine::Function = stationarydist,
     maxroutine::Function=maximizationstep ,verbose::Bool=true )

# wguess=0.65; VFIfunction=firmVFIParallelOmega!; distr_routine = stationarydist; maxroutine=maximizationstep; verbose=true;
  eq = init_equilibirium(wguess,tau,pa);
  pr  = init_firmproblem(pa);

  #Compute the model on first time
  VFIfunction(pr,eq,tau,pa; maximizationroutine=maxroutine , verbose = verbose); #pr is updated, computes Value Function

  #Compute wage such that free entry condition holds
  w=free_entry!(eq, pr, tau, pa, VFIfunction, maxroutine ; verbose = verbose)

  #Extract policies and other idiosyncratic results of interest
  getpolicies!(pr,eq,tau,pa);  #r is updated exctracts policies

  # Compute mass of entrants and stationary distribution
  # both are updated in p.
  mass_of_entrantsGHH!( pr, eq, tau, pa, distr_routine ; verbose = false);

  # Compute aggregate results of interest and moments

  aggregates!(pr, eq, tau, pa);

  return pr,eq
end
