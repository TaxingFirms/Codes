function SolveSteadyState(tau::Taxes,pa::Param;wguess::Float64=0.65, VFItol =10.0^-3.0,
  VFIfunction::Function = firmVFIParallelOmega!, distr_routine::Function = stationarydist,
     maxroutine::Function=maximizationfast, displayit0::Bool=true, displayw::Bool=true, displayitt::Bool=false, firmvalueguess::Matrix = repmat(pa.omega.grid,1,pa.Nz) )

# wguess=0.54; VFItol =10.0^-3.0; VFIfunction=firmVFIParallelOmega!; distr_routine = stationarydist; maxroutine=maximizationstep; verbose=true; firmvalueguess = repmat(pa.omega.grid,1,pa.Nz);
  eq = init_equilibirium(wguess,tau,pa);
  pr  = init_firmproblem(pa, firmvalueguess= firmvalueguess);

  #Compute the model on first time
  VFIfunction(pr,eq,tau,pa; maxroutine=maxroutine , verbose = displayit0, tol = VFItol ); #pr is updated, computes Value Function

  #Compute wage such that free entry condition holds
  w=free_entry!(eq, pr, tau, pa, VFIfunction, maxroutine, VFItol ; displayw = displayw , displayit =displayitt)

  #Extract policies and other idiosyncratic results of interest
  getpolicies!(pr,eq,tau,pa);  #pr is updated exctracts policies

  # Compute mass of entrants and stationary distribution
  # both are updated in eq.
  mass_of_entrantsGHH!( pr, eq, tau, pa, distr_routine ; verbose = false);

  # Compute aggregate results of interest and moments

  aggregates!(pr, eq, tau, pa);

  return pr,eq
end
