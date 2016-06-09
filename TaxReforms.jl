#Variables of interest are GDP, Welfare, TFP, Consumption and Labor.

type TaxEquilibrium
  tau::Taxes
  pr::FirmProblem
  p::Equilibrium
end

function taxreform1(tauc::Float64, p::Equilibrium, tau::Taxes, fp::FirmParam, hp::HouseholdParam)

  #Compute tax base for "revenue neutral" reforms
  C = p.collections.c / tau.c #corporate base
  D = p.collections.d / tau.d #corporate base
  I = p.collections.i / tau.i #corporate base
  G = p.collections.g / tau.g #corporate base

  x= (tauc - tau.c)*C/(D+I+G);
  taunew = Taxes(tau.d-x,tauc,tau.i-x,tau.g-x);
  println("New rates: d = ", taunew.d, " c = ", taunew.c, " i = ", taunew.i, " g = ", taunew.g   )

  #Initiate prices and firm problem, and ultimately, the counterfactual object.
  pnew = deepcopy(p);
  prnew  = init_firmproblem(pnew,taunew,fp,hp);
  counterfactual = TaxEquilibrium(taunew,prnew,pnew)

  taxequilibrium!(counterfactual, fp, hp)

  return counterfactual
end

function taxreform2(taud::Float64, p::Equilibrium, tau::Taxes, fp::FirmParam, hp::HouseholdParam)

  #Compute tax base for "revenue neutral" reforms
  C = p.collections.c / tau.c #corporate base
  D = p.collections.d / tau.d #corporate base
  I = p.collections.i / tau.i #corporate base
  G = p.collections.g / tau.g #corporate base

  x= (tau.d - taud)*C/(D+I+G);
  taunew = Taxes(taud,tau.c-x,tau.i-x,tau.g-x);
  println("New rates: d = ", taunew.d, " c = ", taunew.c, " i = ", taunew.i, " g = ", taunew.g   )

  #Initiate prices and firm problem, and ultimately, the counterfactual object.
  pnew = deepcopy(p);
  prnew  = init_firmproblem(pnew,taunew,fp,hp);
  counterfactual = TaxEquilibrium(taunew,prnew,pnew)

  taxequilibrium!(counterfactual, fp, hp)

  return counterfactual
end

function taxequilibrium!(cf::TaxEquilibrium, fp::FirmParam, hp::HouseholdParam)
#Computes equilibrium with new tax rates and fills pnew and returns prnew

  @time firmVFIParallel!(cf.pr,cf.p,cf.tau,fp); #pr is updated, computes Value Function
  @time w=free_entry!(cf.pr, cf.p, cf.tau, fp, hp); #p is updated

  #Extract policies and other idiosyncratic results of interest
  res=copy_opt_policies(cf.pr);
  getpolicies!(res,cf.pr,cf.p,cf.tau,fp);  #r is updated exctracts policies

  #Compute invariant distribution for E and compute aggregate results of interest
  equilibrium_aggregates!( res, cf.pr, cf.p, cf.tau, fp);

end
