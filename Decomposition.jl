function decomposition(pr0::FirmProblem,eq0::Equilibrium,tau0::Taxes,pr1::FirmProblem,eq1::Equilibrium,tau1::Taxes,pa::Param)
  #Value function under new taxes but initial prices and initial distributions

  labels = ["Output", "Capital", "Labor", "Consumption", "CEV"];
  baseline = [eq0.a.output, eq0.a.capital, eq0.a.laborsupply, eq0.a.consumption, 0.0];
  cev1=cev(eq0,eq1,tau0.l,tau1.l,pa);
  distributioneffect= [eq1.a.output, eq1.a.capital, eq1.a.laborsupply, eq1.a.consumption, cev1];

  eq100 = init_equilibirium(eq0.w,tau,pa);
  eq100.distr= deepcopy(eq0.distr);
  eq100.E=deepcopy(eq0.E);

  valueguess = deepcopy(pr0.firmvaluegrid);
  pr100  = init_firmproblem( pa; firmvalueguess = valueguess);

  firmVFIParallelOmega!(pr100,eq100,tau1,pa; maxroutine=maximizationfast , verbose = true, tol = 10^-3.0 );
  #Extract policies and other idiosyncratic results of interest
  getpolicies!(pr100,eq100,tau1,pa);
  aggregates!(pr100, eq100, tau1, pa; nochecks =true);

  cev100=cev(eq0,eq100,tau0.l,tau1.l,pa);

  policyeffect= [eq100.a.output, eq100.a.capital, eq100.a.laborsupply, eq100.a.consumption, cev100 ];

  println(DataFrame(Variable=labels,Baseline= baseline, Policy = policyeffect, Distribution= distributioneffect ))


  T=50;
  S=500;
  capital0, debt0, networth0, dividends0, investment0, z_history_ind0 = simulation(S, T, pr0,pa ; seed=1234);
  capital100, debt100, networth100, dividends100, investment100, z_history_ind100 = simulation(S, T, pr100,pa ; seed=1234);

  z_history_ind0!=z_history_ind100 && error("Shocks are not the same")

  x=1:T;
  investment= [sum(investment0,2) sum(investment100,2)];
  figure()
  plot(x,investment)
  legend("01", loc="best")

end


function cev(eq0::Equilibrium,eq1::Equilibrium, taul0::Float64, taul1::Float64, pa::Param)
    cev= (eq1.a.consumption - eq0.a.consumption)/eq0.a.consumption - pa.H/(eq0.a.consumption*(1+1/pa.psi))*( (eq1.w*(1-taul1)/pa.H)^(1+pa.psi) - (eq0.w*(1-taul0)/pa.H)^(1+pa.psi) );
end
