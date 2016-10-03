
function decompostion(pr0::FirmProblem, eq0::Equilibrium, tau0::Taxes, pr1::FirmProblem, eq1::Equilibrium, tau1::Taxes, pa::Param)
  #Decompose the effects of a policy change

  #1 Policies
  eq100 = init_equilibirium(eq0.w,tau0,pa);
  eq100.distr= copy(eq0.distr);
  eq100.E=eq0.E;
  pr100  = init_firmproblem(pa, firmvalueguess= pr0.firmvaluegrid);

  firmVFIParallelOmega!(pr100,eq100,tau1,pa; tol = 10.0^-3.0 ); #pr is updated, computes Value Function
  getpolicies!(pr100,eq100,tau1,pa);
  aggregates!(pr100, eq100, tau1, pa; consistencychecks=false);

  #2 Distribution
  eq010 = init_equilibirium(eq0.w,tau0,pa);
  pr010 = deepcopy(pr0);
  eq010.distr= copy(eq1.distr);
  eq010.E= eq1.E;
  aggregates!(pr010, eq010, tau0, pa; consistencychecks=false);

  #3 Prices
  eq001 = init_equilibirium(eq1.w,tau1,pa);
  pr001 = deepcopy(pr0);
  eq001.distr= copy(eq0.distr);
  eq001.E= eq0.E;
  aggregates!(pr001, eq001, tau0, pa; consistencychecks=false);



  labels = ["Output","Labor Demand","Consumption","Turnover", "TFP"];
  ################################################################################
  #Copy after changing params
  initialVec=[eq0.a.output, eq0.a.laborsupply, eq0.a.consumption, eq0.E/sum(eq0.distr), eq0.a.output/(eq0.a.capital^pa.alphak*eq0.a.laborsupply^pa.alphal) ];
  policiesVec=[eq100.a.output, eq100.a.laborsupply, eq100.a.consumption, eq100.E/sum(eq100.distr), eq100.a.output/(eq100.a.capital^pa.alphak*eq100.a.laborsupply^pa.alphal) ]
  distributionVec=[eq010.a.output, eq010.a.laborsupply, eq010.a.consumption, eq010.E/sum(eq010.distr), eq010.a.output/(eq010.a.capital^pa.alphak*eq010.a.laborsupply^pa.alphal) ]
  pricesVec=[eq001.a.output, eq001.a.laborsupply, eq001.a.consumption, eq001.E/sum(eq001.distr), eq001.a.output/(eq001.a.capital^pa.alphak*eq001.a.laborsupply^pa.alphal) ]
  finalVec=[eq1.a.output, eq1.a.laborsupply, eq1.a.consumption, eq1.E/sum(eq1.distr), eq1.a.output/(eq1.a.capital^pa.alphak*eq1.a.laborsupply^pa.alphal) ]


  println(DataFrame(Var=labels, initial=initialVec ,policies =policiesVec , distribution = distributionVec, prices= pricesVec, final=finalVec))

end
