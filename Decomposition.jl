
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
  #eq110 = init_equilibirium(eq0.w,tau1,pa);
  #pr110=deepcopy(pr100);
  #mass_of_entrantsGHH!( pr110, eq110, tau1, pa , stationarydist; verbose = false);
  #aggregates!(pr110, eq110, tau1, pa, ; consistencychecks=false);

  #3 Prices
  #eq101 = init_equilibirium(eq1.w,tau1,pa);
  #pr101 = deepcopy(pr1);
  #eq101.distr= copy(eq0.distr);
  #eq101.E= eq0.E;
  #aggregates!(pr101, eq101, tau1, pa; consistencychecks=false);
  eq101 = init_equilibirium(eq0.w,tau0,pa);
  pr101 = deepcopy(pr0);
  eq101.distr= copy(eq1.distr);
  eq101.E= eq1.E;
  aggregates!(pr101, eq101, tau1, pa; consistencychecks=false);



  labels = ["Output","Labor Demand","Consumption","Turnover", "TFP"];
  ################################################################################
  #Copy after changing params
  initialVec=[eq0.a.output, eq0.a.laborsupply, eq0.a.consumption, eq0.E/sum(eq0.distr), eq0.a.output/(eq0.a.capital^pa.alphak*eq0.a.laborsupply^pa.alphal) ];
  taxesVec=[eq100.a.output, eq100.a.laborsupply, eq100.a.consumption, eq100.E/sum(eq100.distr), eq100.a.output/(eq100.a.capital^pa.alphak*eq100.a.laborsupply^pa.alphal) ]
#  dist0Vec=[eq110.a.output, eq110.a.laborsupply, eq110.a.consumption, eq110.E/sum(eq110.distr), eq110.a.output/(eq110.a.capital^pa.alphak*eq110.a.laborsupply^pa.alphal) ]
  dist1Vec=[eq101.a.output, eq101.a.laborsupply, eq101.a.consumption, eq101.E/sum(eq101.distr), eq101.a.output/(eq101.a.capital^pa.alphak*eq101.a.laborsupply^pa.alphal) ]
  finalVec=[eq1.a.output, eq1.a.laborsupply, eq1.a.consumption, eq1.E/sum(eq1.distr), eq1.a.output/(eq1.a.capital^pa.alphak*eq1.a.laborsupply^pa.alphal) ]


  println(DataFrame(Var=labels, initial=initialVec ,taxes =taxesVec , distribution1 = dist1Vec,final=finalVec))

end
