#Variables of interest are GDP, Welfare, TFP, Consumption and Labor.


function taxreform1(tauc::Float64, p::Equilibrium, tau::Taxes, fp::FirmParam, hp::HouseholdParam)

  #Compute tax base for "revenue neutral" reforms
  C = p.a.collections.c / tau.c #corporate base
  D = p.a.collections.d / tau.d #corporate base
  I = p.a.collections.i / tau.i #corporate base
  CG = p.a.collections.g / tau.g #corporate base

  originalG=p.a.G;


  x= (tauc - tau.c)*C/D;
  taunew = Taxes(tau.d-x,tauc,tau.i,tau.g);
  println("New rates: d = ", taunew.d, " c = ", taunew.c, " i = ", taunew.i, " g = ", taunew.g)

  p,res,pr=SolveModel!(taunew,fp,hp)
  save("/home/dwills/firms/counterfactual1.jld", "pr", pr, "tau", taunew, "fp", fp, "res",res, "p",p);

  newG=p.a.G;

  while abs(originalG -newG)<10^-3.0
    tau=deepcopy(taunew);

    D= p.a.collections.d / tau.d;
    ntaud= (newG - originalG)/D;
    taunew = Taxes(ntaud,tau.c,tau.i,tau.g);
    println("New rates: d = ", taunew.d, " c = ", taunew.c, " i = ", taunew.i, " g = ", taunew.g)

    originalG=newG;
    p,res,pr= SolveModel!(taunew,fp,hp);
    newG=p.a.G;
  end

  save("/home/dwills/firms/counterfactual2.jld", "pr", pr, "tau", taunew, "fp", fp, "res",res, "p",p);

end





function taxreform2(tauc::Float64, p::Equilibrium, tau::Taxes, fp::FirmParam, hp::HouseholdParam)

  #Compute tax base for "revenue neutral" reforms
  C = p.a.collections.c / tau.c #corporate base
  D = p.a.collections.d / tau.d #corporate base
  I = p.a.collections.i / tau.i #corporate base
  G = p.a.collections.g / tau.g #corporate base

  x= (tauc - tau.c)*C/(D+G);
  taunew = Taxes(tau.d-x,tauc,tau.i,tau.g-x);
  println("New rates: d = ", taunew.d, " c = ", taunew.c, " i = ", taunew.i, " g = ", taunew.g)

  #Initiate prices and firm problem, and ultimately, the counterfactual object.

  SolveModel!(tau,fp,hp)
  save("/home/dwills/firms/counterfactual2.jld", "pr", pr, "tau", taunew, "fp", fp, "res",res, "p",p);

  newG=p.a.G;

  while abs(originalG -newG)<10^-3.0
    tau=deepcopy(taunew);

    D= p.a.collections.d / tau.d;
    ntaud= (newG - originalG)/D;
    taunew = Taxes(ntaud,tau.c,tau.i,tau.g);
    println("New rates: d = ", taunew.d, " c = ", taunew.c, " i = ", taunew.i, " g = ", taunew.g)

    originalG=newG;
    SolveModel!(taunew,fp,hp; wguess=p.w);
    newG=p.a.G;
  end
  save("/home/dwills/firms/counterfactual2.jld", "pr", pr, "tau", taunew, "fp", fp, "res",res, "p",p);
end
