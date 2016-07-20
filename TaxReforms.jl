#Variables of interest are GDP, Welfare, TFP, Consumption and Labor.


function taxreform1(tauc::Float64, eq::Equilibrium, tau::Taxes, pa::Param; update =0.7, tol =10.0^-6.0)

  #Compute tax base for "revenue neutral" reforms
  C = eq.a.collections.c / tau.c
  D = eq.a.collections.d / tau.d
  I = eq.a.collections.i / tau.i

  originalG=eq.a.G;
  wguess= eq.w;


  x= (tauc - tau.c)*C/(D);
  taunew = Taxes(tau.d-x,tauc,tau.i,tau.g);
  println("New rates: d = ", taunew.d, " c = ", taunew.c, " i = ", taunew.i, " g = ", taunew.g)

  pr1,eq1=SolveSteadyState(taunew,pa;  wguess = wguess, verbose=false);
  newG=eq1.a.G;


  while abs(originalG -newG)>tol
    println("originalG - newG ", originalG -newG)
    tau=deepcopy(taunew);

    D= eq1.a.collections.d / tau.d;
    ntau= update*tau.d + (1-update)*(tau.d + (originalG -newG)/D)
    taunew = Taxes(ntaud,tau.c,tau.i,tau.g);
    println("New rates: d = ", taunew.d, " c = ", taunew.c, " i = ", taunew.i, " g = ", taunew.g)

    pr1,eq1= SolveSteadyState(taunew,pa; wguess = wguess, verbose=false);
    newG=eq1.a.G;
  end
  println("originalG - newG ", originalG -newG)
 return pr1, eq1, taunew
end





function taxreform2(tauc::Float64, eq::Equilibrium, tau::Taxes, pa::Param; update =0.7, tol =10.0^-6.0)

  #Compute tax base for "revenue neutral" reforms
  C = eq.a.collections.c / tau.c #corporate base
  D = eq.a.collections.d / tau.d #divideend base
  I = eq.a.collections.i / tau.i #interest base

  originalG=eq.a.G;
  wguess= eq.w;

  x= (tauc - tau.c)*C/D;
  taunew = Taxes(tau.d-x,tauc,tau.i,tau.g-x);
  println("New rates: d = ", taunew.d, " c = ", taunew.c, " i = ", taunew.i, " g = ", taunew.g)

  #Initiate prices and firm problem, and ultimately, the counterfactual object.

  pr1,eq1=SolveSteadyState(taunew,pa; wguess = wguess, verbose=false)
  newG=eq1.a.G;

  while abs(originalG -newG)>tol
    println("originalG - newG ", originalG -newG)
    tau=deepcopy(taunew);

    D= eq1.a.collections.d / tau.d;
    ntau= update*tau.d + (1-update)*(tau.d + (originalG -newG)/D);
    taunew = Taxes(ntau,tau.c,tau.i,ntau);
    println("New rates: d = ", taunew.d, " c = ", taunew.c, " i = ", taunew.i, " g = ", taunew.g)

    pr1,eq1=SolveSteadyState(taunew,pa; wguess = wguess, verbose=false);
    newG=eq1.a.G;
  end
  println("originalG - newG ", originalG -newG)
  return pr1, eq1, taunew
end




function taxreform3(tauc::Float64, eq::Equilibrium, tau::Taxes, pa::Param; update =0.7, tol =10.0^-6.0, verbose =false)

  #Compute tax base for "revenue neutral" reforms
  C = eq.a.collections.c / tau.c; #corporate base
  D = eq.a.collections.d / tau.d; #dividend base
  I = eq.a.collections.i / tau.i; #interest base

  originalG=eq.a.G;
  wguess= eq.w;

  tauind= (originalG - tauc*C)/(D + I);
  taunew = Taxes(tauind,tauc,tauind,tauind);
  println("New rates: d = ", taunew.d, " c = ", taunew.c, " i = ", taunew.i, " g = ", taunew.g)

  #Initiate prices and firm problem, and ultimately, the counterfactual object.

  pr1,eq1=SolveSteadyState(taunew,pa; wguess = wguess, verbose=verbose)
  newG=eq1.a.G;

  while abs(originalG -newG)>tol
    println("originalG - newG ", originalG -newG)
    tau=deepcopy(taunew);

    D= eq1.a.collections.d / tau.d;
    tauind= (originalG - tauc*C)/(D + I);
    ntauind= update*tau.d + (1-update)*tauind;

    taunew = Taxes(ntauind,tauc,ntauind,ntauind);
    println("New rates: d = ", taunew.d, " c = ", taunew.c, " i = ", taunew.i, " g = ", taunew.g)

    pr1,eq1=SolveSteadyState(taunew,pa; wguess = wguess, verbose=false);
    newG=eq1.a.G;
  end
  println("originalG - newG ", originalG -newG)
  return pr1, eq1, taunew
end
