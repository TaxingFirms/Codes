#Variables of interest are GDP, Welfare, TFP, Consumption and Labor.


function taxreform1(tauc::Float64, eq::Equilibrium, tau::Taxes, pa::Param; update::Float64 =0.7, tol::Float64 =10.0^-2.0, maxroutine::Function=maxroutine)

  #Compute tax base for "revenue neutral" reforms
  C = eq.a.collections.c / tau.c
  D = eq.a.collections.d / tau.d
  I = eq.a.collections.i / tau.i

  originalG=eq.a.G;
  wguess= eq.w;


  x= (tauc - tau.c)*C/(D);
  taunew = Taxes(tau.d-x,tauc,tau.i,tau.g,tau.l);
  println("New rates: d = ", taunew.d, " c = ", taunew.c, " i = ", taunew.i, " g = ", taunew.g)

  pr1,eq1=SolveSteadyState(taunew,pa;  wguess = wguess,maxroutine=maxroutine, verbose=false);
  newG=eq1.a.G;


  while abs((originalG - newG)/originalG)>tol
    println("(originalG - newG)/originalG ", (originalG - newG)/originalG)
    tau=deepcopy(taunew);

    D= eq1.a.collections.d / tau.d;
    ntau= update*tau.d + (1-update)*(tau.d + (originalG -newG)/D)
    taunew = Taxes(ntaud,tau.c,tau.i,tau.g,tau.l);
    println("New rates: d = ", taunew.d, " c = ", taunew.c, " i = ", taunew.i, " g = ", taunew.g)

    pr1,eq1= SolveSteadyState(taunew,pa; wguess = wguess,maxroutine=maxroutine, verbose=false);
    newG=eq1.a.G;
  end
  println("(originalG - newG)/originalG",(originalG - newG)/originalG)
 return pr1, eq1, taunew
end





function taxreform2(tauc::Float64, eq::Equilibrium, tau::Taxes, pa::Param; update::Float64 =0.7, tol::Float64 =10.0^-2.0, maxroutine::Function=maxroutine)

  #Compute tax base for "revenue neutral" reforms
  C = eq.a.collections.c / tau.c #corporate base
  D = eq.a.collections.d / tau.d #divideend base
  I = eq.a.collections.i / tau.i #interest base

  originalG=eq.a.G;
  wguess= eq.w;

  x= (tauc - tau.c)*C/D;
  if (tauc)*(1-tau.g-x)>(1-tau.i)
    println("No equilibrium under current taxes")
  else
    taunew = Taxes(tau.d-x,tauc,tau.i,tau.g-x,tau.l);
    println("New rates: d = ", taunew.d, " c = ", taunew.c, " i = ", taunew.i, " g = ", taunew.g)

    #Initiate prices and firm problem, and ultimately, the counterfactual object.

    pr1,eq1=SolveSteadyState(taunew,pa; wguess = wguess,maxroutine=maxroutine, verbose=false)
    newG=eq1.a.G;

    while abs((originalG - newG)/originalG)>tol
      println("(originalG - newG)/originalG", (originalG - newG)/originalG)
      tau=deepcopy(taunew);

      D= eq1.a.collections.d / tau.d;
      ntau= update*tau.d + (1-update)*(tau.d + (originalG -newG)/D);

      if (tau.c)*(1-ntau-x)>(1-tau.i)
        println("No equilibrium under current taxes")
        break
      end
      taunew = Taxes(ntau,tau.c,tau.i,ntau,tau.l);
      println("New rates: d = ", taunew.d, " c = ", taunew.c, " i = ", taunew.i, " g = ", taunew.g)

      pr1,eq1=SolveSteadyState(taunew,pa; wguess = wguess,maxroutine=maxroutine, verbose=false);
      newG=eq1.a.G;
    end
    println("(originalG - newG)/originalG",(originalG - newG)/originalG)
    return pr1, eq1, taunew
  end
end




function taxreform3(tauc::Float64, eq::Equilibrium, tau::Taxes, pa::Param; update::Float64 =0.7, tol::Float64 =10.0^-2.0, maxroutine::Function=maxroutine, verbose::Bool =false)

  #Compute tax base for "revenue neutral" reforms
  C = eq.a.collections.c / tau.c; #corporate base
  D = eq.a.collections.d / tau.d; #dividend base
  I = eq.a.collections.i / tau.i; #interest base

  originalG=eq.a.G;
  wguess= eq.w;

  tauind= (originalG - tauc*C)/(D + I);
  taunew = Taxes(tauind,tauc,tauind,tauind,tau.l);
  println("New rates: d = ", taunew.d, " c = ", taunew.c, " i = ", taunew.i, " g = ", taunew.g)

  #Initiate prices and firm problem, and ultimately, the counterfactual object.

  pr1,eq1=SolveSteadyState(taunew,pa; wguess = wguess,maxroutine=maxroutine, verbose=verbose)
  newG=eq1.a.G;

  while abs((originalG - newG)/originalG)>tol
    println("(originalG - newG)/originalG",(originalG - newG)/originalG)
    tau=deepcopy(taunew);

    D= eq1.a.collections.d / tau.d;
    tauind= (originalG - tauc*C)/(D + I);
    ntauind= update*tau.d + (1-update)*tauind;

    taunew = Taxes(ntauind,tauc,ntauind,ntauind,tau.l);
    println("New rates: d = ", taunew.d, " c = ", taunew.c, " i = ", taunew.i, " g = ", taunew.g)

    pr1,eq1=SolveSteadyState(taunew,pa; wguess = wguess,maxroutine=maxroutine, verbose=false);
    newG=eq1.a.G;
  end
  println("(originalG - newG)/originalG ",(originalG - newG)/originalG)
  return pr1, eq1, taunew
end
