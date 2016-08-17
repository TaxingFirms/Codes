#Variables of interest are GDP, Welfare, TFP, Consumption and Labor.


function taxreform1(tauc::Float64, eq::Equilibrium, tau::Taxes, pa::Param; update::Float64 =0.7, tol::Float64 =10.0^-2.0, maxroutine::Function=maximizationfast)

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





function taxreform2(tauc::Float64, eq::Equilibrium, tau::Taxes, pa::Param; update::Float64 =0.7, tol::Float64 =10.0^-2.0, momentsprint::Bool=false, verbose::Bool=false,firmvalueguess::Matrix = repmat(pa.omega.grid,1,pa.Nz) )

#tauc= 0.0; update=0.7; tol= 10.0^-2.0; momentsprint=false; verbose=false; firmvalueguess=copy(pr.firmvaluegrid);


    #Compute tax base for "revenue neutral" reforms
    C = eq.a.collections.c / tau.c #corporate base
    D = eq.a.collections.d / tau.d #divideend base
    I = eq.a.collections.i / tau.i #interest base

    originalG=eq.a.G;
    wguess= eq.w;

    x= (tauc - tau.c)*C/D;

    if (1-tauc)*(1-(tau.g-x))>(1-tau.i)
        error("No equilibrium under current taxes")
    end

    taunew = Taxes(tau.d-x,tauc,tau.i,tau.g-x,tau.l);
    println("New rates: d = ", taunew.d, " c = ", taunew.c, " i = ", taunew.i, " g = ", taunew.g)

    #Initiate prices and firm problem, and ultimately, the counterfactual object.
    pr1,eq1=SolveSteadyState(taunew,pa; wguess = wguess, verbose=verbose, firmvalueguess = firmvalueguess);
    if momentsprint
        moments=computeMomentsCutoff(eq1.E,pr1,eq1,tau,pa,cutoffCapital=0.0,toPrint=momentsprint);
    end

    newG=eq1.a.G;
    newfirmvalueguess=copy(pr1.firmvaluegrid);
    newwguess=eq1.w;

    while abs((originalG - newG)/originalG)>tol
        println("(originalG - newG)/originalG ", (originalG - newG)/originalG)
        tauprime=deepcopy(taunew);

        D= eq1.a.collections.d / tauprime.d;
        ntau= update*tauprime.d + (1-update)*(tauprime.d + (originalG -newG)/D);

    #    if (1-tauprime.c)*(1-ntau)>(1-tau.i)
    #        error("No equilibrium under current taxes")
    #    end
        taunew = Taxes(ntau,tauprime.c,tauprime.i,ntau,tau.l);
        println("New rates: d = ", taunew.d, " c = ", taunew.c, " i = ", taunew.i, " g = ", taunew.g)

        pr1,eq1=SolveSteadyState(taunew,pa; wguess = newwguess, verbose=verbose, firmvalueguess = newfirmvalueguess);
        newG=eq1.a.G;
        newfirmvalueguess=copy(pr1.firmvaluegrid);
        newwguess=eq1.w;
        println("(originalG - newG)/originalG ",(originalG - newG)/originalG)
        println(abs((originalG - newG)/originalG)>tol)
        return pr1, eq1, taunew
    end
end



function taxreform3(tauc::Float64, eq::Equilibrium, tau::Taxes, pa::Param; update::Float64 =0.7, tol::Float64 =10.0^-2.0, momentsprint::Bool=false, verbose::Bool=false,firmvalueguess::Matrix = repmat(pa.omega.grid,1,pa.Nz) )

  #Compute tax base for "revenue neutral" reforms
  C = eq.a.collections.c / tau.c; #corporate base
  D = eq.a.collections.d / tau.d; #dividend base
  I = eq.a.collections.i / tau.i; #interest base

  originalG=eq.a.G;
  wguess= eq.w;

  tauind= (originalG - tauc*C - eq.a.collections.l)/(D + I);
  taunew = Taxes(tauind,tauc,tauind,tauind,tau.l);
  println("New rates: d = ", taunew.d, " c = ", taunew.c, " i = ", taunew.i, " g = ", taunew.g)

  #Initiate prices and firm problem, and ultimately, the counterfactual object.

  pr1,eq1=SolveSteadyState(taunew,pa; wguess = wguess,maxroutine=maxroutine, verbose=verbose)
  newG=eq1.a.G;

  while abs((originalG - newG)/originalG)>tol
    println("(originalG - newG)/originalG",(originalG - newG)/originalG)
    tau=deepcopy(taunew);

    D= eq1.a.collections.d / tau.d;
    tauind= (originalG - tauc*C - eq1.a.collections.l)/(D + I);
    ntauind= update*tau.d + (1-update)*tauind;

    taunew = Taxes(ntauind,tauc,ntauind,ntauind,tau.l);
    println("New rates: d = ", taunew.d, " c = ", taunew.c, " i = ", taunew.i, " g = ", taunew.g)

    pr1,eq1=SolveSteadyState(taunew,pa; wguess = wguess,maxroutine=maxroutine, verbose=false);
    newG=eq1.a.G;
  end
  println("(originalG - newG)/originalG ",(originalG - newG)/originalG)
  return pr1, eq1, taunew
end
