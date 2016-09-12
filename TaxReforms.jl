
function taxreform1(tauc::Float64, eq::Equilibrium, tau::Taxes, pa::Param; update::Float64 =0.7, tol::Float64 =10.0^-3.0, momentsprint::Bool=false, verbose::Bool=false,firmvalueguess::Matrix = repmat(pa.omega.grid,1,pa.Nz))
# Gets a new level of tauc. CLoses using the dividend tax

  #Compute tax base for "revenue neutral" reforms
  C = eq.a.collections.c / tau.c
  D = eq.a.collections.d / tau.d
  I = eq.a.collections.i / tau.i

  originalG=eq.a.G;
  wguess= eq.w;

  x= (tauc - tau.c)*C/D;
  taunew = Taxes(tau.d-x,tauc,tau.i,tau.g,tau.l);
  println("New rates: d = ", taunew.d, " c = ", taunew.c, " i = ", taunew.i, " g = ", taunew.g)

  pr1,eq1=SolveSteadyState(taunew,pa;  wguess = wguess,verbose=false);
  if momentsprint
      moments=computeMomentsCutoff(eq1.E,pr1,eq1,tau,pa,cutoffCapital=0.0,toPrint=momentsprint);
  end

  newG=eq1.a.G;

  while abs((originalG - newG)/originalG)>tol
    println("(originalG - newG)/originalG ", (originalG - newG)/originalG)
    tauprime=deepcopy(taunew);

    D= eq1.a.collections.d / tauprime.d;
    ntaud= update*tauprime.d + (1-update)*(tauprime.d + (originalG -newG)/D)
    taunew = Taxes(ntaud,tauprime.c,tauprime.i,tauprime.g,tauprime.l);
    println("New rates: d = ", taunew.d, " c = ", taunew.c, " i = ", taunew.i, " g = ", taunew.g)

    pr1,eq1= SolveSteadyState(taunew,pa; wguess = wguess, verbose=false);
    if momentsprint
        moments=computeMomentsCutoff(eq1.E,pr1,eq1,tau,pa,cutoffCapital=0.0,toPrint=momentsprint);
    end

    newG=eq1.a.G;
  end
  println("(originalG - newG)/originalG",(originalG - newG)/originalG)
 return pr1, eq1, taunew
end



function taxreform2(tauc::Float64, eq::Equilibrium, tau::Taxes, pa::Param; update::Float64 =0.7, tol::Float64 =10.0^-3.0, momentsprint::Bool=false, verbose::Bool=false,firmvalueguess::Matrix = repmat(pa.omega.grid,1,pa.Nz) )
# Gets a new level of tauc. CLoses using the dividend and capital gains taxes at the same time.

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
    pr1,eq1=SolveSteadyState(taunew,pa; wguess = wguess, verbose=verbose);
    if momentsprint
        moments=computeMomentsCutoff(eq1.E,pr1,eq1,tau,pa,cutoffCapital=0.0,toPrint=momentsprint);
    end

    newG=eq1.a.G;
    while abs((originalG - newG)/originalG)>tol
        println("(originalG - newG)/originalG ", (originalG - newG)/originalG)
        tauprime=deepcopy(taunew);

        D= eq1.a.collections.d / tauprime.d;
        ntau= update*tauprime.d + (1-update)*(tauprime.d + (originalG -newG)/D);

        if (1-tauprime.c)*(1-ntau)>(1-tau.i)
            error("No equilibrium under current taxes")
        end
        taunew = Taxes(ntau,tauprime.c,tauprime.i,ntau,tau.l);
        println("New rates: d = ", taunew.d, " c = ", taunew.c, " i = ", taunew.i, " g = ", taunew.g)

        pr1,eq1=SolveSteadyState(taunew,pa; wguess = wguess, verbose=verbose);
        if momentsprint
            moments=computeMomentsCutoff(eq1.E,pr1,eq1,tau,pa,cutoffCapital=0.0,toPrint=momentsprint);
        end
        newG=eq1.a.G;
    end
  println("(originalG - newG)/originalG ",(originalG - newG)/originalG)
  return pr1, eq1, taunew
end



function taxreform3(tauc::Float64, eq::Equilibrium, tau::Taxes, pa::Param; update::Float64 =0.7, tol::Float64 =10.0^-3.0, momentsprint::Bool=false, verbose::Bool=false,firmvalueguess::Matrix = repmat(pa.omega.grid,1,pa.Nz) )
# Gets a new level of tauc. Closes using dividend = capital gains = interest.


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
  pr1,eq1=SolveSteadyState(taunew,pa; wguess = wguess, verbose=verbose);
  if momentsprint
      moments=computeMomentsCutoff(eq1.E,pr1,eq1,tau,pa,cutoffCapital=0.0,toPrint=momentsprint);
  end
  newG=eq1.a.G;

  while abs((originalG - newG)/originalG)>tol
    println("(originalG - newG)/originalG ",(originalG - newG)/originalG)
    tauprime=deepcopy(taunew);

    D= eq1.a.collections.d / tauprime.d;
    I= eq1.a.collections.i / tauprime.i;
    tauind= (originalG - tauc*C - eq1.a.collections.l)/(D + I);
    ntauind= update*tauprime.d + (1-update)*tauind;

    taunew = Taxes(ntauind,tauc,ntauind,ntauind,tau.l);
    println("New rates: d = ", taunew.d, " c = ", taunew.c, " i = ", taunew.i, " g = ", taunew.g)

    pr1,eq1=SolveSteadyState(taunew,pa; wguess = wguess, verbose=false);
    if momentsprint
        moments=computeMomentsCutoff(eq1.E,pr1,eq1,tau,pa,cutoffCapital=0.0,toPrint=momentsprint);
    end
    newG=eq1.a.G;
  end
  println("(originalG - newG)/originalG ",(originalG - newG)/originalG)
  return pr1, eq1, taunew
end


function taxreform4(tauc::Float64, eq::Equilibrium, tau::Taxes, pa::Param; update::Float64 =0.7, tol::Float64 =10.0^-3.0, momentsprint::Bool=false, verbose::Bool=false,firmvalueguess::Matrix = repmat(pa.omega.grid,1,pa.Nz) )

  #Compute tax base for "revenue neutral" reforms
  C = eq.a.collections.c / tau.c; #corporate base
  D = eq.a.collections.d / tau.d; #dividend base
  I = eq.a.collections.i / tau.i; #interest base

  originalG=eq.a.G;
  wguess= eq.w;

  tauinew=1-(1-tau.g)*(1-tauc) - 10.0^-3.0
  taudnew = (originalG - tauc*C - eq.a.collections.l - tauinew*I)/D ;
  taunew = Taxes(taudnew,tauc,tauinew,tau.g,tau.l);
  println("New rates: d = ", taunew.d, " c = ", taunew.c, " i = ", taunew.i, " g = ", taunew.g)


  #Initiate prices and firm problem, and ultimately, the counterfactual object.
  pr1,eq1=SolveSteadyState(taunew,pa; wguess = wguess, verbose=verbose);
  if momentsprint
      moments=computeMomentsCutoff(eq1.E,pr1,eq1,tau,pa,cutoffCapital=0.0,toPrint=momentsprint);
  end
  newG=eq1.a.G;

  while abs((originalG - newG)/originalG)>tol
    println("(originalG - newG)/originalG ",(originalG - newG)/originalG)
    tauprime=deepcopy(taunew);

    D= eq1.a.collections.d / tauprime.d;
    I= eq1.a.collections.i / tauprime.i;
    tauinew=1-(1-tau.g)*(1-tauc) - 10.0^-3.0
    taudnew = (originalG - tauc*C - eq.a.collections.l - tauinew*I)/D ;

    taudnew2= update*tauprime.d + (1-update)*taudnew;

    taunew = Taxes(taudnew2,tauc,tauinew,tau.g,tau.l);
    println("New rates: d = ", taunew.d, " c = ", taunew.c, " i = ", taunew.i, " g = ", taunew.g)

    pr1,eq1=SolveSteadyState(taunew,pa; wguess = wguess, verbose=false);
    if momentsprint
        moments=computeMomentsCutoff(eq1.E,pr1,eq1,tau,pa,cutoffCapital=0.0,toPrint=momentsprint);
    end
    newG=eq1.a.G;
  end
  println("(originalG - newG)/originalG ",(originalG - newG)/originalG)
  return pr1, eq1, taunew
end


function taxreform5(tauc::Float64, eq::Equilibrium, tau::Taxes, pa::Param; update::Float64 =0.7, tol::Float64 =10.0^-3.0, momentsprint::Bool=false, verbose::Bool=false,firmvalueguess::Matrix = repmat(pa.omega.grid,1,pa.Nz) )

  #Compute tax base for "revenue neutral" reforms
  C = eq.a.collections.c / tau.c; #corporate base
  D = eq.a.collections.d / tau.d; #dividend base
  I = eq.a.collections.i / tau.i; #interest base

  originalG=eq.a.G - eq.a.collections.l;
  wguess= eq.w;

  tauind= (originalG - tauc*C)/(D + I);
  taunew = Taxes(tauind,tauc,tauind,tauind,tau.l);
  println("New rates: d = ", taunew.d, " c = ", taunew.c, " i = ", taunew.i, " g = ", taunew.g)


  #Initiate prices and firm problem, and ultimately, the counterfactual object.
  pr1,eq1=SolveSteadyState(taunew,pa; wguess = wguess, verbose=verbose);
  if momentsprint
      moments=computeMomentsCutoff(eq1.E,pr1,eq1,tau,pa,cutoffCapital=0.0,toPrint=momentsprint);
  end
  newG=eq1.a.G - eq1.a.collections.l ;

  while abs((originalG - newG)/originalG)>tol
    println("(originalG - newG)/originalG ",(originalG - newG)/originalG)
    tauprime=deepcopy(taunew);

    D= eq1.a.collections.d / tauprime.d;
    I= eq1.a.collections.i / tauprime.i;
    tauind= (originalG - tauc*C)/(D + I);
    ntauind= update*tauprime.d + (1-update)*tauind;

    taunew = Taxes(ntauind,tauc,ntauind,ntauind,tau.l);
    println("New rates: d = ", taunew.d, " c = ", taunew.c, " i = ", taunew.i, " g = ", taunew.g)

    pr1,eq1=SolveSteadyState(taunew,pa; wguess = wguess, verbose=false);
    if momentsprint
        moments=computeMomentsCutoff(eq1.E,pr1,eq1,tau,pa,cutoffCapital=0.0,toPrint=momentsprint);
    end
    newG=eq1.a.G - eq1.a.collections.l ;;
  end
  println("(originalG - newG)/originalG ",(originalG - newG)/originalG)
  return pr1, eq1, taunew
end


function taxreform6(tauc::Float64, eq::Equilibrium, tau::Taxes, pa::Param; update::Float64 =0.7, tol::Float64 =10.0^-3.0, momentsprint::Bool=false, verbose::Bool=false,firmvalueguess::Matrix = repmat(pa.omega.grid,1,pa.Nz) )

  #Compute tax base for "revenue neutral" reforms
  C = eq.a.collections.c / tau.c; #corporate base
  D = eq.a.collections.d / tau.d; #dividend base
  I = eq.a.collections.i / tau.i; #interest base

  originalG=eq.a.G - 0.7*eq.a.collections.l;
  wguess= eq.w;

  tauind= (originalG - tauc*C - 0.3*eq.a.collections.l)/(D + I);
  taunew = Taxes(tauind,tauc,tauind,tauind,tau.l);
  println("New rates: d = ", taunew.d, " c = ", taunew.c, " i = ", taunew.i, " g = ", taunew.g)


  #Initiate prices and firm problem, and ultimately, the counterfactual object.
  pr1,eq1=SolveSteadyState(taunew,pa; wguess = wguess, verbose=verbose);
  if momentsprint
      moments=computeMomentsCutoff(eq1.E,pr1,eq1,tau,pa,cutoffCapital=0.0,toPrint=momentsprint);
  end
  newG=eq1.a.G- 0.7*eq.a.collections.l;;

  while abs((originalG - newG)/originalG)>tol
    println("(originalG - newG)/originalG ",(originalG - newG)/originalG)
    tauprime=deepcopy(taunew);

    D= eq1.a.collections.d / tauprime.d;
    I= eq1.a.collections.i / tauprime.i;
    tauind= (originalG - tauc*C - 0.3*eq1.a.collections.l)/(D + I);
    ntauind= update*tauprime.d + (1-update)*tauind;

    taunew = Taxes(ntauind,tauc,ntauind,ntauind,tau.l);
    println("New rates: d = ", taunew.d, " c = ", taunew.c, " i = ", taunew.i, " g = ", taunew.g)

    pr1,eq1=SolveSteadyState(taunew,pa; wguess = wguess, verbose=false);
    if momentsprint
        moments=computeMomentsCutoff(eq1.E,pr1,eq1,tau,pa,cutoffCapital=0.0,toPrint=momentsprint);
    end
    newG=eq1.a.G- 0.7*eq.a.collections.l;;
  end
  println("(originalG - newG)/originalG ",(originalG - newG)/originalG)
  return pr1, eq1, taunew
end
