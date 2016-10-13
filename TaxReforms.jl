function taxreform1(tauc::Float64, govexp::Float64, pr::FirmProblem, eq::Equilibrium, tau::Taxes, pa::Param; update::Float64 =0.7, tol::Float64 =10.0^-3.0,taxtol::Float64 =10.0^-6.0, momentsprint::Bool=false,firmvalueguess::Matrix = repmat(pa.omega.grid,1,pa.Nz) )
# Gets a new level of tauc. CLoses using the dividend and capital gains taxes at the same time.

#tauc= 0.0; update=0.7; tol= 10.0^-2.0; momentsprint=false; verbose=false; firmvalueguess=copy(pr.firmvaluegrid);


    #Compute tax base for "revenue neutral" reforms
    C = eq.a.collections.c / tau.c #corporate base
    D = eq.a.collections.d / tau.d #divideend base
    I = eq.a.collections.i / tau.i #interest base

    originalG=govexp;
    wguess= eq.w;

    x= (tauc - tau.c)*C/D;

    if (1-tauc)*(1-tau.g)>(1-tau.i)
         error("No equilibrium under current taxes")
    end

    taud= tau.d - (1-update)*x;
    taunew = Taxes(taud,tauc,tau.i,tau.g,tau.l,tau.exit);
    println("New rates: d = ", taunew.d, " c = ", taunew.c, " i = ", taunew.i, " g = ", taunew.g, " exit = ", taunew.exit)

    #Initiate prices and firm problem, and ultimately, the counterfactual object.
    pr1,eq1=SolveSteadyState(taunew,pa; wguess = wguess, firmvalueguess = pr.firmvaluegrid, displayit0=false, displayw = false);
    if momentsprint
        moments=computeMomentsCutoff(eq1.E,pr1,eq1,tau,pa,cutoffCapital=0.0,toPrint=momentsprint);
    end

    newG=eq1.a.G;
    deficit = (newG - originalG)/originalG;
    newDeficit = deficit;
    tau0 = taunew.d; taxdif=Inf;
    while abs(newDeficit)>tol && taxdif > taxtol
      println("pct surplus ",newDeficit)
      tauprime=deepcopy(taunew);

      D= eq1.a.collections.d / tauprime.d;

      if newDeficit*deficit<0 #if deficit changes sign, decrease updating speed
        update+=(1-update)/2;
      end

      ntau= tauprime.d + (1-update)*(originalG -newG)/D;

      if (1-tauprime.c)*(1-tau.g)>(1-tau.i)
          error("No equilibrium under current taxes")
      end

      taunew = Taxes(ntau,tauprime.c,tauprime.i,tauprime.g,tauprime.l,tauprime.exit));
      println("New rates: d = ", taunew.d, " c = ", taunew.c, " i = ", taunew.i, " g = ", taunew.g, " exit = ", taunew.exit)

      initialradius = min(abs(eq1.w-wguess),10.0^-2.0);
      wguess=eq1.w;
      pr1,eq1=SolveSteadyState(taunew,pa; wguess = wguess , firmvalueguess = pr1.firmvaluegrid, displayit0=false, displayw = false , initialradius = initialradius);
      if momentsprint
          moments=computeMomentsCutoff(eq1.E,pr1,eq1,tau,pa,cutoffCapital=0.0,toPrint=true);
      end
      newG=eq1.a.G;
      deficit= newDeficit;
      newDeficit= (newG - originalG)/originalG;
      taxdif=abs(tau0 - taunew.d); tau0=taunew.d;
    end
    println("pct surplus ",newDeficit)
  return pr1, eq1, taunew
end




function taxreform2(tauc::Float64, govexp::Float64, pr::FirmProblem, eq::Equilibrium, tau::Taxes, pa::Param; update::Float64 =0.7, tol::Float64 =10.0^-3.0,taxtol::Float64 =10.0^-6.0, momentsprint::Bool=false,firmvalueguess::Matrix = repmat(pa.omega.grid,1,pa.Nz) )
# Gets a new level of tauc. CLoses using the dividend and capital gains taxes at the same time.

#tauc= 0.0; update=0.7; tol= 10.0^-2.0; momentsprint=false; verbose=false; firmvalueguess=copy(pr.firmvaluegrid);


    #Compute tax base for "revenue neutral" reforms
    C = eq.a.collections.c / tau.c #corporate base
    D = eq.a.collections.d / tau.d #divideend base
    I = eq.a.collections.i / tau.i #interest base

    originalG=govexp;
    wguess= eq.w;

    x= (tauc - tau.c)*C/D;

    if (1-tauc)*(1-(tau.g-x))>(1-tau.i)
         error("No equilibrium under current taxes")
    end

    taue= tau.d - (1-update)*x;
    taunew = Taxes(taue,tauc,tau.i,taue,tau.l,tau.exit);
    println("New rates: d = ", taunew.d, " c = ", taunew.c, " i = ", taunew.i, " g = ", taunew.g, " exit = ", taunew.exit)

    #Initiate prices and firm problem, and ultimately, the counterfactual object.
    pr1,eq1=SolveSteadyState(taunew,pa; wguess = wguess, firmvalueguess = pr.firmvaluegrid, displayit0=false, displayw = false);
    if momentsprint
        moments=computeMomentsCutoff(eq1.E,pr1,eq1,tau,pa,cutoffCapital=0.0,toPrint=momentsprint);
    end

    newG=eq1.a.G;
    deficit = (newG - originalG)/originalG;
    newDeficit = deficit;
    tau0 = taunew.d; taxdif=Inf;
    while abs(newDeficit)>tol && taxdif > taxtol
      println("pct surplus ",newDeficit)
      tauprime=deepcopy(taunew);

      D= eq1.a.collections.d / tauprime.d;

      if newDeficit*deficit<0 #if deficit changes sign, decrease updating speed
        update+=(1-update)/2;
      end

      ntau= tauprime.d + (1-update)*(originalG -newG)/D;

      if (1-tauprime.c)*(1-ntau)>(1-tau.i)
          error("No equilibrium under current taxes")
      end

      taunew = Taxes(ntau,tauprime.c,tauprime.i,ntau,tauprime.l,tauprime.exit);
      println("New rates: d = ", taunew.d, " c = ", taunew.c, " i = ", taunew.i, " g = ", taunew.g, " exit = ", taunew.exit)

      initialradius = min(abs(eq1.w-wguess),10.0^-2.0);
      wguess=eq1.w;
      pr1,eq1=SolveSteadyState(taunew,pa; wguess = wguess , firmvalueguess = pr1.firmvaluegrid, displayit0=false, displayw = false , initialradius = initialradius);
      if momentsprint
          moments=computeMomentsCutoff(eq1.E,pr1,eq1,tau,pa,cutoffCapital=0.0,toPrint=true);
      end
      newG=eq1.a.G;
      deficit= newDeficit;
      newDeficit= (newG - originalG)/originalG;
      taxdif=abs(tau0 - taunew.d); tau0=taunew.d;
    end
    println("pct surplus ",newDeficit)
  return pr1, eq1, taunew
end



function taxreform3(tauc::Float64, govexp::Float64, pr::FirmProblem, eq::Equilibrium, tau::Taxes, pa::Param; update::Float64 =0.7, tol::Float64 =10.0^-3.0,taxtol::Float64 =10.0^-6.0, momentsprint::Bool=false, verbose::Bool=false,firmvalueguess::Matrix = repmat(pa.omega.grid,1,pa.Nz) )
# Gets a new level of tauc. Closes using dividend = capital gains = interest.


  #Compute tax base for "revenue neutral" reforms
  C = eq.a.collections.c / tau.c; #corporate base
  D = eq.a.collections.d / tau.d; #dividend base
  I = eq.a.collections.i / tau.i; #interest base

  originalG=govexp;
  wguess= eq.w;

  tauind= update*(tau.d*D +tau.i*I)/(D+I)+(1-update)*(originalG - tauc*C - eq.a.collections.l)/(D + I);
  taunew = Taxes(tauind,tauc,tauind,tauind,tau.l,tau.exit);
  println("New rates: d = ", taunew.d, " c = ", taunew.c, " i = ", taunew.i, " g = ", taunew.g, " exit = ", taunew.exit)


  #Initiate prices and firm problem, and ultimately, the counterfactual object.
  pr1,eq1=SolveSteadyState(taunew,pa; wguess = wguess, firmvalueguess = pr.firmvaluegrid, displayit0=false, displayw = false);
  if momentsprint
      moments=computeMomentsCutoff(eq1.E,pr1,eq1,tau,pa,cutoffCapital=0.0,toPrint=true);
  end
  newG=eq1.a.G;
  deficit = (newG - originalG)/originalG;
  newDeficit = deficit;
  tau0 = taunew.d; taxdif= Inf;
  while abs(newDeficit)>tol && taxdif > taxtol
    println("pct surplus ",newDeficit)
    tauprime=deepcopy(taunew);

    D= eq1.a.collections.d / tauprime.d;
    I= eq1.a.collections.i / tauprime.i;

    if newDeficit*deficit<0 #if deficit changes sign, decrease updating speed
      update+=(1-update)/2;
    end
    ntau= tauprime.d + (1-update)*(originalG -newG)/(D+I);

    taunew = Taxes(ntau,tauc,ntau,ntau,tau.l,tau.exit);
    println("New rates: d = ", taunew.d, " c = ", taunew.c, " i = ", taunew.i, " g = ", taunew.g, " exit = ", taunew.exit)

    initialradius = min(abs(eq1.w-wguess),10.0^-2.0);
    wguess=eq1.w;
    pr1,eq1=SolveSteadyState(taunew,pa; wguess = wguess , firmvalueguess = pr1.firmvaluegrid, displayit0=false, displayw = false , initialradius = initialradius);
    if momentsprint
        moments=computeMomentsCutoff(eq1.E,pr1,eq1,tau,pa,cutoffCapital=0.0,toPrint=true);
    end
    newG=eq1.a.G;
    deficit= newDeficit;
    newDeficit= (newG - originalG)/originalG;
    taxdif=abs(tau0 - taunew.d); tau0=taunew.d;
  end
  println("pct surplus ",newDeficit)
  return pr1, eq1, taunew
end


function taxreform2_taui(taui::Float64, govexp::Float64, pr::FirmProblem, eq::Equilibrium, tau::Taxes, pa::Param; update::Float64 =0.0, tol::Float64 =10.0^-3.0,taxtol::Float64 =10.0^-6.0, momentsprint::Bool=false,firmvalueguess::Matrix = repmat(pa.omega.grid,1,pa.Nz) )


    #Compute tax base for "revenue neutral" reforms
    C = eq.a.collections.c / tau.c #corporate base
    D = eq.a.collections.d / tau.d #divideend base
    I = eq.a.collections.i / tau.i #interest base

    originalG=govexp;
    wguess= eq.w;

    x= (taui - tau.i)*I/D;

    if (1-tau.c)*(1-(tau.g-x))>(1-taui)
         error("No equilibrium under current taxes")
    end

    taue= tau.d - (1-update)*x;
    taunew = Taxes(taue,tau.c,taui,taue,tau.l);
    println("New rates: d = ", taunew.d, " c = ", taunew.c, " i = ", taunew.i, " g = ", taunew.g)

    #Initiate prices and firm problem, and ultimately, the counterfactual object.
    pr1,eq1=SolveSteadyState(taunew,pa; wguess = wguess, firmvalueguess = pr.firmvaluegrid, displayit0=false, displayw = false);
    if momentsprint
        moments=computeMomentsCutoff(eq1.E,pr1,eq1,tau,pa,cutoffCapital=0.0,toPrint=true);
    end

    newG=eq1.a.G;
    deficit = (newG - originalG)/originalG;
    newDeficit = deficit;
    tau0 = taunew.d; taxdif= Inf;
    while abs(newDeficit)>tol && taxdif > taxtol
      println("pct surplus ",newDeficit)
      tauprime=deepcopy(taunew);

      D= eq1.a.collections.d / tauprime.d;
      if newDeficit*deficit<0.0  #if deficit changes sign, decrease updating speed
        update+=(1-update)/2;
      end

      ntau= tauprime.d + (1-update)*(originalG -newG)/D;

      if (1-tauprime.c)*(1-ntau)>(1-tau.i)
          error("No equilibrium under current taxes")
      end
      taunew = Taxes(ntau,tauprime.c,tauprime.i,ntau,tau.l);
      println("New rates: d = ", taunew.d, " c = ", taunew.c, " i = ", taunew.i, " g = ", taunew.g)

      initialradius = min(abs(eq1.w-wguess),10.0^-2.0);
      wguess=eq1.w;
      pr1,eq1=SolveSteadyState(taunew,pa; wguess = wguess , firmvalueguess = pr1.firmvaluegrid, displayit0=false, displayw = false , initialradius = initialradius);
      if momentsprint
          moments=computeMomentsCutoff(eq1.E,pr1,eq1,tau,pa,cutoffCapital=0.0,toPrint=momentsprint);
      end
      newG=eq1.a.G;
      deficit= newDeficit;
      newDeficit= (newG - originalG)/originalG;
      taxdif=abs(tau0 - taunew.d); tau0=taunew.d;
    end
    println("pct surplus ",newDeficit)
  return pr1, eq1, taunew
end



function taxreform4(tauc::Float64, govexp::Float64, pr::FirmProblem, eq::Equilibrium, tau::Taxes, pa::Param; update::Float64 =0.7, tol::Float64 =10.0^-3.0,taxtol::Float64 =10.0^-6.0, momentsprint::Bool=false, verbose::Bool=false,firmvalueguess::Matrix = repmat(pa.omega.grid,1,pa.Nz) )
# Gets a new level of tauc. Closes using dividend = capital gains = interest.


  #Compute tax base for "revenue neutral" reforms
  C = eq.a.collections.c / tau.c; #corporate base
  D = eq.a.collections.d / tau.d; #dividend base
  I = eq.a.collections.i / tau.i; #interest base

  originalG=govexp;
  wguess= eq.w;

  tauind= update*(tau.d*D+tau.i*I)/(D+I) +(1-update)*(originalG - tauc*C - eq.a.collections.l)/(D + I);
  taunew = Taxes(tauind,tauc,tauind,tau.g,tau.l);
  println("New rates: d = ", taunew.d, " c = ", taunew.c, " i = ", taunew.i, " g = ", taunew.g)


  #Initiate prices and firm problem, and ultimately, the counterfactual object.
  pr1,eq1=SolveSteadyState(taunew,pa; wguess = wguess, firmvalueguess = pr.firmvaluegrid, displayit0=false, displayw = false);
  if momentsprint
      moments=computeMomentsCutoff(eq1.E,pr1,eq1,tau,pa,cutoffCapital=0.0,toPrint=true);
  end
  newG=eq1.a.G;
  deficit = (newG - originalG)/originalG;
  newDeficit = deficit;
  tau0 = taunew.d; taxdif= Inf;
  while abs(newDeficit)>tol && taxdif > taxtol
    println("pct surplus ",newDeficit)
    tauprime=deepcopy(taunew);

    D= eq1.a.collections.d / tauprime.d;
    I= eq1.a.collections.i / tauprime.i;

    if newDeficit*deficit<0 #if deficit changes sign, decrease updating speed
      update+=(1-update)/2;
    end
    ntau= tauprime.d + (1-update)*(originalG -newG)/(D+I);

    taunew = Taxes(ntau,tauc,ntau,tau.g,tau.l);
    println("New rates: d = ", taunew.d, " c = ", taunew.c, " i = ", taunew.i, " g = ", taunew.g)

    initialradius = min(abs(eq1.w-wguess),10.0^-2.0);
    wguess=eq1.w;
    pr1,eq1=SolveSteadyState(taunew,pa; wguess = wguess , firmvalueguess = pr1.firmvaluegrid, displayit0=false, displayw = false , initialradius = initialradius);
    if momentsprint
        moments=computeMomentsCutoff(eq1.E,pr1,eq1,tau,pa,cutoffCapital=0.0,toPrint=true);
    end
    newG=eq1.a.G;
    deficit= newDeficit;
    newDeficit= (newG - originalG)/originalG;
    taxdif=abs(tau0 - taunew.d); tau0=taunew.d;
  end
  println("pct surplus ",newDeficit)
  return pr1, eq1, taunew
end


function taxreform5(tauc::Float64, govexp::Float64, pr::FirmProblem, eq::Equilibrium, tau::Taxes, pa::Param; update::Float64 =0.7, tol::Float64 =10.0^-3.0,taxtol::Float64 =10.0^-6.0, momentsprint::Bool=false, verbose::Bool=false,firmvalueguess::Matrix = repmat(pa.omega.grid,1,pa.Nz) )
# Sets taui to zero and decreases tauc by increaseing taud

  #Compute tax base for "revenue neutral" reforms
  C = eq.a.collections.c / tau.c; #corporate base
  D = eq.a.collections.d / tau.d; #dividend base
  I = eq.a.collections.i / tau.i; #interest base

  originalG=govexp;
  wguess= eq.w;

  tauind=  tau.d; #update*tau.d + (1-update)*(originalG - tauc*C - eq.a.collections.l)/D;
  taunew = Taxes(tauind,tauc,0.0,tau.g,tau.l,tau.exit);
  println("New rates: d = ", taunew.d, " c = ", taunew.c, " i = ", taunew.i, " g = ", taunew.g, " exit = ", taunew.exit )


  #Initiate prices and firm problem, and ultimately, the counterfactual object.
  pr1,eq1=SolveSteadyState(taunew,pa; wguess = wguess, firmvalueguess = pr.firmvaluegrid, displayit0=false, displayw = false);
  println(@sprintf(" Dividend Base \t DivsPrevDistr \t  Consumption \t Mass Entrants \t      Wage     \t     Welfare   \t      TFP     \t       G     "))
  println(@sprintf("%9.4f \t %9.4f \t %9.4f \t %9.4f \t %9.4f  \t %9.4f \t %9.4f  \t %9.4f ", eq1.a.collections.d/taunew.d, sum(eq.distr.*pr1.grossdividends), eq1.a.consumption, eq1.E, eq1.w, eq1.a.welfare , eq1.a.output/(eq1.a.capital^pa.alphak*eq1.a.laborsupply^pa.alphal), eq1.a.G))
  if momentsprint
      moments=computeMomentsCutoff(eq1.E,pr1,eq1,tau,pa,cutoffCapital=0.0,toPrint=true);
  end
  newG=eq1.a.G;
  deficit = (newG - originalG)/originalG;
  newDeficit = deficit;
  tau0 = taunew.d; taxdif= Inf;
  while abs(newDeficit)>tol && taxdif > taxtol
    println("pct surplus ",newDeficit)
    tauprime=deepcopy(taunew);

    D= eq1.a.collections.d / tauprime.d;

    if newDeficit*deficit<0 #if deficit changes sign, decrease updating speed
      update+=(1-update)/2;
    end
    ntau= tauprime.d + (1-update)*(originalG -newG)/D;

    taunew = Taxes(ntau,tauc,0.0,tau.g,tau.l,tau.exit);
    println("New rates: d = ", taunew.d, " c = ", taunew.c, " i = ", taunew.i, " g = ", taunew.g, " exit = ", taunew.exit )

    initialradius = min(abs(eq1.w-wguess),10.0^-2.0);
    wguess=eq1.w;
    pr1,eq1=SolveSteadyState(taunew,pa; wguess = wguess , firmvalueguess = pr1.firmvaluegrid, displayit0=false, displayw = false , initialradius = initialradius);
    println(@sprintf(" Dividend Base \t DivsPrevDistr \t  Consumption \t Mass Entrants \t      Wage     \t     Welfare   \t      TFP     \t       G     "))
    println(@sprintf("%9.4f \t %9.4f \t %9.4f \t %9.4f \t %9.4f  \t %9.4f \t %9.4f  \t %9.4f ", eq1.a.collections.d/taunew.d,  sum(eq.distr.*pr1.grossdividends), eq1.a.consumption, eq1.E, eq1.w, eq1.a.welfare , eq1.a.output/(eq1.a.capital^pa.alphak*eq1.a.laborsupply^pa.alphal), eq1.a.G))

    if momentsprint
        moments=computeMomentsCutoff(eq1.E,pr1,eq1,tau,pa,cutoffCapital=0.0,toPrint=true);
    end
    newG=eq1.a.G;
    deficit= newDeficit;
    newDeficit= (newG - originalG)/originalG;
    taxdif=abs(tau0 - taunew.d); tau0=taunew.d;
  end
  println("pct surplus ",newDeficit)

  return pr1, eq1, taunew
end



function taxreform6(tauc::Float64, govexp::Float64, pr::FirmProblem, eq::Equilibrium, tau::Taxes, pa::Param; update::Float64 =0.7, tol::Float64 =10.0^-3.0,taxtol::Float64 =10.0^-6.0, momentsprint::Bool=false, verbose::Bool=false,firmvalueguess::Matrix = repmat(pa.omega.grid,1,pa.Nz) )
# Sets taui = tauc and decreases tauc by increaseing taud

  #Compute tax base for "revenue neutral" reforms
  C = eq.a.collections.c / tau.c; #corporate base
  D = eq.a.collections.d / tau.d; #dividend base
  I = eq.a.collections.i / tau.i; #interest base

  originalG=govexp;
  wguess= eq.w;

  tauind= update*tau.d + (1-update)*(originalG - tauc*(C+I) - eq.a.collections.l)/D;
  taunew = Taxes(tauind,tauc,tauc-eps(),tau.g,tau.l,tau.exit);
  println("New rates: d = ", taunew.d, " c = ", taunew.c, " i = ", taunew.i, " g = ", taunew.g)


  #Initiate prices and firm problem, and ultimately, the counterfactual object.
  pr1,eq1=SolveSteadyState(taunew,pa; wguess = wguess, firmvalueguess = pr.firmvaluegrid, displayit0=false, displayw = false);
  if momentsprint
      moments=computeMomentsCutoff(eq1.E,pr1,eq1,tau,pa,cutoffCapital=0.0,toPrint=true);
  end
  newG=eq1.a.G;
  deficit = (newG - originalG)/originalG;
  newDeficit = deficit;
  tau0 = taunew.d; taxdif= Inf;
  while abs(newDeficit)>tol && taxdif > taxtol
    println("pct surplus ",newDeficit)
    tauprime=deepcopy(taunew);

    D= eq1.a.collections.d / tauprime.d;

    if newDeficit*deficit<0 #if deficit changes sign, decrease updating speed
      update+=(1-update)/2;
    end
    ntau= tauprime.d + (1-update)*(originalG -newG)/D;

    taunew = Taxes(ntau,tauc,tauc-eps(),tau.g,tau.l,tau.exit);
    println("New rates: d = ", taunew.d, " c = ", taunew.c, " i = ", taunew.i, " g = ", taunew.g)

    initialradius = min(abs(eq1.w-wguess),10.0^-2.0);
    wguess=eq1.w;
    pr1,eq1=SolveSteadyState(taunew,pa; wguess = wguess , firmvalueguess = pr1.firmvaluegrid, displayit0=false, displayw = false , initialradius = initialradius);
    if momentsprint
        moments=computeMomentsCutoff(eq1.E,pr1,eq1,tau,pa,cutoffCapital=0.0,toPrint=true);
    end
    newG=eq1.a.G;
    deficit= newDeficit;
    newDeficit= (newG - originalG)/originalG;
    taxdif=abs(tau0 - taunew.d); tau0=taunew.d;
  end
  println("pct surplus ",newDeficit)
  return pr1, eq1, taunew
end


function taxreform1old(tauc::Float64,  govexp::Float64, pr::FirmProblem, eq::Equilibrium, tau::Taxes, pa::Param; update::Float64 =0.7, tol::Float64 =10.0^-3.0,taxtol::Float64 =10.0^-6.0, momentsprint::Bool=false, verbose::Bool=false,firmvalueguess::Matrix = repmat(pa.omega.grid,1,pa.Nz))
# Gets a new level of tauc. CLoses using the dividend tax
#Compute tax base for "revenue neutral" reforms
C = eq.a.collections.c / tau.c #corporate base
D = eq.a.collections.d / tau.d #divideend base
I = eq.a.collections.i / tau.i #interest base

originalG=govexp;
wguess= eq.w;

x= (tauc - tau.c)*C/D;

if (1-tauc)*(1-(tau.g))>(1-tau.i)
     error("No equilibrium under current taxes")
end

taud=tau.d - (1-update)*x
taunew = Taxes(taud,tauc,tau.i,tau.g,tau.l);
println("New rates: d = ", taunew.d, " c = ", taunew.c, " i = ", taunew.i, " g = ", taunew.g)

#Initiate prices and firm problem, and ultimately, the counterfactual object.
pr1,eq1=SolveSteadyState(taunew,pa; wguess = wguess, firmvalueguess = pr.firmvaluegrid, displayit0=false, displayw = false);
if momentsprint
    moments=computeMomentsCutoff(eq1.E,pr1,eq1,tau,pa,cutoffCapital=0.0,toPrint=momentsprint);
end

newG=eq1.a.G;
deficit = (newG - originalG)/originalG;
newDeficit = deficit;
tau0 = taunew.d; taxdif=Inf;
while abs(newDeficit)>tol && taxdif > taxtol
  println("pct surplus ",newDeficit)
  tauprime=deepcopy(taunew);

  D= eq1.a.collections.d / tauprime.d;

  if newDeficit*deficit<0 #if deficit changes sign, decrease updating speed
    update+=(1-update)/2;
  end

  ntau= tauprime.d + (1-update)*(originalG -newG)/D;

  taunew = Taxes(ntau,tauprime.c,tauprime.i,tauprime.g,tauprime.l);
  println("New rates: d = ", taunew.d, " c = ", taunew.c, " i = ", taunew.i, " g = ", taunew.g)

  initialradius = min(abs(eq1.w-wguess),10.0^-2.0);
  wguess=eq1.w;
  pr1,eq1=SolveSteadyState(taunew,pa; wguess = wguess , firmvalueguess = pr1.firmvaluegrid, displayit0=false, displayw = false , initialradius = initialradius);
  if momentsprint
      moments=computeMomentsCutoff(eq1.E,pr1,eq1,tau,pa,cutoffCapital=0.0,toPrint=true);
  end
  newG=eq1.a.G;
  deficit= newDeficit;
  newDeficit= (newG - originalG)/originalG;
  taxdif=abs(tau0 - taunew.d); tau0=taunew.d;
end
println("pct surplus ",newDeficit)
return pr1, eq1, taunew
end
