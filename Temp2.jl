function taxreform2(tauc::Float64, eq::Equilibrium, tau::Taxes, pa::Param; update::Float64 =0.7, tol::Float64 =10.0^-2.0, momentsprint::Bool=false, verbose::Bool=false,firmvalueguess::Matrix = repmat(pa.omega.grid,1,pa.Nz) )

#tauc= 0.0; update=0.7; tol= 10.0^-2.0; momentsprint::Bool=false; verbose::Bool=false; initialguess=copy(pr.firmvaluegrid);


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
    pr1,eq1=SolveSteadyState(taunew,pa; wguess = wguess, verbose=verbose, firmvalueguess = firmvalueguess)
    if momentsprint
        moments=computeMomentsCutoff(eq1.E,pr1,eq1,tau,pa,cutoffCapital=0.0,toPrint=momentsprint);
    end

    newG=eq1.a.G;
    newfirmvalueguess=copy(pr1.firmvaluegrid);
    newwguess=eq1.w;

    while abs((originalG - newG)/originalG)>tol
        println("(originalG - newG)/originalG", (originalG - newG)/originalG)
        tau=deepcopy(taunew);

        D= eq1.a.collections.d / tau.d;
        ntau= update*tau.d + (1-update)*(tau.d + (originalG -newG)/D);

        if (1-tau.c)*(1-ntau)>(1-tau.i)
            error("No equilibrium under current taxes")
        end
        taunew = Taxes(ntau,tau.c,tau.i,ntau,tau.l);
        println("New rates: d = ", taunew.d, " c = ", taunew.c, " i = ", taunew.i, " g = ", taunew.g)

        pr1,eq1=SolveSteadyState(taunew,pa; wguess = newwguess, verbose=verbose, firmvalueguess = newfirmvalueguess);
        newG=eq1.a.G;
        newfirmvalueguess=copy(pr1.firmvaluegrid);
        newwguess=eq1.w;
        println("(originalG - newG)/originalG",(originalG - newG)/originalG)
        return pr1, eq1, taunew
    end
end
