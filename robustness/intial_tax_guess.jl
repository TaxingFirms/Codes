
@everywhere using Grid:CoordInterpGrid, BCnan, BCnearest, InterpLinear
using JLD
using DataFrames

include("markov_approx.jl")
include("mc_tools.jl")

@everywhere include("Main.jl")
@everywhere include("Firms.jl")
@everywhere include("FreeEntry.jl")
@everywhere include("Distribution.jl")
@everywhere include("Aggregation.jl")
@everywhere include("SolveSteadyState.jl")
@everywhere include("calibrate.jl")
@everywhere include("Transitions.jl")



function taxreformRobustness(newTaxGuess::Taxes,tauc::Float64, govexp::Float64, pr::FirmProblem, eq::Equilibrium, tau::Taxes, pa::Param; update::Float64 =0.7, tol::Float64 =10.0^-3.0,taxtol::Float64 =10.0^-6.0, momentsprint::Bool=false,firmvalueguess::Matrix = repmat(pa.omega.grid,1,pa.Nz) )
# Fixes a new level of corporate income taxes, tauc. Closes the deficit using dividend taxes.

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

    taunew = Taxes(newTaxGuess.d,tauc,newTaxGuess.i,newTaxGuess.g,newTaxGuess.l);
    println("New rates: d = ", taunew.d, " c = ", taunew.c, " i = ", taunew.i, " g = ", taunew.g)

    #Initiate prices and firm problem, and ultimately, the counterfactual object.
    pr1,eq1=SolveSteadyState(taunew,pa; wguess = wguess, firmvalueguess = pr.firmvaluegrid, displayit0=false, displayw = false);
    if momentsprint
        moments = computeMomentsCutoff(eq1.E,pr1,eq1,taunew,pa,cutoffCapital=0.0,toPrint=momentsprint);
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

      taunew = Taxes(ntau,tauprime.c,tauprime.i,tauprime.g,tauprime.l);
      println("New rates: d = ", taunew.d, " c = ", taunew.c, " i = ", taunew.i, " g = ", taunew.g)

      initialradius = min(abs(eq1.w-wguess),10.0^-2.0);
      wguess=eq1.w;
      pr1,eq1=SolveSteadyState(taunew,pa; wguess = wguess , firmvalueguess = pr1.firmvaluegrid, displayit0=false, displayw = false , initialradius = initialradius);
      if momentsprint
          moments = computeMomentsCutoff(eq1.E,pr1,eq1,taunew,pa,cutoffCapital=0.0,toPrint=true);
      end
      newG=eq1.a.G;
      deficit= newDeficit;
      newDeficit= (newG - originalG)/originalG;
      taxdif=abs(tau0 - taunew.d); tau0=taunew.d;
    end
    println("pct surplus ",newDeficit)
  return pr1, eq1, taunew
end

  




function compareTaxesToBenchmark(;initial_guess::Taxes = init_taxes(ttaud =0.15, ttauc= 0.35, ttaui= 0.28, ttaug= 0.0, ttaul=0.28),
    
    closing_with_default_guess::Array{Economy,1}=Array(Economy,1),
    tauc_to_try=[0.33 0.31 0.29 0.27 0.25 0.23 0.21 0.19 0.17 0.15 0.13 0.11 0.09 0.07 0.05 0.03 0.01 0.00])
    ~,Nv=size(tauc_to_try);

    pr,eq,tau,pa = load("ModelResults.jld","pr","eq","tau","pa");
    closing_with_new_guess = Array(Economy,length(closing_with_default_guess))
    
    bctol = 5*10.0^-3.0
    update = 0.75
    printMoments = true


    closing_with_new_guess = Array(Economy,(Nv+1,));
    closing_with_new_guess[1] = Economy(pr,eq,tau,0.0);

    for j=1:Nv
        #initialguess = copy(ref[j].pr.firmvaluegrid)
        rpr,req,rtau = taxreformRobustness(initial_guess,tauc_to_try[j], ref[1].eq.a.G, ref[j].pr, ref[j].eq, ref[j].tau, pa; tol=bctol,update=update,momentsprint=printMoments);
        cev= (req.a.consumption - eq.a.consumption)/eq.a.consumption - pa.H/(eq.a.consumption*(1+pa.psi))*( (req.w*(1-rtau.l)/pa.H)^(1+pa.psi) - (eq.w*(1-tau.l)/pa.H)^(1+pa.psi) );
        closing_with_new_guess[j+1]=Economy(rpr,req,rtau,cev);
        save("robustness_check.jld","ref",closing_with_new_guess,"pa",pa);
    end
    
    ref

end


ref,pa = load("reform1ZeroTauI.jld","ref","pa");
newRev = compareTaxesToBenchmark(closing_with_default_guess=ref)

# Print results

robustness_ref,robustness_pa = load("reform1ZeroTauI.jld","ref","pa");
nReformsCompleted = length(robustness_ref)
resultingTaxes = zeros(nReformsCompleted,5)

for i in 1:nReformsCompleted
    resultingTaxes[i,1] = robustness_ref[i].tau.d
    resultingTaxes[i,2] = robustness_ref[i].tau.c
    resultingTaxes[i,3] = robustness_ref[i].tau.i
    resultingTaxes[i,4] = robustness_ref[i].tau.g
    resultingTaxes[i,5] = robustness_ref[i].tau.l
end





