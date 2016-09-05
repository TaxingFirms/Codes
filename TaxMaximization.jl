

function close_gov_bc!(tauhat::Taxes,tau::Taxes, eq::Equilibrium, pr:: FirmProblem, pa::Param; update::Float64 =0.0, verbose = true, tol::Float64 =10.0^-3.0, returnall::Bool = false, updateVFIguess::Bool = true)
  # INPUT: rates tauc and taui in tauhat, G in eq.a.G and an initial equilibrium guess.
  # It MODIFIES taud such that the government budget constraint is satisfied
  # OUTPUT steady state welfare

  #0. Check there is an equilibrium under current taxes
  (1-tauhat.c)*(1-tauhat.g)>(1-tauhat.i) && return -Inf

  #1. Save government expenditure target
  govexp = eq.a.G;

  #2 Guess new tax vector
  #2.0 Compute tax bases for taxes that will change
  corpbase = eq.a.collections.c / tau.c;
  divbase = eq.a.collections.d / tau.d;
  intbase = eq.a.collections.i / tau.i
  #2.1 Guess taud
  taud =  update*tau.d + (1-update)*(govexp - tauhat.c*corpbase - tauhat.i*intbase - eq.a.collections.g - eq.a.collections.l)/divbase;
  #2.2 Generate new tax vector tau0
  tau0 = Taxes(taud,tauhat.c,tauhat.i,tau.g,tau.l);
  verbose && println("it=0"," New rates: d = ", tau0.d, " c = ", tau0.c, " i = ", tau0.i, " g = ", tau0.g);

  #3 Solve equilibrium for candidate taxes
  guessVFI = repmat(pa.omega.grid,1,pa.Nz);
  if updateVFIguess
    guessVFI = deepcopy(pr.firmvaluegrid)
  end
  pr0, eq0 = SolveSteadyState(tau0, pa; wguess = eq.w, VFItol =10.0^-3.0,  VFIfunction = firmVFIParallelOmega!, displayit0 = true, displayw = true, firmvalueguess = guessVFI);
  govrevenue = eq0.a.G;

  #4 While government revenue is different from its expenditures, keep iterating
  itcount=1; maxit = 500; relgap = (govrevenue - govexp)/govexp;
  verbose && println(" rev - exp = ", relgap);
  while abs(relgap)>tol && itcount < maxit;
    #4.1 Update taxes
    corpbase = eq.a.collections.c / tau.c; divbase = eq.a.collections.d / tau.d; intbase = eq.a.collections.i / tau.i;
    taud = update*tau.d + (1-update)*(govexp - tau0.c*corpbase - tau0.i*intbase - eq0.a.collections.g - eq0.a.collections.l)/divbase;
    tau0 = Taxes(taud,tau0.c,tau0.i,tau0.g,tau0.l);
    verbose && println("it=",itcount," New rates: d = ", tau0.d, " c = ", tau0.c, " i = ", tau0.i, " g = ", tau0.g);
    #4.2 Solve equilibrium
    if updateVFIguess
      guessVFI = deepcopy(pr0.firmvaluegrid)
    end
    pr0, eq0 = SolveSteadyState(tau0, pa; wguess = eq0.w, VFItol =10.0^-3.0,  VFIfunction = firmVFIParallelOmega!, displayit0 =true, displayw = true, firmvalueguess = guessVFI);
    relgap = (eq0.a.G - govexp)/govexp;
    verbose && println(" rev - exp = ", relgap);
    #4.3 Update government revenue
    govrevenue = eq0.a.G; itcount +=1;
  end

  #5.1 The function allows to return the entire equilibrium
  returnall && return eq0.a.welfare, pr0, eq0, tau0
  #5.2 Default is to just return welafare
  eq0.a.welfare
end


function maximize_welfare(tau_g::Float64, tau_l::Float64, pa::Param)
  # INPUT: economy parameters. Labor and capital gains taxes, wich won't be changed.
  # OUTPUT: welfare maximizing tax vector (taud, tauc, taui) given

  #1. Construct vector of candidate taxes
        ## Note. I think parallelizing at this point is a bad idea as prices should move smoothly
  Ntau=1000;
  #1.1 Create Sobol Matrix
  s=SobolSeq(2); skip(s,Ntau);
  tauivec = Array()
  #1.2 Sort array taui, array tauc
  #1.3 Create algorithm to move through the space continuously


  #2. Initialize vector of objects (tau, eq, pr, pa) where eq,pr,pa are initialized but empty


  #3. Initial equilibirium: under current taxes. This fixes the value of G


  #4. Loop over taxes moving them slowly to neighboring combinations


end



function maximize_welfare_parallel(tau_g::Float64, tau_l::Float64, pa::Param)
    # INPUT: economy parameters. Labor and capital gains taxes, wich won't be changed.
    # OUTPUT: welfare maximizing tax vector (taud, tauc, taui) given

    #1. Construct vector of candidate taxes (perhaps using Sobol)


    #2. Initialize vector of objects (tau, eq, pr, pa) where eq,pr,pa are initialized but empty


    #3. Initial equilibirium: under current taxes. This fixes the value of G
        #This will be the COMMON initial guess for all iterations


    #4. Pmap cadidate taxes @tau close_gov_bc(tau,eq,pr,pa)

end
