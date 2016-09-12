using Sobol


function close_gov_tauc!(govexp::Float64,tau0::Taxes, pa::Param; update::Float64 =0.9, verbose = true, tol::Float64 =10.0^-3.0, returnall::Bool = false, wguess::Float64 = 0.55, updateVFIguess::Bool =true, outsideparallel::Bool =true)
  # INPUT: Taxes on dividends and interests + economy constants including government expenditure target, and tau.l, tau.g in tauhat.
  # It MODIFIES tauc in placesuch that the government budget constraint is satisfied
  # OUTPUT steady state welfare

  #0. Set level of parallelization
  if outsideparallel
    VFIfunc = firmVFI!
  else
    VFIfunc = firmVFIParallelOmega!
  end

  #1. Start at the lower bound for tauc
  tau0.c = max(1 - (1-tau0.i)/(1-tau0.g) , 0)

  #2. Compute Equilibirum under current taxes
  verbose && println("it=0"," New rates: d = ", tau0.d, " c = ", tau0.c, " i = ", tau0.i, " g = ", tau0.g);
  pr0, eq0 = SolveSteadyState(tau0, pa; wguess = wguess, VFItol =10.0^-3.0,  VFIfunction = VFIfunc, displayit0 = false, displayw = false);

  #3. If goverment contraint is satisfied we are done
  verbose && println("revenue = ", eq0.a.G, " expenditure = ", govexp);
  eq0.a.G >= govexp && return eq0.a.welfare

  #4. Keep track of how the revenue is changing
    flag=0; previousrevenue=eq0.a.G;

  #5. While the government constraint is not satisfied,
  itcount=1; maxit = 500; relgap = (eq0.a.G- govexp)/govexp;
  verbose && println("Relative budget deficit = ", relgap);
  while abs(relgap)>tol && itcount < maxit;
    #5.1 Compute tax bases for taxes that will change
    corpbase = eq0.a.collections.c / tau0.c;
    divbase = eq0.a.collections.d / tau0.d;
    intbase = eq0.a.collections.i / tau0.i
    #5.2 Update tauc
    tau0.c =  update*tau0.c + (1-update)*(govexp - tau0.d*divbase - tau0.i*intbase - eq0.a.collections.g - eq0.a.collections.l)/corpbase;
    verbose && println("it= ",itcount," New rates: d = ", tau0.d, " c = ", tau0.c, " i = ", tau0.i, " g = ", tau0.g);

    #5.3 Solve equilibrium for candidate taxes
    if updateVFIguess
      guessVFI = deepcopy(pr0.firmvaluegrid);
    else
      guessVFI = repmat(pa.omega.grid,1,pa.Nz);
    end
    initialradius = min(abs(eq0.w-wguess),10.0^-2.0); #How far to look for wages
    wguess=eq0.w
    pr0, eq0 = SolveSteadyState(tau0, pa; wguess = wguess, VFItol =10.0^-3.0,  VFIfunction = VFIfunc, displayit0 = false, displayw = false, firmvalueguess = guessVFI, initialradius = initialradius);

    #5.4 Update bugdet deficit
    relgap = (eq0.a.G - govexp)/govexp;
    verbose && println("it= ", itcount, " Relative budget deficit = ", relgap);
    #5.5 Update iteration counter
    itcount +=1;
    eq0.a.G< previousrevenue && flag+=1;
    flag == 3 && return -Inf
  end

  #6.1 The function allows to return the entire equilibrium
  returnall && return eq0.a.welfare, pr0, eq0
  #6.2 Default is to just return welafare
  eq0.a.welfare
end



function close_gov_bc(govexp::Float64,tauhat::Taxes,tau::Taxes, eq::Equilibrium, pr:: FirmProblem, pa::Param; update::Float64 =0.9, verbose = true, tol::Float64 =10.0^-3.0, returnall::Bool = false, updateVFIguess::Bool = true)
  # INPUT: rates tauc and taui in tauhat, G in govexp and an initial equilibrium guess.
  # It MODIFIES taud such that the government budget constraint is satisfied
  # OUTPUT steady state welfare

  #1. Check there is an equilibrium under current taxes
  (1-tauhat.c)*(1-tauhat.g)>(1-tauhat.i) && return -Inf

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
  if updateVFIguess
    guessVFI = deepcopy(pr.firmvaluegrid);
  else
    guessVFI = repmat(pa.omega.grid,1,pa.Nz);
  end
  wguess=eq.w
  pr0, eq0 = SolveSteadyState(tau0, pa; wguess = wguess, VFItol =10.0^-3.0,  VFIfunction = firmVFIParallelOmega!, displayit0 = false, displayw = false, firmvalueguess = guessVFI);
  govrevenue = eq0.a.G;

  #4 While government revenue is different from its expenditures, keep iterating
  itcount=1; maxit = 500; relgap = (govrevenue - govexp)/govexp;
  verbose && println(" rev - exp = ", relgap);
  while abs(relgap)>tol && itcount < maxit;
    #4.1 Update taxes
    divbase = eq0.a.collections.d / tau0.d;
    taud = update*tau0.d + (1-update)*(govexp - eq0.a.collections.c - eq0.a.collections.i - eq0.a.collections.g - eq0.a.collections.l)/divbase;
    tau0 = Taxes(taud,tau0.c,tau0.i,tau0.g,tau0.l);
    verbose && println("it=",itcount," New rates: d = ", tau0.d, " c = ", tau0.c, " i = ", tau0.i, " g = ", tau0.g);
    #4.2 Solve equilibrium
    if updateVFIguess
      guessVFI = deepcopy(pr0.firmvaluegrid)
    end
    initialradius = min(abs(eq0.w-wguess),10.0^-2.0); #How far to look for wages
    wguess=eq0.w;
    pr0, eq0 = SolveSteadyState(tau0, pa; wguess = eq0.w, VFItol =10.0^-3.0,  VFIfunction = firmVFIParallelOmega!, displayit0 =false, displayw = false, firmvalueguess = guessVFI, initialradius=initialradius);
    relgap = (eq0.a.G - govexp)/govexp;
    verbose && println(" rev - exp = ", relgap);
    println(govexp," ", eq0.a.collections)
    #4.3 Update government revenue
    govrevenue = eq0.a.G; itcount +=1;
  end

  #5.1 The function allows to return the entire equilibrium
  returnall && return eq0.a.welfare, pr0, eq0, tau0
  #5.2 Default is to just return welafare
  eq0.a.welfare
end


function sortsobol(initialpoint::Array{Float64,1}, space = Array{Float64,1})
  Ntau,dimtau = size(space)
  # Initialize Array sorted space
  sortedspace = Array(Float64,(Ntau,dimtau));
  distances = Array(Float64,(Ntau))
  #Loop over each point in the space
  for j=1:Ntau
    #1. Compute a vector of distances to initial point
    for k=1:Ntau
      distances[k] =  ((space[k,1] - initialpoint[1])^2.0 + (space[k,2] - initialpoint[2])^2.0 )^0.5
    end
    #2. Minimize the distances
    closest = indmin(distances)
    #3. Save closest point in the sorted Array
    sortedspace[j,:]= space[closest,:]
    #3. Update initial point to new point
    initialpoint = space[closest,:]
    #4. Set that point to infinity so it is not taken again
    space[closest,:]= [Inf Inf]
  end
  sortedspace
end

function maximize_welfare_parallel(tau_g::Float64, tau_l::Float64, pa::Param)
    # INPUT: economy parameters. Labor and capital gains taxes, wich won't be changed.
    # OUTPUT: welfare maximizing tax vector (taud, tauc, taui) given

    #1. Construct vector of candidate taxes (perhaps using Sobol)
    Ntau=500;
    dimtau=2;
    #1.1 Create tax vectors: first component is taui, second is tauc
    s=SobolSeq(dimtau, [0.0 0.0], [0.5 0.5]); skip(s,Ntau);
    sobolspace = Array(Float64,(Ntau,dimtau));
    for j=1:Ntau
      sobolspace[j,:]= next(s);
    end

    #2. Initialize vector of objects (tau, eq, pr, pa) where eq,pr,pa are initialized but empty
    taxvec = Array(Taxes,(Ntau,));
    for i=1:Ntau
      taxvec[i]= Taxes(sobolspace[1],0.0,sobolspace[2],tau.g,tau.l)
    end
    welfarevec = Array(Float64,(Ntau,));

    #3. Initial equilibirium: under current taxes. This fixes the value of G
        #This will be the COMMON initial guess for all iterations
    pr,eq,tau,pa=load("ModelResults.jld", "pr","eq","tau","pa");
    govexp=eq.a.G;

    #4. Pmap cadidate taxes @tau close_gov_bc(tau,eq,pr,pa)
    map(objectivefun, input)

end

#objectivefun(arg)
#  close_gov_tauc!(arg.govexp,arg.tau, arg.pa; update=0.9, verbose = false, tol =10.0^-3.0, wguess = 0.55)
#end
