using Sobol

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
  Ntau=100;
  dimtau=2;
  #1.1 Create tax vectors: first component is taui, second is tauc
  s=SobolSeq(dimtau); skip(s,Ntau);
  sobolspace = Array(Float64,(Ntau,dimtau));
  for j=1:Ntau
    sobolspace[j,:]= next(s);
  end
  #1.3 Create algorithm to move through the space continuously
  initialpoint = [0.28,0.35];
  taxspace =sortsobol(initialpoint, sobolspace);

  #2. Initialize vector of objects (tau, eq, pr, pa) where eq,pr,pa are initialized but empty
  welfarevec = Array(Float64,(Ntau,));

  #3. Initial equilibirium: under current taxes. This fixes the value of G
  pr,eq,tau,pa=load("ModelResults.jld", "pr","eq","tau","pa");
  pr0=deepcopy(pr); eq0=deepcopy(eq); tau0=deepcopy(tau);
  govexp = eq.a.G;

  #4. Loop over taxes moving them slowly to neighboring combinations
  for j=1:Ntau
    tauhat=Taxes(tau0.d,taxspace[j,2],taxspace[j,1],tau0.g,tau0.l)
    welfarevec[j] = close_gov_bc!(govexp,tauhat,tau0, eq0, pr0, pa; update=0.9, verbose = true, tol=10.0^-3.0)
    #close_gov_vec! updates tau0, eq0, pr0, in place to the new equilibrium under taxes tauhat
    println(taxspace[j,2], " ", welfarevec[j])
  end


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


    #2. Initialize vector of objects (tau, eq, pr, pa) where eq,pr,pa are initialized but empty


    #3. Initial equilibirium: under current taxes. This fixes the value of G
        #This will be the COMMON initial guess for all iterations


    #4. Pmap cadidate taxes @tau close_gov_bc(tau,eq,pr,pa)

end

function close_gov_bc!(govexp::Float64,tauhat::Taxes,tau0::Taxes, eq0::Equilibrium, pr0:: FirmProblem, pa::Param; update::Float64 =0.9, verbose = true, tol::Float64 =10.0^-3.0, returnall::Bool = false, updateVFIguess::Bool = true)
  # INPUT: rates tauc and taui in tauhat, G in govexp and an initial equilibrium guess.
  # It MODIFIES tau0, eq0, pr0 such that the government budget constraint is satisfied under the new taxes tauhat (now in tau)
  # OUTPUT steady state welfare

  #1. Check there is an equilibrium under current taxes
  (1-tauhat.c)*(1-tauhat.g)>(1-tauhat.i) && return -Inf

  #2 Guess new tax vector
  #2.0 Compute tax bases for taxes that will change
  corpbase = eq0.a.collections.c / tau0.c;
  divbase = eq0.a.collections.d / tau0.d;
  intbase = eq0.a.collections.i / tau0.i
  #2.1 Guess taud
  taud =  update*tau0.d + (1-update)*(govexp - tauhat.c*corpbase - tauhat.i*intbase - eq0.a.collections.g - eq0.a.collections.l)/divbase;
  #2.2 Update tax vector tau
  tau0.d = taud;
  tau0.c = tauhat.c;
  tau0.i = tauhat.i;
  verbose && println("it=0"," New rates: d = ", tau0.d, " c = ", tau0.c, " i = ", tau0.i, " g = ", tau0.g);

  #3 Solve equilibrium for candidate taxes
  if updateVFIguess
    guessVFI = deepcopy(pr.firmvaluegrid);
  else
    guessVFI = repmat(pa.omega.grid,1,pa.Nz);
  end
  wguess=eq.w
   SolveSteadyState!( eq0,pr0, tau0, pa; wguess = wguess, VFItol =10.0^-3.0,  VFIfunction = firmVFIParallelOmega!, displayit0 = false, displayw = false, firmvalueguess = guessVFI);
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
    pr0, eq0 = SolveSteadyState!(eq0, pr0, tau0, pa; wguess = eq0.w, VFItol =10.0^-3.0,  VFIfunction = firmVFIParallelOmega!, displayit0 =false, displayw = false, firmvalueguess = guessVFI, initialradius=initialradius);
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
