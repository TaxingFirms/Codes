

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
  taud =  (govexp - tauhat.c*corpbase - tauhat.i*intbase - eq0.a.collections.g - eq0.a.collections.l)/divbase;
  #2.2 Update tax vector tau
  tau0.d = taud;
  tau0.c = tauhat.c;
  tau0.i = tauhat.i;
  verbose && println("it=0"," New rates: d = ", tau0.d, " c = ", tau0.c, " i = ", tau0.i, " g = ", tau0.g);

  #3 Solve equilibrium for candidate taxes
  if updateVFIguess
    guessVFI = deepcopy(pr0.firmvaluegrid);
  else
    guessVFI = repmat(pa.omega.grid,1,pa.Nz);
  end
  wguess=eq0.w
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


function maximize_welfare()
  # INPUT: economy parameters. Labor and capital gains taxes, wich won't be changed.
  # OUTPUT: welfare maximizing tax vector (taud, tauc, taui) given

  #1. Construct vector of candidate taxes
        ## Note. I think parallelizing at this point is a bad idea as prices should move smoothly
  Ntau=500;
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
  outmaximization=open("sobolmaximization.txt","a")
  println(outmaximization, "-------------------------------------------------------------------------------------------------------")
  flush(outmaximization)
  for j=1:Ntau
    tauhat=Taxes(tau0.d,taxspace[j,2],taxspace[j,1],tau0.g,tau0.l)
    welfarevec[j] = close_gov_bc!(govexp,tauhat,tau0, eq0, pr0, pa; update=0.95, verbose = true, tol=10.0^-3.0)
    #close_gov_vec! updates tau0, eq0, pr0, in place to the new equilibrium under taxes tauhat
    println(outmaximization, taxspace[j,2], " ", welfarevec[j])
    if mod1(j,50)==50
      flush(outmaximization)
    end
  end
  println(outmaximization, "=======================================================================================================")
  close(outmaximization)
end
