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





function maximize_welfare_alltaxes(initialpoint::Array{Float64})
  # INPUT: economy parameters. Labor and capital gains taxes, wich won't be changed.
  # OUTPUT: welfare maximizing tax vector (taud, tauc, taui) given

  #1. Construct vector of candidate taxes
        ## Note. I think parallelizing at this point is a bad idea as prices should move smoothly
  Ntau=4500;
  dimtau=3;
  #1.1 Create tax vectors: first component is taud, second is taui, third is taug
  s=SobolSeq(dimtau, [0.0 0.0 0.0], [0.6 0.6 0.6]); skip(s,Ntau);
  sobolspace = Array(Float64,(Ntau,dimtau));
  for j=1:Ntau
    sobolspace[j,:]= next(s);
  end
  #1.3 Create algorithm to move through the space continuously
  taxspace =sortsobol(initialpoint, sobolspace);

  #2. Initialize vector of objects (tau, eq, pr, pa) where eq,pr,pa are initialized but empty
  welfarevec = Array(Float64,(Ntau,));

  #3. Initial equilibirium: under current taxes. This fixes the value of G

  #4. Loop over taxes moving them slowly to neighboring combinations
end





function maximize_welfare_parallel(tau_g::Float64, tau_l::Float64, pa::Param)
    # INPUT: economy parameters. Labor and capital gains taxes, wich won't be changed.
    # OUTPUT: welfare maximizing tax vector (taud, tauc, taui) given

    #1. Construct vector of candidate taxes (perhaps using Sobol)
    Ntau=1200;
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


function sortsobol(initialpoint::Array{Float64,1}, space = Array{Float64,2})
  Ntau,dimtau = size(space)
  # Initialize Array sorted space
  sortedspace = Array(Float64,(Ntau,dimtau));
  distances = Array(Float64,(Ntau,))
  #Loop over each point in the space
  for j=1:Ntau
    #1. Compute a vector of distances to initial point
    for k=1:Ntau
      sum = 0.0;
      for l=1:dimtau
        sum+= (space[k,l] - initialpoint[l])^2.0
      end
      distances[k] = sum^0.5
    end
    #2. Minimize the distances
    closest = indmin(distances)
    #3. Save closest point in the sorted Array
    sortedspace[j,:]= space[closest,:]
    #3. Update initial point to new point
    initialpoint = space[closest,:]
    #4. Set that point to infinity so it is not taken again
    space[closest,:]= Inf*ones(1,dimtau)
  end
  sortedspace
end


function create_taxspace(initialpoint::Array{Float64,1}, Ntau::Int64; jldfile::Bool=true)
  # INPUT: Initial point.
  # OUTPUT: JLD file with a sequence o

  dimtau=size(initialpoint)[1];
  s=SobolSeq(dimtau, zeros(dimtau), ones(dimtau)); skip(s,Ntau);
  sobolspace = Array(Float64,(Ntau,dimtau));
  for j=1:Ntau
    sobolspace[j,:]= next(s);
  end
  taxspace =sortsobol(initialpoint, sobolspace);
  if jldfile
    save("taxspace.jld","taxspace",taxspace)
  end

  return taxspace
end

######################################################################################################################

function create_half_taxspace(initialpoint::Array{Float64,1}, Ntau::Int64)
  # INPUT: Initial point.
  # OUTPUT: JLD file with a sequence o

  dimtau=2;
  s=SobolSeq(dimtau, zeros(dimtau), ones(dimtau)); skip(s,Ntau);
  sobolspace = Array(Float64,(Ntau,dimtau));
  i=1;
  for j=1:Ntau
    aux =next(s);
    if aux[1]< aux[2];
      sobolspace[i,:]= aux;
      i+=1;
    end
  end
  sobol2=Array(Float64,(i-1,dimtau));

  sobol2 = deepcopy(sobolspace[1:i-1,:])

  return sobol2
end

######################################################################################################################



function close_gov_taue!(govexp::Float64,tau0::Taxes, pa::Param; update::Float64 =0.9, verbose = true, tol::Float64 = 5.0*10.0^-3.0, returnall::Bool = false, wguess::Float64 = 0.55, updateVFIguess::Bool =true, outsideparallel::Bool =true)
  # INPUT: Taxes on dividends and interests + economy constants including government expenditure target, and tau.l, tau.g in tauhat.
  # It MODIFIES tauc in placesuch that the government budget constraint is satisfied
  # OUTPUT steady state welfare

  #0. Set level of parallelization
  if outsideparallel
    VFIfunc = firmVFI!
  else
    VFIfunc = firmVFIParallelOmega!
  end

  #1. Start at the lower bound for taue
  tau0.d = max(1 - (1-tau0.i)/(1-tau0.c) , 0)
  tau0.g = max(1 - (1-tau0.i)/(1-tau0.c) , 0)


  #2. Compute Equilibirum under current taxes
  verbose && println("it=0"," New rates: d = ", tau0.d, " c = ", tau0.c, " i = ", tau0.i, " g = ", tau0.g);
  pr0, eq0 = SolveSteadyState(tau0, pa; wguess = wguess, VFItol =10.0^-3.0,  VFIfunction = VFIfunc, displayit0 = false, displayw = false);

  #3. If goverment contraint is satisfied decrease taxes
  verbose && println("revenue = ", eq0.a.G, " expenditure = ", govexp);
  if eq0.a.G >= govexp  #WE CAN DECREASE TAXES
    #3.0 While budget is not balanced,
    itcount=1; maxit = 500; relgap = (eq0.a.G- govexp)/govexp;
    verbose && println("Relative budget deficit = ", relgap);
    while abs(relgap)>tol && itcount < maxit;
      #3.1 Guess common rate decrease that balances the budget
      corpbase = eq0.a.collections.c / tau0.c;
      intbase = eq0.a.collections.i / tau0.i

      x = (govexp - eq0.a.G)/(corpbase + intbase)
      tau0.c =  update*tau0.c + (1-update)*(tau0.c+x);
      tau0.i =  update*tau0.i + (1-update)*(tau0.i+x);

      verbose && println("it= ",itcount," New rates: d = ", tau0.d, " c = ", tau0.c, " i = ", tau0.i, " g = ", tau0.g);

      #3.2 Compute equilibrium
      if updateVFIguess
        guessVFI = deepcopy(pr0.firmvaluegrid);
      else
        guessVFI = repmat(pa.omega.grid,1,pa.Nz);
      end
      initialradius = min(abs(eq0.w-wguess),10.0^-2.0); #How far to look for wages
      wguess=eq0.w
      pr0, eq0 = SolveSteadyState(tau0, pa; wguess = wguess, VFItol =10.0^-3.0,  VFIfunction = VFIfunc, displayit0 = false, displayw = false, firmvalueguess = guessVFI, initialradius = initialradius);

      #3.3 Update bugdet deficit
      relgap = (eq0.a.G - govexp)/govexp;
      verbose && println("it= ", itcount, " Relative budget deficit = ", relgap);
      #5.5 Update iteration counter
      itcount +=1;
    end
  else #WE NEED TO RAISE EXTRA MONEY
    #4. Keep track of how the revenue is changing
    flag=0; previousrevenue=eq0.a.G;

    #5. While the government constraint is not satisfied,
    itcount=1; maxit = 500; relgap = (eq0.a.G- govexp)/govexp;
    verbose && println("Relative budget deficit = ", relgap);
    while abs(relgap)>tol && itcount < maxit;
      #5.1 Compute tax bases for taxes that will change
      divbase = eq0.a.collections.d / tau0.d;
      #5.2 Update tauc
      tau0.d = tau0.d +  (1-update)*(govexp - eq0.a.G)/divbase;
      tau0.g = tau0.d;
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
      verbose && println("it= ", itcount, " Relative budget deficit = ", relgap, " Revenue = ", eq0.a.G);
      #5.5 Update iteration counter
      itcount +=1;

      #If revenue decreases for three consecutive tax increases, we assume we hit the top of the Laffer curve.
      if eq0.a.G< previousrevenue
        flag+=1;
      else
        flag=0;
      end

      flag == 3 && return -Inf
      previousrevenue = eq0.a.G

    end
  end
  #6.1 The function allows to return the entire equilibrium
  returnall && return eq0.a.welfare, pr0, eq0
  #6.2 Default is to just return welafare
  eq0.a.welfare, eq0.w
end



function maximize_welfare_tesla(taxspace::Array{Float64,2},govexp::Float64, taul::Float64, wguess0::Float64, pa::Param)
  # INPUT:
  # OUTPUT: welfare maximizing tax vector (taud, tauc, taui) given
  Ntau,dimtau=size(taxspace)

  #1. Allocate memory
  welfarevec = Array(Float64,(Ntau,));
  taxvec = Array(Taxes,(Ntau,));
  #4. Loop over taxes moving them slowly to neighboring combinations
  for j=1:Ntau
    println(j,j,j,j,j,j,j,j,j,j,j,j,j,j,j,j,j,j,j,j,j,j,j,j,j,j,j,j,j,j,j,j,j,j,j,j,j,j,j,j,j,j,j,j,j,j,j,j,j,j,j,j)
    tau0 = Taxes(0.0,taxspace[j,1],taxspace[j,2],0.0,taul)
    welfarevec[j],wguess = close_gov_taue!(govexp,tau0, pa; update=0.75, verbose = true, wguess= wguess0, updateVFIguess = true, outsideparallel= false)
    taxvec[j]= deepcopy(tau0)
    #close_gov_vec! updates tau0, eq0, pr0, in place to the new equilibrium under taxes tauhat
    if mod1(j,50)==1
      save("maxwelfare.jld","welfarevec",welfarevec,"taxvec",taxvec)
    end
  end
end
