

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
  tau0.d = max(1 - (1-tau0.i)/(1-tau0.g) , 0)
  tau0.g = max(1 - (1-tau0.i)/(1-tau0.g) , 0)


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
