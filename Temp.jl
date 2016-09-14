
#govexp = ref[end].eq.a.G
#update =0.9; verbose = true; tol =10.0^-3.0; returnall= false; wguess = ref[end].eq.w; updateVFIguess=true; outsideparallel=false

function close_gov_tauc!(govexp::Float64,tau0::Taxes, pa::Param; update::Float64 =0.9, verbose = true, tol::Float64 = 5.0*10.0^-3.0, returnall::Bool = false, wguess::Float64 = 0.55, updateVFIguess::Bool =true, outsideparallel::Bool =true)
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

  #3. If goverment contraint is satisfied decrease taxes
  verbose && println("revenue = ", eq0.a.G, " expenditure = ", govexp);
  if eq0.a.G >= govexp  #WE CAN DECREASE TAXES
    #3.0 While budget is not balanced,
    itcount=1; maxit = 500; relgap = (eq0.a.G- govexp)/govexp;
    verbose && println("Relative budget deficit = ", relgap);
    while abs(relgap)>tol && itcount < maxit;
      #3.1 Guess common rate decrease that balances the budget
      gainsbase = eq0.a.collections.g / tau0.g;
      divbase = eq0.a.collections.d / tau0.d;
      intbase = eq0.a.collections.i / tau0.i

      x = (govexp - eq0.a.G)/(gainsbase + divbase + intbase)
      tau0.d =  update*tau0.d + (1-update)*(tau0.d+x);
      tau0.i =  update*tau0.i + (1-update)*(tau0.i+x);
      tau0.g =  update*tau0.g + (1-update)*(tau0.g+x);

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
      corpbase = eq0.a.collections.c / tau0.c;
      #5.2 Update tauc
      tau0.c =  tau0.c + (1-update)*(govexp - eq0.a.G)/corpbase;
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
  eq0.a.welfare
end


######################################################################################################################


function close_gov_taud(govexp::Float64,tauhat::Taxes, pa::Param; update::Float64 =0.9, verbose = true, tol::Float64 = 5.0*10.0^-3.0, returnall::Bool = false, wguess::Float64 = 0.55, updateVFIguess::Bool =true, outsideparallel::Bool =true)
  # INPUT: Taxes on dividends and interests + economy constants including government expenditure target, and tau.l, tau.g in tauhat.
  # It MODIFIES tauc in placesuch that the government budget constraint is satisfied
  # OUTPUT steady state welfare

  #0. Set level of parallelization
  if outsideparallel
    VFIfunc = firmVFI!
  else
    VFIfunc = firmVFIParallelOmega!
  end

  #1. Start at the lower bound for taud
  tau0 = deepcopy(tauhat) #so tauhat is not modified
  tau0.i >  tau0.g && return - Inf
  tau0.d = 0.0;
  tau0.c = 0.0;

  #2. Compute Equilibirum under current taxes
  verbose && println("it=0"," New rates: d = ", tau0.d, " c = ", tau0.c, " i = ", tau0.i, " g = ", tau0.g);
  pr0, eq0 = SolveSteadyState(tau0, pa; wguess = wguess, VFItol =10.0^-3.0,  VFIfunction = VFIfunc, displayit0 = true, displayw = true);

  #3. If goverment contraint is satisfied decrease taxes taui and taug
  verbose && println("revenue = ", eq0.a.G, " expenditure = ", govexp);
  if eq0.a.G >= govexp  #WE CAN DECREASE TAXES
    #3.0 While budget is not balanced,
    itcount=1; maxit = 500; relgap = (eq0.a.G- govexp)/govexp;
    verbose && println("Relative budget deficit = ", relgap);
    while abs(relgap)>tol && itcount < maxit;
      #3.1 Guess common rate decrease that balances the budget
      gainsbase = eq0.a.collections.g / tau0.g;
      intbase = eq0.a.collections.i / tau0.i

      x = (govexp - eq0.a.G)/(gainsbase + intbase)
      tau0.i =  update*tau0.i + (1-update)*(tau0.i+x);
      tau0.g =  update*tau0.g + (1-update)*(tau0.g+x);

      verbose && println("it= ",itcount," New rates: d = ", tau0.d, " c = ", tau0.c, " i = ", tau0.i, " g = ", tau0.g);

      #3.2 Compute equilibrium
      if updateVFIguess
        guessVFI = deepcopy(pr0.firmvaluegrid);
      else
        guessVFI = repmat(pa.omega.grid,1,pa.Nz);
      end
      initialradius = min(abs(eq0.w-wguess),10.0^-2.0); #How far to look for wages
      wguess=eq0.w
      pr0, eq0 = SolveSteadyState(tau0, pa; wguess = wguess, VFItol =10.0^-3.0,  VFIfunction = VFIfunc, displayit0 = true, displayw = true, firmvalueguess = guessVFI, initialradius = initialradius);

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
      #5.2 Update taud
      tau0.d =  tau0.d + (1-update)*(govexp - eq0.a.G)/corpbase;
      verbose && println("it= ",itcount," New rates: d = ", tau0.d, " c = ", tau0.c, " i = ", tau0.i, " g = ", tau0.g);

      #5.3 Solve equilibrium for candidate taxes
      if updateVFIguess
        guessVFI = deepcopy(pr0.firmvaluegrid);
      else
        guessVFI = repmat(pa.omega.grid,1,pa.Nz);
      end
      initialradius = min(abs(eq0.w-wguess),10.0^-2.0); #How far to look for wages
      wguess=eq0.w
      pr0, eq0 = SolveSteadyState(tau0, pa; wguess = wguess, VFItol =10.0^-3.0,  VFIfunction = VFIfunc, displayit0 = true, displayw = true, firmvalueguess = guessVFI, initialradius = initialradius);

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

  #6 Return welfare and taxes
  [eq0.a.welfare , tau0.d, tau0.i, tau0.g]
end

immutable Argument
  govexp::Float64
  tau0::Taxes
  pa::Param
end


function taxequilibrium(arg::Argument)
  close_gov_tauc!(arg.govexp,arg.tau0, arg.pa; update=0.75, verbose = true, wguess= 0.53, updateVFIguess = true, outsideparallel= true)
end


function maximize_welfare_hawk(taxspace::Array{Float64,2},govexp::Float64, taul::Float64, wguess0::Float64, pa::Param)
  # INPUT:
  # OUTPUT: welfare maximizing tax vector (taud, tauc, taui) given
  initialpoint = [0.0, 0.28];
  taxspace=create_taxspace(initialpoint, 10; jldfile = false)

  Nexp,dimtau=size(taxspace)

  #Allocate memory for vector of arguments
  args = Array(Argument,(Nexp,))
  for j=1:Nexp
    tau0 = Taxes(0.0,taxspace[j,1],taxspace[j,2],0.0,taul)
    args[j] = Argument(govexp,tau0,pa)
  end


  #Call map function
  welf_and_taxes = map(taxequilibrium,args);


end
