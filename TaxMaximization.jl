
function close_gov_tauc(govexp::Float64,tauhat::Taxes, pa::Param; update::Float64 =0.75, verbose = true, tol::Float64 = 5.0*10.0^-3.0,taxtol::Float64 =10.0^-6.0, returnall::Bool = false, wguess::Float64 = 0.55, updateVFIguess::Bool =true, outsideparallel::Bool =true)
  # INPUT: Taxes on dividends and interests + economy constants including government expenditure target, and tau.l, tau.g in tauhat.
  # OUTPUT steady state welfare and tax vector such that BC is satisfied

  #0. Set level of parallelization
  if outsideparallel
    VFIfunc = firmVFI!
  else
    VFIfunc = firmVFIParallelOmega!
  end

  #1. Start at the lower bound for taue
  tau0 =deepcoy(tauhat) #We don't modify taxes in place
  tau0.c = max(1 - (1-tau0.i)/(1-tau0.g) , 0)

  #2. Compute Equilibirum under current taxes
  verbose && println("it=0"," New rates: d = ", tau0.d, " c = ", tau0.c, " i = ", tau0.i, " g = ", tau0.g);
  pr0, eq0 = SolveSteadyState(tau0, pa; wguess = wguess, VFItol =10.0^-3.0,  VFIfunction = VFIfunc, displayit0 = false, displayw = false);

  #3. If goverment contraint is satisfied decrease taxes
  verbose && println("revenue = ", eq0.a.G, " expenditure = ", govexp);
  if eq0.a.G >= govexp  #WE CAN DECREASE TAXES
    #3.0 While budget is not balanced,
    itcount=1; maxit = 500; relgap = (eq0.a.G- govexp)/govexp;
    oldgap = relgap; oldtau = tau0.d; taxdif= Inf;
    verbose && println("Relative budget surplus = ", relgap);
    while abs(relgap)>tol && itcount < maxit;
      #3.1 Guess common rate decrease that balances the budget
      divbase = eq0.a.collections.d / tau0.d;
      intbase = eq0.a.collections.i / tau0.i;

      if oldgap*relgap<0 #if deficit changes sign, decrease updating speed
        update+=(1-update)/2;
      end

      x = (govexp - eq0.a.G)/(corpbase + divbase)
      tau0.d = max(0, tau0.d + (1-update)*x);
      tau0.g = max(0, tau0.g + (1-update)*x);
      tau0.i = max(0, tau0.i + (1-update)*x);

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
      oldgap=relgap;
      relgap = (eq0.a.G - govexp)/govexp;
      taxdif=abs(oldtau - tau0.d); oldtau=tau0.d;
      verbose && println("it= ", itcount, " Relative budget deficit = ", relgap);
      #5.5 Update iteration counter
      itcount +=1;
    end
  else #WE NEED TO RAISE EXTRA MONEY
    #4. While the government constraint is not satisfied,
    itcount=1; maxit = 500; relgap = (eq0.a.G- govexp)/govexp;
    verbose && println("Relative budget deficit = ", relgap);
    oldgap = relgap; oldtau = tau0.d; taxdif= Inf;
    while abs(relgap)>tol && itcount < maxit;
      #4.1 Compute tax bases for taxes that will change
      corpbase = eq0.a.collections.c / tau0.c;
      #4.2 Update tauc
      if oldgap*relgap<0 #if deficit changes sign, decrease updating speed
        update+=(1-update)/2;
      end

      tau0.c = max(0, tau0.c + (1-update)*(govexp - eq0.a.G)/corpbase);
      verbose && println("it= ",itcount," New rates: d = ", tau0.d, " c = ", tau0.c, " i = ", tau0.i, " g = ", tau0.g);

      #4.2 Check if below peek of Laffer
      if taxdif < taxtol
        taxvec = [tau0.d tau0.c tau0.i tau0.g tau0.l]
        if abs(relgap)>10*tol
          return -Inf, taxvec
        else
          return eq0.a.welfare, taxvec
        end
      end

      #4.3 Solve equilibrium for candidate taxes
      if updateVFIguess
        guessVFI = deepcopy(pr0.firmvaluegrid);
      else
        guessVFI = repmat(pa.omega.grid,1,pa.Nz);
      end
      initialradius = min(abs(eq0.w-wguess),10.0^-2.0); #How far to look for wages
      wguess=eq0.w
      pr0, eq0 = SolveSteadyState(tau0, pa; wguess = wguess, VFItol =10.0^-3.0,  VFIfunction = VFIfunc, displayit0 = false, displayw = false, firmvalueguess = guessVFI, initialradius = initialradius);

      #4.4 Update bugdet deficit
      oldgap=relgap;
      relgap = (eq0.a.G - govexp)/govexp;
      verbose && println("it= ", itcount, " Relative budget deficit = ", relgap, " Revenue = ", eq0.a.G);
      #4.5 Update iteration counter
      itcount +=1;

    end
  end
  #5.1 The function allows to return the entire equilibrium
  returnall && return eq0.a.welfare, pr0, eq0
  #5.2 Default is to just return welafare and taxes
  taxvec = [tau0.d tau0.c tau0.i tau0.g tau0.l]

  return eq0.a.welfare, taxvec
end
