# This script has the functions needed to run the value function iteration.
# The value function iteration takes prices from eq::Equilibrium,
# taxes from tau::Taxes and parameters from pa::Param and it computes the
# value and policy functions in pr::FirmProblem
#
# We have 2 functions that can be called to perform the VFI: firmVFIParallel!,
# firmVFIParallel!
#
# We have 2 functions that compute the maximization step in the Bellman
# operator: maximizationbf, maximizationstep


#Interpolate current grid for value function
function firmvaluefunction(omegaprime::Real,i_z::Int,pr::FirmProblem)
  pr.InterpolationGrid[i_z][omegaprime]
end

# compute continuation value, given controls and z
function continuation(kprime::Real, qprime::Real, i_z::Int, pr::FirmProblem, eq::Equilibrium, tau::Taxes, pa::Param)
  cont =0.0;
  exitvalue = (1-pr.taudtilde)*(pa.kappa*(1-pa.delta)*kprime - (1+eq.r)*qprime) ;
  for (i_zprime, zprime) in enumerate(pa.zgrid)
    omegaprime = omegaprimefun(kprime,qprime,i_zprime,eq,tau,pa);
    cont += max(exitvalue, firmvaluefunction(omegaprime,i_zprime,pr))*pa.ztrans[i_zprime,i_z];
    end
  return cont
end


#Objective functions for maximization step
# Positive distributions
function objectivefun(kprime::Real, qprime::Real, omega::Real, i_z::Int, pr::FirmProblem, eq::Equilibrium, tau::Taxes, pa::Param)
  grossdistributions(omega,kprime,qprime,pa)>=0 ?
    (1-pr.taudtilde)*grossdistributions(omega,kprime,qprime,pa) + pr.betatilde*continuation(kprime, qprime, i_z, pr, eq, tau, pa):
    (1+pa.lambda1)*grossdistributions(omega,kprime,qprime,pa) - pa.lambda0 + pr.betatilde*continuation(kprime, qprime, i_z, pr, eq, tau, pa);
end

#Maximization brute force
function maximizationbf(omega::Real, i_z::Int, eq::Equilibrium, pr::FirmProblem, tau::Taxes, pa::Param, extractpolicies::Bool)
  max = -Inf;
  kprimestar::Real = NaN;
  qprimestar::Real = NaN;

  for kprime in pa.kprime.grid
    for qprime in pa.qprime.lb:pa.qprime.step:pa.collateral_factor*kprime
      objective = objectivefun(kprime, qprime, omega, i_z, pr, eq, tau, pa);
      if objective > max
        max = objective;
        kprimestar = kprime;
        qprimestar = qprime;
      end
    end
    qprime = pa.collateral_factor*kprime;
    objective = objectivefun(kprime, qprime, omega, i_z,pr, eq, tau, pa);
    if objective > max
      max = objective;
      kprimestar = kprime;
      qprimestar = qprime;
    end
  end
    extractpolicies?
    (max, kprimestar, qprimestar):
    max
end


# Maximization step
function maximizationstep(omega::Real, i_z::Int, eq::Equilibrium, pr::FirmProblem, tau::Taxes, pa::Param, extractpolicies::Bool)
  max = -Inf;
  kprimestar::Real = NaN;
  qprimestar::Real = NaN;
  #CASE 1: Firm never distributes dividends
  if pr.discounted_interest>1 && pr.taudtilde >= 0
    error("Firm never distributes dividends")
  #CASE 2: When issuing equity, firm always use as much debt as possible
  elseif pr.discounted_interest <= 1
    #2.1 Capital below max leverage
    for kprime in pa.kprime.lb:pa.kprime.step:(omega*pa.leverageratio)
      for qprime in (kprime-omega):pa.qprime.step:pa.collateral_factor*kprime
        objective = objectivefun(kprime, qprime, omega, i_z,pr, eq, tau, pa);
        if objective > max
          max = objective;
          kprimestar = kprime;
          qprimestar = qprime;
        end
      end
      #Evaluate at collateral (may not be on the grid)
      qprime = pa.collateral_factor*kprime;
      objective = objectivefun(kprime, qprime, omega, i_z, pr, eq, tau, pa);
      if objective > max
        max = objective;
        kprimestar = kprime;
        qprimestar = qprime;
      end
    end
    #2.2 At zero dividends  (may not be on the grid)
    kprime= omega*pa.leverageratio;
    qprime = pa.collateral_factor*kprime;
    objective = objectivefun(kprime, qprime, omega, i_z, pr, eq, tau, pa);
    if objective > max
      max = objective;
      kprimestar = kprime;
      qprimestar = qprime;
    end
    #2.3 Capital above maximum leverage, debt is always at constraint
    for kprime in (omega*pa.leverageratio):pa.kprime.step:pa.kprime.ub
      qprime = pa.collateral_factor*kprime;
      objective = objectivefun(kprime, qprime, omega, i_z, pr, eq, tau, pa);
      if objective > max
        max = objective;
        kprimestar = kprime;
        qprimestar = qprime;
      end
    end

  #CASE 3: Optimal debt may be interior for both distribution and equity issuances (this will be important if lambda2>0)
  else
    println("Shouldn't be in this loop")
    #3.1 Capital below max leverage
    for kprime in pa.kprime.lb:pa.kprime.step:(omega*pa.leverageratio)
      for qprime in pa.qprime.lb:pa.qprime.step:pa.collateral_factor*kprime;
        objective = objectivefun(kprime, qprime, omega, i_z, pr, eq, tau, pa);
        if objective > max
          max = objective;
          kprimestar = kprime;
          qprimestar = qprime;
        end
      end
      #Evaluate at contraint (may not be on the grid)
      qprime = pa.collateral_factor*kprime;
      objective = objectivefun(kprime, qprime, omega, i_z, pr, eq, tau, pa);
      if objective > max
        max = objective;
        kprimestar = kprime;
        qprimestar = qprime;
      end
    end
    #3.2 At zero dividends
    kprime= omega*pa.leverageratio;
    qprime = pa.collateral_factor*kprime;
    objective = objectivefun(kprime, qprime, omega, i_z, pr, eq, tau, pa);
    if objective > max
      max = objective;
      kprimestar = kprime;
      qprimestar = qprime;
    end
    #3.3 Capital above max leverage
    for kprime in (omega*leverageratio):pa.kprime.step:pa.kprime.ub
      #3.3.1 Interior debt and equity
      for qprime in pa.qprime.lb:pa.qprime.step:pa.collateral_factor*kprime
        objective = objectivefun(kprime, qprime, omega, i_z, pr, eq, tau, pa)
        if objective > max
          max = objective;
          kprimestar = kprime;
          qprimestar = qprime;
        end
      end
      #3.3.2 Debt at contraint (may not be on the grid)
      qprime = pa.collateral_factor*kprime;
      objective = objectivefun(kprime, qprime, omega, i_z, pr, eq, tau, pa);
      if objective > max
        max = objective;
        kprimestar = kprime;
        qprimestar = qprime;
      end
    end
  end
  extractpolicies?
    (max, kprimestar, qprimestar):
    max
end


function firmVFIParallel!(pr::FirmProblem, eq::Equilibrium, tau::Taxes, pa::Param; tol=10.0^-3, maxit=500, mp=false, maximizationroutine=maximizationstep)
  dist=Inf;
  dif=similar(pr.firmvalueguess);
  it=1;
  while dist > tol
    println("it=", it);
    firmbellmanParallel!(pr,eq,tau,pa, maximizationroutine);
    dif = pr.firmvalueguess - pr.firmvaluegrid;
    dist= norm(dif, Inf);
    println("it=", it, "   dist=", dist);

    # McQueen - Porteus Accelerating thing
    if mp
    bUnder = minimum(pr.firmvaluegrid - pr.firmvalueguess)
    bOver  = maximum(pr.firmvaluegrid - pr.firmvalueguess)
    pr.firmvalueguess = deepcopy(pr.firmvaluegrid) .+ (pr.betatilde/(1-pr.betatilde))*(bUnder + bOver)/2
    else
    pr.firmvalueguess = deepcopy(pr.firmvaluegrid);
    end

    pr.InterpolationGrid = map(x->CoordInterpGrid(pa.omega.grid,pr.firmvalueguess[:,x],BCnearest, InterpLinear),1:pa.Nz)
    it>= maxit ? error("maximum number of iterations reached"):it+=1;
  end
end




function firmbellmanParallel!(pr::FirmProblem, eq::Equilibrium, tau::Taxes, pa::Param, maximizationroutine::Function)

  #Update value function

    # For every state
    input = [modelState(z,pr,eq,tau,pa, maximizationroutine) for z in 1:pa.Nz];
    resultVector = pmap(maximizeParallel,input)

    for i in 1:length(resultVector)
      # three dimensional linear indexing
      #q,z = ind2sub((model.nQ,model.nZ),i)
      #anArray = results[thisIndex] # two dimensional linear indexing
      # print(resultVector)
      pr.firmvaluegrid[:,i]  = resultVector[i][:,1]
      pr.kpolicy[:,i]    = resultVector[i][:,2]
      pr.qpolicy[:,i]    = resultVector[i][:,3]
    end
end

type modelState
  i_z::Int64
  pr::FirmProblem
  eq::Equilibrium
  tau::Taxes
  pa::Param
  maxroutine::Function
end

function maximizeParallel(currentState::modelState)

  # First column of results has the value function
  # Second has K choice
  # Third has Q Choice
  results = Array(Float64,currentState.pa.Nomega,3)
  routine=currentState.maxroutine;

  for (i_omega,omega) in enumerate(currentState.pa.omega.grid)
    results[i_omega,1], results[i_omega,2], results[i_omega,3] = routine(omega, currentState.i_z,currentState.eq ,currentState.pr, currentState.tau, currentState.pa,true);
    # i_omega+=1;
  end

  results

end


##### Parallelize over omega


type modelStateOmega
  i_omega::Int64
  pr::FirmProblem
  eq::Equilibrium
  tau::Taxes
  pa::Param
  maxroutine::Function
end

function maximizeParallelOmega(currentState::modelStateOmega)

  # First column of results has the value function
  # Second has K choice
  # Third has Q Choice

  results = Array(Float64,currentState.pa.Nz,3)
  routine=currentState.maxroutine;

  for i_z in 1:currentState.pa.Nz
    results[i_z,1], results[i_z,2], results[i_z,3] = routine(currentState.pa.omega.grid[currentState.i_omega], i_z,currentState.eq ,currentState.pr, currentState.tau, currentState.pa,true);
    # i_omega+=1;
  end

  results

end

function firmbellmanParallelOmega!(pr::FirmProblem, eq::Equilibrium, tau::Taxes, pa::Param, maximizationroutine::Function)

  #Update value function

    # For every state
    input = [modelStateOmega(i_omega,pr,eq,tau,pa,maximizationroutine) for i_omega in 1:pa.Nomega];
    resultVector = pmap(maximizeParallelOmega,input);

    for i in 1:length(resultVector)
      # three dimensional linear indexing
      #q,z = ind2sub((model.nQ,model.nZ),i)
      #anArray = results[thisIndex] # two dimensional linear indexing
      # print(resultVector)
      pr.firmvaluegrid[i,:]  = resultVector[i][:,1]
      pr.kpolicy[i,:]    = resultVector[i][:,2]
      pr.qpolicy[i,:]    = resultVector[i][:,3]
    end
end


function firmVFIParallelOmega!(pr::FirmProblem, eq::Equilibrium, tau::Taxes, pa::Param; maximizationroutine::Function=maximizationstep, tol=10.0^-3, maxit=500, mp=false)
  dist=Inf;
  dif=similar(pr.firmvalueguess);
  it=1;
  while dist > tol
    println("it=", it);
    firmbellmanParallelOmega!(pr,eq,tau,pa,maximizationroutine);
    dif = pr.firmvalueguess - pr.firmvaluegrid;
    dist= norm(dif, Inf);
    println("it=", it, "   dist=", dist);

    # McQueen - Porteus Accelerating thing
    if mp
    bUnder = minimum(pr.firmvaluegrid - pr.firmvalueguess)
    bOver  = maximum(pr.firmvaluegrid - pr.firmvalueguess)
    pr.firmvalueguess = deepcopy(pr.firmvaluegrid) .+ (pr.betatilde/(1-pr.betatilde))*(bUnder + bOver)/2
    else
    pr.firmvalueguess = deepcopy(pr.firmvaluegrid);
    end

    pr.InterpolationGrid = map(x->CoordInterpGrid(pa.omega.grid,pr.firmvalueguess[:,x],BCnearest, InterpLinear),1:pa.Nz)
    it>= maxit ? error("maximum number of iterations reached"):it+=1;
  end
end
