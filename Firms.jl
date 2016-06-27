#Initialize firm problem

function init_firmproblem(p::Equilibrium, tau::Taxes, fp::FirmParam, hp::HouseholdParam ; guessvalue = false, firmvalueguess::Matrix= ones(fp.Nomega,fp.Nz))

  #Nk, Nq, Nz, Nomega = fp.Nk, fp.Nq, fp.Nz, fp.Nomega

  betatilde = (1.0 + (1-tau.i)/(1-tau.g)*p.r )^(-1);
  taudtilde = 1-(1-tau.d)/(1-tau.g);

  #guess firm value
  if !guessvalue
    firmvalueguess = repmat(omega.grid,1,fp.Nz);
  end
  firmvaluegrid  = copy(firmvalueguess);
  kpolicygrid    = similar(firmvaluegrid);
  qpolicygrid    = similar(firmvaluegrid);

  #Preallocate results
  distributions = Array(Float64,(pr.Nomega,pr.Nz));
  financialcosts = zeros(Float64,(pr.Nomega,pr.Nz));
  grossdividends = zeros(Float64,(pr.Nomega,pr.Nz));
  grossequityis = zeros(Float64,(pr.Nomega,pr.Nz));
  exitprobability = Array(Float64,(pr.Nomega,pr.Nz));
  exitrule = falses(pr.Nomega,pr.Nz,pr.Nz);
  positivedistributions = falses(pr.Nomega,pr.Nz);

  #Initialize the firm problem object
  FirmProblem( betatilde, taudtilde, omega, firmvalueguess, kprime, qprime, betatilde*(1+(1-tau.c)*p.r),
   firmvaluegrid, kpolicygrid, qpolicygrid ,
   map(x->CoordInterpGrid(omega.grid,firmvalueguess[:,x],BCnearest, InterpLinear),1:fp.Nz),
   distributions, financialcosts, grossdividends, grossequityis, exitprobability, exitrule, positivedistributions);
end



#Interpolate current grid for value function
function firmvaluefunction(omegaprime::Real,i_z::Int,pr::FirmProblem)
  pr.InterpolationGrid[i_z][omegaprime]
end

# compute continuation value, given controls and z
function continuation(kprime::Real, qprime::Real, i_z::Int, p::Equilibrium, pr::FirmProblem, tau::Taxes, fp::FirmParam)
  cont =0.0;
  exitvalue = (1-pr.taudtilde)*(fp.kappa*(1-fp.delta)*kprime - (1+p.r)*qprime) ;
  for (i_zprime, zprime) in enumerate(fp.zgrid)
    omegaprime = omegaprimefun(kprime,qprime,i_zprime,p,tau,fp);
    cont += max(exitvalue, firmvaluefunction(omegaprime,i_zprime,pr))*fp.ztrans[i_zprime,i_z];
    end
  return cont
end


#Objective functions for maximization step
# Positive distributions
function objectivefun(kprime::Real, qprime::Real, omega::Real, i_z::Int,p::Equilibrium, pr::FirmProblem, tau::Taxes, fp::FirmParam)
  grossdistributions(omega,kprime,qprime,fp)>=0 ?
    (1-pr.taudtilde)*grossdistributions(omega,kprime,qprime,fp) + pr.betatilde*continuation(kprime, qprime, i_z, p, pr, tau, fp):
    (1+fp.lambda1)*grossdistributions(omega,kprime,qprime,fp) - fp.lambda0 + pr.betatilde*continuation(kprime, qprime, i_z, p, pr, tau, fp) ;
end

#Maximization brute force
function maximizationbf(omega::Real, i_z::Int, p::Equilibrium, pr::FirmProblem, tau::Taxes, fp::FirmParam, extractpolicies::Bool)
  max = -Inf;
  kprimestar::Real = NaN;
  qprimestar::Real = NaN;

  for kprime in pr.kprime.grid
    for qprime in pr.qprime.lb:pr.qprime.step:fp.collateral_factor*kprime
      objective = objectivefun(kprime, qprime, omega, i_z, p, pr, tau, fp);
      if objective > max
        max = objective;
        kprimestar = kprime;
        qprimestar = qprime;
      end
    end
    qprime = fp.collateral_factor*kprime;
    objective = objectivefun(kprime, qprime, omega, i_z,p, pr, tau, fp);
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
function maximizationstep(omega::Real, i_z::Int, p::Equilibrium, pr::FirmProblem, tau::Taxes, fp::FirmParam, extractpolicies::Bool)
  max = -Inf;
  kprimestar::Real = NaN;
  qprimestar::Real = NaN;
  #CASE 1: Firm never distributes dividends
  if pr.discounted_interest>1 && pr.taudtilde >= 0
    error("Firm never distributes dividends")
  #CASE 2: When issuing equity, firm always use as much debt as possible
  elseif pr.discounted_interest <= 1
    #2.1 Capital below max leverage
    for kprime in pr.kprime.lb:pr.kprime.step:(omega*fp.leverageratio)
      for qprime in (kprime-omega):pr.qprime.step:fp.collateral_factor*kprime
        objective = objectivefun(kprime, qprime, omega, i_z,p, pr, tau, fp);
        if objective > max
          max = objective;
          kprimestar = kprime;
          qprimestar = qprime;
        end
      end
      #Evaluate at collateral (may not be on the grid)
      qprime = fp.collateral_factor*kprime;
      objective = objectivefun(kprime, qprime, omega, i_z, p, pr, tau, fp);
      if objective > max
        max = objective;
        kprimestar = kprime;
        qprimestar = qprime;
      end
    end
    #2.2 At zero dividends  (may not be on the grid)
    kprime= omega*fp.leverageratio;
    qprime = fp.collateral_factor*kprime;
    objective = objectivefun(kprime, qprime, omega, i_z, p, pr, tau, fp);
    if objective > max
      max = objective;
      kprimestar = kprime;
      qprimestar = qprime;
    end
    #2.3 Capital above maximum leverage, debt is always at constraint
    for kprime in (omega*fp.leverageratio):pr.kprime.step:pr.kprime.ub
      qprime = fp.collateral_factor*kprime;
      objective = objectivefun(kprime, qprime, omega, i_z, p, pr, tau, fp);
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
    for kprime in pr.kprime.lb:pr.kprime.step:(omega*fp.leverageratio)
      for qprime in pr.qprime.lb:pr.qprime.step:fp.collateral_factor*kprime;
        objective = objectivefun(kprime, qprime, omega, i_z, p, pr, tau, fp);
        if objective > max
          max = objective;
          kprimestar = kprime;
          qprimestar = qprime;
        end
      end
      #Evaluate at contraint (may not be on the grid)
      qprime = fp.collateral_factor*kprime;
      objective = objectivefun(kprime, qprime, omega, i_z, p, pr, tau, fp);
      if objective > max
        max = objective;
        kprimestar = kprime;
        qprimestar = qprime;
      end
    end
    #3.2 At zero dividends
    kprime= omega*fp.leverageratio;
    qprime = fp.collateral_factor*kprime;
    objective = objectivefun(kprime, qprime, omega, i_z, p, pr, tau, fp);
    if objective > max
      max = objective;
      kprimestar = kprime;
      qprimestar = qprime;
    end
    #3.3 Capital above max leverage
    for kprime in (omega*leverageratio):pr.kprime.step:pr.kprime.ub
      #3.3.1 Interior debt and equity
      for qprime in pr.qprime.lb:pr.qprime.step:fp.collateral_factor*kprime
        objective = objectivefun(kprime, qprime, omega, i_z, p, pr, tau, fp)
        if objective > max
          max = objective;
          kprimestar = kprime;
          qprimestar = qprime;
        end
      end
      #3.3.2 Debt at contraint (may not be on the grid)
      qprime = fp.collateral_factor*kprime;
      objective = objectivefun(kprime, qprime, omega, i_z, p, pr, tau, fp);
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


function firmVFIParallel!(pr::FirmProblem, p::Equilibrium, tau::Taxes, fp::FirmParam, tol=10.0^-3, maxit=500; mp=false)
  dist=Inf;
  dif=similar(pr.firmvalueguess);
  it=1;
  while dist > tol
    println("it=", it);
    firmbellmanParallel!(pr,p,tau,fp);
    #firmbellmanParallelOmega!(pr,p,tau,fp);
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

    pr.InterpolationGrid = map(x->CoordInterpGrid(pr.omega.grid,pr.firmvalueguess[:,x],BCnearest, InterpLinear),1:pr.Nz)
    it>= maxit ? error("maximum number of iterations reached"):it+=1;
  end
end




function firmbellmanParallel!(pr::FirmProblem, p::Equilibrium, tau::Taxes, fp::FirmParam)

  #Update value function

    # For every state
    input = [modelState(z,p,pr,tau,fp) for z in 1:pr.Nz];
    resultVector = pmap(maximizeParallel,input)

    for i in 1:length(resultVector)
      # three dimensional linear indexing
      #q,z = ind2sub((model.nQ,model.nZ),i)
      #anArray = results[thisIndex] # two dimensional linear indexing
      # print(resultVector)
      pr.firmvaluegrid[:,i]  = resultVector[i][:,1]
      pr.kpolicygrid[:,i]    = resultVector[i][:,2]
      pr.qpolicygrid[:,i]    = resultVector[i][:,3]
    end
end

type modelState
  i_z::Int64
  p::Equilibrium
  pr::FirmProblem
  tau::Taxes
  fp::FirmParam
end

function maximizeParallel(currentState::modelState)

  # First column of results has the value function
  # Second has K choice
  # Third has Q Choice

  results = Array(Float64,currentState.pr.Nomega,3)

  for (i_omega,omega) in enumerate(currentState.pr.omega.grid)
    results[i_omega,1], results[i_omega,2], results[i_omega,3] = maximizationstep(omega, currentState.i_z,currentState.p ,currentState.pr, currentState.tau, currentState.fp,true);
    # i_omega+=1;
  end

  results

end


##### Parallelize over omega


type modelStateOmega
  i_omega::Int64
  p::Equilibrium
  pr::FirmProblem
  tau::Taxes
  fp::FirmParam
end

function maximizeParallelOmega(currentState::modelStateOmega)

  # First column of results has the value function
  # Second has K choice
  # Third has Q Choice

  results = Array(Float64,currentState.pr.Nz,3)

  for i_z in 1:currentState.pr.Nz
    results[i_z,1], results[i_z,2], results[i_z,3] = maximizationstep(currentState.pr.omega.grid[currentState.i_omega], i_z,currentState.p ,currentState.pr, currentState.tau, currentState.fp,true);
    # i_omega+=1;
  end

  results

end

function firmbellmanParallelOmega!(pr::FirmProblem, p::Equilibrium, tau::Taxes, fp::FirmParam)

  #Update value function

    # For every state
    input = [modelStateOmega(i_omega,p,pr,tau,fp) for i_omega in 1:pr.Nomega];
    resultVector = pmap(maximizeParallelOmega,input);

    for i in 1:length(resultVector)
      # three dimensional linear indexing
      #q,z = ind2sub((model.nQ,model.nZ),i)
      #anArray = results[thisIndex] # two dimensional linear indexing
      # print(resultVector)
      pr.firmvaluegrid[i,:]  = resultVector[i][:,1]
      pr.kpolicygrid[i,:]    = resultVector[i][:,2]
      pr.qpolicygrid[i,:]    = resultVector[i][:,3]
    end
end


function firmVFIParallelOmega!(pr::FirmProblem, p::Equilibrium, tau::Taxes, fp::FirmParam, tol=10.0^-3, maxit=500; mp=false)
  dist=Inf;
  dif=similar(pr.firmvalueguess);
  it=1;
  while dist > tol
    println("it=", it);
    firmbellmanParallelOmega!(pr,p,tau,fp);
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

    pr.InterpolationGrid = map(x->CoordInterpGrid(pr.omega.grid,pr.firmvalueguess[:,x],BCnearest, InterpLinear),1:pr.Nz)
    it>= maxit ? error("maximum number of iterations reached"):it+=1;
  end
end
