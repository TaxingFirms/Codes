type basicState
  omega::Float64
  i_z::Int64
  pr::FirmProblem
  eq::Equilibrium
  tau::Taxes
  pa::Param
end


function maximize_aux(state::basicState)
  maximizationfast(state.omega, state.i_z, state.eq, state.pr, state.tau, state.pa,true);
end


function firmbellman!(pr::FirmProblem, eq::Equilibrium, tau::Taxes, pa::Param, maxroutine::Function)

  #Update value function

    # For every state
    input = Array(basicState,(pa.Nomega,pa.Nz))
    for (i_omega, omega) in enumerate(pa.omega.grid)
      for i_z=1:pa.Nz
        input[i_omega, i_z] = basicState(omega,i_z,pr,eq,tau,pa)
      end
    end

    resultVector = map(maximize_aux,input);

    for i in 1:length(resultVector)
      pr.firmvaluegrid[i], pr.kpolicy[i], pr.qpolicy[i] = resultVector[i];
    end
end


function firmVFI!(pr::FirmProblem, eq::Equilibrium, tau::Taxes, pa::Param; maxroutine::Function=maximizationstep, tol=10.0^-5, maxit=5000, mp=false, verbose=true )
  dist=Inf;
  dif=similar(pr.firmvalueguess);
  it=1;
  while dist > tol
    verbose && println("it=", it);
    firmbellman!(pr,eq,tau,pa,maxroutine);
    dif = pr.firmvalueguess - pr.firmvaluegrid;
    dist= norm(dif, Inf);
    verbose && println("it=", it, "   dist=", dist);

    if mp
    betatilde = (1.0 + (1-tau.i)/(1-tau.g)*eq.r )^(-1);
    bUnder = minimum(pr.firmvaluegrid - pr.firmvalueguess)
    bOver  = maximum(pr.firmvaluegrid - pr.firmvalueguess)
    pr.firmvalueguess = deepcopy(pr.firmvaluegrid) .+ (betatilde/(1-betatilde))*(bUnder + bOver)/2
    else
    pr.firmvalueguess = deepcopy(pr.firmvaluegrid);
    end

    pr.InterpolationGrid = map(x->CoordInterpGrid(pa.omega.grid,pr.firmvalueguess[:,x],BCnearest, InterpLinear),1:pa.Nz)
    it>= maxit ? error("maximum number of iterations reached"):it+=1;
  end
end
