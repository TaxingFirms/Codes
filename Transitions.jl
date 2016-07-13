# This script takes an initial equilibrium, a final equilibrium and a path for
# taxes and return a path for equilibrium prices and allocations

type PeriodSolution
  fpr::FirmProblem
  eq::Equilibrium
end

function init_transitions(T::Int64)
  transitions= Array(PeriodSolution,(T,));
end

function transitions!(tr:: Array(PeriodSolution,(T,)), T::Int64, ss0::PeriodSolution, ssT::PeriodSolution, tauseq::Array{Taxes,1}, pa::Param; update::Float64=0.9, tol = 10^-4.0, maxit=10^4.0 , maxroutine::Function=maximizationstep )


    # Guess paths for r, w, E
    for t=1:T-1
      wguess= t*(ssT.eq.w - ss0.eq.w)/T + ss0.eq.w;
      rguess= t*(ssT.eq.r - ss0.eq.r)/T + ss0.eq.r;
      eq = init_equilibirium(wguess,tau,pa; r=rguess);
      fpr = init_firmproblem(eq,tau,pa);
      tr[t] = PeriodSolution(fpr,eq);
      tr[t].eq.E= t*(ssT.eq.E - ss0.eq.E)/T + ss0.eq.E;
    end

    # Assume new steady state is reached after T periods
    tr[T] = ssT;






    #Initialize array
    expvalentry= Array(Float64,(T,));

    distance=+Inf;
    it=1;
    while distance>tol  && it<maxit
      # Solve for path for prices using a shooting algorithm
      updateprices!(tr:: Array(PeriodSolution,(T,)), T::Int64);

      # Free entry condition
      for t=1:T
        expvalentry[t]=compute_expvalentry(tr[t].fpr, pa);
      end

      distance = norm(expvalentry,Inf);
      it+=1;
      if distance >tol
        for t=1:T
          tr[t].eq.E = (expvalentry[t] - tol)/tol
        end
    end


end



##############################################


##############################################



function firm_aggregates_transitions(crt_fpr::FirmProblem, prev_fpr::FirmProblem, crt_eq::Equilibrium, prev_eq::Equilibrium, crt_tau::Taxes, pa::Param)
# Computes labor demand and net supply of goods at time t. Where time t is crt.

  prev_distr=prev_eq.distr;

  crt_distr= crt_eq.distr
  wage=crt_eq.w;
  irate= crt_eq.r;

  # Compute aggregates depending on current K
  lprime_d=0.0;
  gdp=0.0;
  corptax=0.0;
  for i_zprime in 1:pa.Nz
    zprime= pa.zgrid[i_zprime];
    for i_z in 1:pa.Nz
      probzprime=pa.ztrans[i_zprime,i_z];

      for i_omega in 1:pa.Nomega
        kprime=prev_fpr.kpolicy[i_omega,i_z];
        qprime=prev_fpr.qpolicy[i_omega,i_z];
        omega = pa.omega.grid[i_omega];
        mass = prev_distr[i_omega,i_z];

        if !prev_fpr.exitrule[i_omega,i_z,i_zprime]
          lprime    = (zprime*pa.alphal*kprime^pa.alphak / wage)^(1/(1-pa.alphal));
          lprime_d += lprime*mass*probzprime;
          gdp      += (zprime*kprime^pa.alphak*lprime^pa.alphal - pa.f)*mass*probzprime;
          corptax  += tau.c*(zprime*kprime^pa.alphak*lprime^pa.alphal -wage*lprime - pa.delta*kprime - irate*qprime - pa.f)*mass*probzprime;
        end
      end
    end
  end

  # Compute aggregates that do not depend on K
  capital= sum(prev_distr.*prev_fpr.kpolicy);
  investment = sum(crt_distr.*crt_fpr.kpolicy) - (1-pa.delta)*capital;
  grossdividends = sum(crt_distr.*crt_fpr.grossdividends);
  financialcosts= - sum(crt_distr.*crt_fpr.financialcosts);
  debt = sum(crt_distr.*crt_fpr.qpolicy);
  divtax= tau.d*grossdividends;
  inctax = tau.i*eq.r*debt;
  G = divtax + corptax + inctax;

  goodsnet_s = gdp - G - financialcosts - investment;

  return lprime_d,goodsnet_s;
end


##########################################


##########################################

function updateprices!(tr:: Array(PeriodSolution,(T,)), T::Int64)

  what=Array(Float64,(T,));
  rhat=Array(Float64,(T,));

  distance=+Inf;
  it=1;
  while distance>tol  && it<maxit

    #Solve value function backwards
    norm=0.0;

    for t=T:2
      # Firms
      pr= deepcopy(tr[t].fpr);
      firmbellmanParallelOmega!(pr,tr[t].eq,tausec[t],pa,maxroutine);
      getpolicies!(pr,tr[t].eq,tr[t].tau,pa);
      tr[t-1].fpr = pr;
      #Households
      ## householdProblem!(tr[t-1].hpr, aprime::tr[t].hpr.a,aprime::tr[t].hpr.mu, tr[t].eq.r, tr[t].eq.w, tausec[t], pa);
    end

    #Compute cross-sectional distribution forward (this can be made faster if needed)
    # and aggregate firm problem

    tr[1].eq.distr = ss0.eq.distr;
    tr[1].eq.a.labor, tr[1].eq.a.consumption = firm_aggregates_transitions(tr[1].fpr, SS0.fpr, tr[1].eq, SS0.eq, tr[1].tau, pa);

    for t=2:T
      tr[t].eq.distr = transitionrule(tr[t-1].eq.distr,tr[t].eq.E, tr[t-1].fpr, tr[t].eq , tausec[t], pa);
      tr[t].eq.a.labor, tr[t].eq.a.consumption = firm_aggregates_transitions(tr[t].fpr, tr[t-1].fpr, tr[t].eq, tr[t-1].eq, tr[t].tau, pa);

      # Compute candidate wage
      what[t] = pa.H*(tr[t].eq.a.labor)^1/pa.psi;
      # Compute contribution to price discrepancy
      norm += (what[t] - tr[t].eq.w)^2.0;
      # Update wage (note that updating at this point is fine)
      tr[t].eq.w = update*tr[t].eq.w + (1-update)*what[t];
      if t>2
        Cprime = tr[t].eq.a.consumption;
        C = tr[t-1].eq.a.consumption;
        Lprime = tr[t].eq.a.labor;
        L = tr[t-1].eq.a.labor;
        rhat[t-1] = (1-tr[t].tau.i)^(-1)* (   pa.beta^(-1) *  ( (Cprime - pa.H/(1+1/pa.sigma)*Lprime)/(C - pa.H/(1+1/pa.sigma)*L) )^pa.sigma  -1   );
        # Compute contribution to price discrepancy
        norm += (rhat[t-1] - tr[t-1].eq.r)^2.0;
        # Update interest rate
        tr[t-1].eq.r = update*tr[t-1].eq.r + (1-update)*rhat[t-1];
      end
      if t=T
        C = tr[T].eq.a.consumption;
        L = tr[T].eq.a.labor;
        rhat[T] = (1-tr[T].tau.i)^(-1)* (   pa.beta^(-1) *  ( (C- pa.H/(1+1/pa.sigma)*L)/(C - pa.H/(1+1/pa.sigma)*L) )^pa.sigma  -1   );
        norm += (rhat[T] - tr[T].eq.r)^2.0;
        tr[T].eq.r = update*tr[T].eq.r + (1-update)*rhat[T];
      end
    end
    what[1] = pa.H*(tr[1].eq.a.labor)^1/pa.psi;
    norm += (what[1] - tr[1].eq.w)^2.0;
    tr[1].eq.w = update*tr[1].eq.w + (1-update)*what[1];

    Cprime = tr[2].eq.a.consumption;
    C = tr[1].eq.a.consumption;
    Lprime = tr[2].eq.a.labor;
    L = tr[1].eq.a.labor;
    rhat[1] = (1-tr[2].tau.i)^(-1)* (   pa.beta^(-1) *  ( (Cprime - pa.H/(1+1/pa.sigma)*Lprime)/(C - pa.H/(1+1/pa.sigma)*L) )^pa.sigma  -1   );
    norm += (rhat[1] - tr[1].eq.r)^2.0;
    tr[1].eq.r = update*tr[1].eq.r + (1-update)*rhat[1];


    distance = norm;
    it+=1;
  end
end
