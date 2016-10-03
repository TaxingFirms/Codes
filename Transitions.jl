# This script takes an initial equilibrium, a final equilibrium and a path for
# taxes and return a path for equilibrium prices and allocations


function updateprices!(T::Int64, tr:: Array{PeriodSolution,1},  ss0::PeriodSolution; verbose=true, update=0.9 , maxroutine=maximizationstep, tol = 10^-5.0, maxit = 10.0^5.0)
  #INPUT: T, length of transition period.
  #       tr, vector of solutions with a guess of w,r for each t. It also has the final steady state as the last element.
  #       ss0: inital steady state
  #OUTPUT: tr is updated in place with the transition dynamics

  #1. Initialize arrays for wage, interest rate and consumption
  what=Array(Float64,(T,));
  rhat=Array(Float64,(T,));
  rdisplay=Array(Float64,(T,));
  cons=Array(Float64,(T,));
  consprime = Array(Float64,(T,));

  #2. Iterate value function backwards,
  distance=+Inf; it=1;
  while distance>tol  && it<maxit
    #2.1 Solve value function backwards, given the guess prices. (This needs to be updated to correct prices)
    maxnorm=0.0;

    for dummy=1:(T-1)
      t = T + 1 - dummy;
      # Firms
      pr= deepcopy(tr[t].fpr);
      firmbellmanParallelOmega!(pr,tr[t].eq,tausec[t],pa,maxroutine);
      getpolicies!(pr,tr[t].eq,tausec[t],pa);
      tr[t-1].fpr = deepcopy(pr);
    println(  norm(pr.firmvaluegrid - ssT.fpr.firmvaluegrid)  )
      #Households
      ## householdProblem!(tr[t-1].hpr, aprime::tr[t].hpr.a,aprime::tr[t].hpr.mu, tr[t].eq.r, tr[t].eq.w, tausec[t], pa);
    end

    #2.2 Compute cross-sectional distribution forward (this can be made faster if needed)
    # and aggregate firm problem (I think I want to first check the free entry condition)


    # Check if value at 1 is zero. Then set E to zero or 1.

    # Compute distribution.
    aux = transitionrule(ss0.eq.distr,tr[1].eq.E, ss0.fpr, tr[1].eq , tausec[1], pa);
    tr[1].eq.distr = deepcopy(aux);
    tr[1].eq.a.laborsupply, tr[1].eq.a.consumption = firm_aggregates_transitions(tr[1].fpr, ss0.fpr, tr[1].eq, ss0.eq, tausec[1], pa);

    # Wage consistent with household problem
    what[1] = pa.H*(tr[1].eq.a.laborsupply)^1/pa.psi;
    # Compute contribution to price discrepancy
    maxnorm = max(abs(what[1] - tr[1].eq.w),maxnorm);
    # Update wage (note that updating at this point is fine)
    tr[1].eq.w = update*tr[1].eq.w + (1-update)*what[1];

    rhat[1] = ss0.eq.r; #since interest rate at time t is r_{t+1} we dont care/use rhat[1]

    for t=2:T
      #Compute firm value, set E to zero or 1

      aux = transitionrule(tr[t-1].eq.distr,tr[t].eq.E, tr[t-1].fpr, tr[t].eq , tausec[t], pa);
      tr[t].eq.distr = deepcopy(aux);
      tr[t].eq.a.laborsupply, tr[t].eq.a.consumption = firm_aggregates_transitions(tr[t].fpr, tr[t-1].fpr, tr[t].eq, tr[t-1].eq, tausec[t], pa);
      # Compute candidate wage
      what[t] = pa.H*(tr[t].eq.a.laborsupply)^(1.0/pa.psi);
      # Compute contribution to price discrepancy
      maxnorm = max(abs(what[t] - tr[t].eq.w),maxnorm);
      # Update wage (note that updating at this point is fine)
      tr[t].eq.w = update*tr[t].eq.w + (1.0-update)*what[t];

      Cprime = tr[t].eq.a.consumption;
      consprime[t]=Cprime;
      C = tr[t-1].eq.a.consumption;
      cons[t]=C;
      Lprime = tr[t].eq.a.laborsupply;
      L = tr[t-1].eq.a.laborsupply;
      rhat[t] = (1.0-tausec[t].i)^(-1.0)* (   pa.beta^(-1.0) *  ( (Cprime - pa.H/(1.0+1.0/pa.sigma)*Lprime)/(C - pa.H/(1.0+1.0/pa.sigma)*L) )^pa.sigma  -1.0   );
      # Compute contribution to price discrepancy
      maxnorm =max( abs(rhat[t] - tr[t].eq.r), maxnorm);
      # Update interest rate
      tr[t].eq.r = update*tr[t].eq.r + (1-update)*rhat[t];
    end

    distance = maxnorm;
    println("it = ", it, " dist = ", distance)
    println(rhat)
    it+=1;
  end

end



##############################################


##############################################



function firm_aggregates_transitions(crt_fpr::FirmProblem, prev_fpr::FirmProblem, crt_eq::Equilibrium, prev_eq::Equilibrium, crt_tau::Taxes, pa::Param)
# Computes labor demand and net supply of goods at time t. Where time t is crt/prime variables.

  prev_distr = prev_eq.distr;
  crt_distr  = crt_eq.distr;
  wage  = crt_eq.w;
  irate = crt_eq.r;

  # Compute aggregates depending on current K
  lprime_d = 0.0;
  gdp = 0.0;
  corptax = 0.0;
  investment_check = 0.0;

  for i_zprime in 1:pa.Nz
    zprime= pa.zgrid[i_zprime];
    # Entrants into current period
    kprime = 0.0;
    qprime = 0.0;
    mass = crt_eq.E*pa.invariant_distr[i_zprime];

    lprime    = (zprime*pa.alphal*kprime^pa.alphak / wage)^(1/(1-pa.alphal));
    lprime_d += lprime*mass;
    gdp      += (zprime*kprime^pa.alphak*lprime^pa.alphal - pa.f)*mass;
    corptax  += tau.c*(zprime*kprime^pa.alphak*lprime^pa.alphal - wage*lprime - pa.delta*kprime - irate*qprime - pa.f)*mass;

    investment_check+=crt_fpr.kpolicy[1,i_zprime]*mass;

    for i_z in 1:pa.Nz
      probzprime=pa.ztrans[i_zprime,i_z];

      # Next consider incumbent that had state omega{t-1} z_{t-1} previous period
      for i_omega in 1:pa.Nomega
        kprime = prev_fpr.kpolicy[i_omega,i_z];
        qprime = prev_fpr.qpolicy[i_omega,i_z];
        mass = prev_distr[i_omega,i_z];

        if !prev_fpr.exitrule[i_omega,i_z,i_zprime]
          lprime    = (zprime*pa.alphal*kprime^pa.alphak / wage)^(1/(1-pa.alphal));
          lprime_d += lprime*mass*probzprime;
          gdp      += (zprime*kprime^pa.alphak*lprime^pa.alphal - pa.f)*mass*probzprime;
          corptax  += tau.c*(zprime*kprime^pa.alphak*lprime^pa.alphal - wage*lprime - pa.delta*kprime - irate*qprime - pa.f)*mass*probzprime;

          omegaprime, i_omegaprime = predict_state(i_zprime, i_omega, i_z, prev_fpr, crt_eq, crt_tau, pa);
          investment_check+=crt_fpr.kpolicy[i_omegaprime,i_zprime]*mass*probzprime;
        else
          corptax  += tau.c*( - pa.delta*kprime - irate*qprime )*mass*probzprime;
        end
      end
    end
  end

  # Compute aggregates that do not depend on K
  capital= sum(prev_distr.*prev_fpr.kpolicy);
  investment = sum(crt_distr.*crt_fpr.kpolicy) - (1-pa.delta)*capital;
  if abs(sum(crt_distr.*crt_fpr.kpolicy) - investment_check)>10.0^-10.0
    error(" Investment doesn't add up by ",sum(crt_distr.*crt_fpr.kpolicy) - investment_check )
  end
  grossdividends = sum(crt_distr.*crt_fpr.grossdividends);
  financialcosts= - sum(crt_distr.*crt_fpr.financialcosts);
  debt = sum(crt_distr.*crt_fpr.qpolicy);
  divtax= tau.d*grossdividends;
  inctax = tau.i*eq.r*debt;
  G = divtax + corptax + inctax;

  goodsnet_s = gdp - G - financialcosts - investment;
  return lprime_d,goodsnet_s;
end
