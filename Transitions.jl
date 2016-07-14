# This script takes an initial equilibrium, a final equilibrium and a path for
# taxes and return a path for equilibrium prices and allocations


function init_transitions(T::Int64)
  T=30;
  ss0=PeriodSolution(pr,eq);
  ssT=PeriodSolution(pr2,eq2);
  tausec=Array(Taxes,(T,));
  for t=1:T
  	tausec[t]=tau2;
  end
  tr= Array(PeriodSolution,(T,));

  return ss0, blah
end

function transitions!(tr:: Array(PeriodSolution,(T,)), T::Int64, ss0::PeriodSolution, ssT::PeriodSolution, tauseq::Array{Taxes,1}, pa::Param; update::Float64=0.9, tol::Float64 = 10^-4.0, maxit=10^4.0::Float64 , maxroutine::Function=maximizationstep )


    # Guess paths for r, w, E
    for t=1:T-1
      wguess= t*(ssT.eq.w - ss0.eq.w)/T + ss0.eq.w;
      rguess= (t-1)*(ssT.eq.r - ss0.eq.r)/T + ss0.eq.r;
      println("w = ",wguess, " r = ", rguess)
      eq = init_equilibirium(wguess,tausec[t],pa; r=rguess);
      fpr = init_firmproblem(eq,tausec[t],pa);
      tr[t] = PeriodSolution(fpr,eq);
      tr[t].eq.E= t*(ssT.eq.E - ss0.eq.E)/T + ss0.eq.E;
    end

    # Assume new steady state is reached after T periods
    tr[T] = ssT;


    #Fix a path for entry and solve for prices



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
  crt_distr= crt_eq.distr;
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
          corptax  += tau.c*(zprime*kprime^pa.alphak*lprime^pa.alphal - wage*lprime - pa.delta*kprime - irate*qprime - pa.f)*mass*probzprime;
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

function updateprices!(tr:: Array(PeriodSolution,(T,)), T::Int64, ss0::PeriodSolution; verbose=true, update=0.9)

  what=Array(Float64,(T,));
  rhat=Array(Float64,(T,));

  distance=+Inf;
  it=1;
  while distance>tol  && it<maxit

    #Solve value function backwards
    norm=0.0;

    for dummy=1:(T-1)
      t = T + 1 - dummy
      println(t)
      # Firms
      pr= deepcopy(tr[t].fpr);
      firmbellmanParallelOmega!(pr,tr[t].eq,tausec[t],pa,maxroutine);
      getpolicies!(pr,tr[t].eq,tausec[t],pa);
      tr[t-1].fpr = deepcopy(pr);
      #Households
      ## householdProblem!(tr[t-1].hpr, aprime::tr[t].hpr.a,aprime::tr[t].hpr.mu, tr[t].eq.r, tr[t].eq.w, tausec[t], pa);
    end

    #Compute cross-sectional distribution forward (this can be made faster if needed)
    # and aggregate firm problem

    tr[1].eq.distr = ss0.eq.distr;
    tr[1].eq.a.laborsupply, tr[1].eq.a.consumption = firm_aggregates_transitions(tr[1].fpr, ss0.fpr, tr[1].eq, ss0.eq, tausec[1], pa);
    what[1] = pa.H*(tr[1].eq.a.laborsupply)^1/pa.psi;
    # Compute contribution to price discrepancy
    norm += max(abs(what[1] - tr[1].eq.w),norm);
    # Update wage (note that updating at this point is fine)
    tr[1].eq.distr = ss0.eq.distr;
    tr[1].eq.a.laborsupply, tr[1].eq.a.consumption = firm_aggregates_transitions(tr[1].fpr, ss0.fpr, tr[1].eq, ss0.eq, tausec[1], pa);
    what[1] = pa.H*(tr[1].eq.a.laborsupply)^1/pa.psi;
    # Compute contribution to price discrepancy
    norm += max(abs(what[1] - tr[1].eq.w),norm);
    # Update wage (note that updating at this point is fine)
    tr[1].eq.w = update*tr[1].eq.w + (1-update)*what[1];

    rhat[1] = ss0.eq.r #since interest rate at time t is r_{t+1} we dont care/use rhat[1]

    for t=2:T
      tr[t].eq.distr = transitionrule(tr[t-1].eq.distr,tr[t].eq.E, tr[t-1].fpr, tr[t].eq , tausec[t], pa);
      tr[t].eq.a.laborsupply, tr[t].eq.a.consumption = firm_aggregates_transitions(tr[t].fpr, tr[t-1].fpr, tr[t].eq, tr[t-1].eq, tausec[t], pa);

      # Compute candidate wage
      what[t] = pa.H*(tr[t].eq.a.laborsupply)^(1.0/pa.psi);
      # Compute contribution to price discrepancy
      norm += max(abs(what[t] - tr[t].eq.w),norm);
      # Update wage (note that updating at this point is fine)
      tr[t].eq.w = update*tr[t].eq.w + (1.0-update)*what[t];

      Cprime = tr[t].eq.a.consumption;
      C = tr[t-1].eq.a.consumption;
      Lprime = tr[t].eq.a.laborsupply;
      L = tr[t-1].eq.a.laborsupply;
      rhat[t] = (1.0-tausec[t].i)^(-1.0)* (   pa.beta^(-1.0) *  ( (Cprime - pa.H/(1.0+1.0/pa.sigma)*Lprime)/(C - pa.H/(1.0+1.0/pa.sigma)*L) )^pa.sigma  -1.0   );
      # Compute contribution to price discrepancy
      norm =max( abs(rhat[t] - tr[t].eq.r), norm);
      # Update interest rate
      tr[t].eq.r = update*tr[t].eq.r + (1-update)*rhat[t];
    end



    distance = norm;
    it+=1;
  end
end
