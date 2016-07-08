# This script takes an initial equilibrium, a final equilibrium and a path for
# taxes and return a path for equilibrium prices and allocations

type SteadyState
  fpr::FirmProblem
  hpr::HouseholdProblem
  eq::Equilibrium
end

type HouseholdProblem
  c::Float64
  l::Float64
  a::Float64
  mu::Float64
end

function init_transitions(T::Int64)
  transitions= Array(SteadyState,(T,));
end

function transitions!(tr:: Array(SteadyState,(T,)), T::Int64, ss0::SteadyState, ssT::SteadyState , tauseq::Array{Taxes,1}, pa::Param; update::Float64=0.9, tol = 10^-4.0, maxit=10^4.0 , maxroutine::Function=maximizationstep )

    # Guess paths for r, w, E
    for t=1:T-1
      tr[t].eq.w= t*(ssT.eq.w - ss0.eq.w)/T + ss0.eq.w;
      tr[t].eq.r= t*(ssT.eq.r - ss0.eq.r)/T + ss0.eq.r;
      tr[t].eq.E= t*(ssT.eq.E - ss0.eq.E)/T + ss0.eq.E;
    end

    # Assume new steady state is reached after T periods
    tr[T] = ssT;



    # Solve for path for prices using a shooting algorithm
    while excess>tol  && it<maxit

      #Solve value function backwards
      for t=T:2
        # Firms
        pr= deepcopy(tr[t].fpr);
        firmbellmanParallelOmega!(pr,tr[t].eq,tausec[t],pa,maxroutine);
        getpolicies!(pr,tr[t].eq,tr[t].tau,pa);
        tr[t-1].fpr = pr;

        #Households
        householdProblem!(tr[t-1].hpr, aprime::tr[t].hpr.a,aprime::tr[t].hpr.mu, tr[t].eq.r, tr[t].eq.w, tausec[t], pa);
      end

      #Compute cross-sectional distribution forward (this can be made faster if needed)
      tr[1].eq.distr = transitionrule(ss0.eq.distr, tr[1].eq.E, tr[1].fpr, tr[1].eq , tausec[1], pa);
      for t=2:T
        tr[t].eq.distr = transitionrule(tr[t-1].eq.distr,tr[t].eq.E, tr[t].fpr, tr[t].eq , tausec[t], pa);
      end

      #Check market clearing conditions
#Need to write an aggregation file for labor demand (maybe within transition rule) but bonds and distributions
#can be computed directly
        #Asset market
        #Labor market
        #Free Entry

      #Update prices

    end


end



function householdProblem!(hpr::HouseholdProblem, aprime::Float64, mgutilityprime::Float64, r::Float64, w::Float64, tau::Taxes, pa::Param)
  labor= (w/pa.H)^pa.psi;
  dsc = (1+ (1-tau.i)*r)*pa.beta;
  hpr.a= (mgutilityprime/dsc)^(-1/pa.sigma) - w*labor + aprime + (pa.H/(1+1/pa.psi))*labor^(1+pa.psi);
  hpr.c= assets +  w*labor - aprime/dsc;
  hpr.mu= (consumption - (pa.H/(1+1/pa.psi))*labor^(1+pa.psi))^(-pa.sigma);
  hpr.l=labor;
end
