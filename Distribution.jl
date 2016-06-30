# This script has the functions needed to update the stationary distribution,
# given policy functions from pr:FirmProblem and a mass of entrants eq.E
# getpolicies! needs to be run before any of the functions to extract the exit
# rule
#
# There are three functions to do this (and all should return the same
# distribution): stationarydist, stationarydist_iterate and distributionStupid



#Find the closest gridpoint to omega'
function closestindex(omegaprime::Real, step::Real)
  aproxindex=omegaprime/step;
  round(Int,aproxindex) +1
end

#Convert the indexes from z and omega into an index from the distribution
#The inner variable is always omega
function find_dist_ind(i_omega::Int, i_z::Int, Nomega::Int)
  Nomega*(i_z-1)+i_omega
end

#Convert the index of the transition matrix into the corresponding indexes from z and omega
function find_inds(i::Int, Nomega::Int)
  if mod(i,Nomega)==0
    i_omega= Nomega;
    i_z= i/Nomega;
  else
    i_omega = mod(i,Nomega);
    i_z= (i-i_omega)/Nomega +1;
  end
  (convert(Int,i_omega), convert(Int,i_z))
end

# STATIONARY DISTRIBUTION,  E is the mass of entrants
function stationarydist(E::Real, pr::FirmProblem, eq::Equilibrium, tau::Taxes, pa::Param;  tol=eps())
  #Initiate Variables
  ii=zeros(Int,pa.Nomega*pa.Nz*pa.Nz);
  jj=zeros(Int,pa.Nomega*pa.Nz*pa.Nz);
  vv=zeros(pa.Nomega*pa.Nz*pa.Nz);
  #Construct the transition matrix
  k=1;
  for i_z=1:pa.Nz
    for i_omega = 1:pa.Nomega
      kprime= pr.kpolicy[i_omega,i_z];
      qprime= pr.qpolicy[i_omega,i_z];
      for i_zprime= 1:pa.Nz
        if !pr.exitrule[i_omega, i_z,i_zprime]
          omegaprime = omegaprimefun(kprime, qprime, i_zprime, eq, tau, pa)
          i_omegaprime = closestindex(omegaprime, pa.omega.step)
          if i_omegaprime<1 || i_omegaprime>pa.Nomega
            i_omega == pa.Nomega?
              i_omegaprime =pa.Nomega:
              error("omega' out of the grid ", " i_z = ", i_z, " i_omega = ", i_omega, " i_z' = ",i_zprime, " i_omega' = ",i_omegaprime)
          end
          ii[k]= find_dist_ind(i_omega, i_z, pa.Nomega);
          jj[k]= find_dist_ind(i_omegaprime, i_zprime, pa.Nomega);
          vv[k]= pa.ztrans[i_zprime,i_z];
          k+=1;
        end
      end
    end
  end
  #The initial vector was larger than needed because of exit. We trim it.
  rows=ii[1:k-1];
  cols=jj[1:k-1];
  vals=vv[1:k-1];
  transmat = sparse(rows,cols,vals,pa.Nz*pa.Nomega,pa.Nz*pa.Nomega);

  #Entrants
  ie=zeros(Int64,pa.Nz);
  ve=zeros(Float64,pa.Nz);
  k=1;
  for i_z in 1:pa.Nz
        ie[k]= find_dist_ind(1, i_z, pa.Nomega);
        ve[k]= E*pa.invariant_distr[i_z];
        k+=1;
  end

  entrants= sparsevec(ie,ve,pa.Nz*pa.Nomega);
  I=speye(Float64,pa.Nz*pa.Nomega); #Identity matrix

  distrguess = ones(pa.Nomega* pa.Nz);
  distrguess=  distrguess/sum(distrguess);

  distr_vectorized=deepcopy(distrguess);
  dist=Inf;
  it=1; maxit = 10.0^5.0;
  while dist >tol && it<maxit
    println("it = ",it," dist = ", dist)
    #println(it)
    distrprime= transmat'*distr_vectorized +entrants;

    dist= norm(distrprime -distr_vectorized);
    distr_vectorized=deepcopy(distrprime);
    it+=1;
  end
  if it==maxit
    error("Distribution did not converge after maximum number of iterations")
  end


  #=
  A= full(I-transmat');
  b=full(entrants);
  det(A) != 0 ?
    dist_vectorized = \(A,b):
    error("I - T is not invertible");
  =#
  distr=reshape(distr_vectorized, (pa.Nomega,pa.Nz))

end

function stationarydist_iterate( E::Float64, pr::FirmProblem, eq::Equilibrium, tau::Taxes, pa::Param; tol=eps() )
  distrguess = ones(pa.Nomega, pa.Nz);
  distrguess=  distrguess/sum(distrguess);

  distr=deepcopy(distrguess);
  dist=Inf;
  it=1; maxit = 10.0^5.0;
  while dist >tol && it<maxit
    println("it = ",it," dist = ", dist)
    distrprime = transitionrule(distr,E,pr,eq,tau,pa)
    dist= norm(distrprime -distr);

    distr=deepcopy(distrprime);
    it+=1;
  end
  if it==maxit
    error("Distribution did not converge after maximum number of iterations")
  end

  distr
end



function transitionrule(distr::Matrix,E::Real, pr::FirmProblem, eq::Equilibrium, tau:: Taxes, pa::Param)

  distrprime= zeros(pa.Nomega, pa.Nz);

  for i_z in 1:pa.Nz
    #First we consider entrants
    distrprime[1,i_z]+=E*pa.invariant_distr[i_z];
    #Non-exiting incumbents
    for i_omega in 1:pa.Nomega
      kprime=pr.kpolicy[i_omega,i_z];
      qprime=pr.qpolicy[i_omega,i_z];
      for i_zprime in 1:pa.Nz
        if !pr.exitrule[i_omega,i_z,i_zprime]
          omegaprime = omegaprimefun(kprime, qprime, i_zprime, eq, tau, pa)
          i_omegaprime = closestindex(omegaprime, pa.omega.step);

          #The block below checks that the index is within reasonable bounds
          if i_omegaprime<1 || i_omegaprime>pa.Nomega
            i_omega == pa.Nomega?
              i_omegaprime =pa.Nomega:
              error("omega' out of the grid ", "i_z = ", i_z, "i_omega = ", i_omega, "i_z' = ",i_zprime, "i_omega' = ",i_omegaprime)
          end
#=          if i_omegaprime<1 || i_omegaprime>pa.Nomega
            if i_omega==pa.Nomega || i_omegaprime < (pa.Nomega + 3)
              i_omegaprime =pa.Nomega
            elseif i_omega==1 || i_omegaprime > -3
              i_omegaprime =1;
            else
              error("omega' out of the grid ", "i_z = ", i_z, "i_omega = ", i_omega, "i_z' = ",i_zprime, "i_omega' = ",i_omegaprime)
            end
          end
=#
          distrprime[i_omegaprime,i_zprime] += pa.ztrans[i_zprime,i_z]*distr[i_omega,i_z];
        end

      end
    end
  end
distrprime
end

function distributionStupid(E::Real,pr::FirmProblem, eq::Equilibrium, tau::Taxes, pa::Param;  tol=10.0^-4.0)

  initialDistr = zeros(Float64,pa.Nomega,pa.Nz)
  nextDistr = zeros(Float64,pa.Nomega,pa.Nz)

  initialDistr[1,:] = pa.invariant_distr
  diffValue = 10.0
  while diffValue > tol
    for i_z=1:pa.Nz, i_omega = 1:pa.Nomega
      kprime= pr.kpolicy[i_omega,i_z];
      qprime= pr.qpolicy[i_omega,i_z];
      for i_zprime= 1:pa.Nz
        if !pr.exitrule[i_omega, i_z,i_zprime]
          omegaprime = omegaprimefun(kprime, qprime, i_zprime, eq, tau, pa)
          i_omegaprime = closestindex(omegaprime, pa.omega.step)
          if i_omegaprime>pa.Nomega  && i_omega == pa.Nomega
            i_omegaprime = pa.Nomega
          elseif i_omegaprime < 1
              error("omega' out of the grid ", " i_z = ", i_z, " i_omega = ", i_omega, " i_z' = ",i_zprime, " i_omega' = ",i_omegaprime)
          end
          nextDistr[i_omegaprime,i_zprime] += pa.ztrans[i_zprime,i_z]*initialDistr[i_omega,i_z]
        end
      end # end i_z'

      if i_omega == 1
        nextDistr[i_omega,i_z] += E*pa.invariant_distr[i_z]
      end
    end # end i_omega, i_z

    diffValue    = maxabs(nextDistr - initialDistr)
    initialDistr = deepcopy(nextDistr)
    nextDistr    = zeros(Float64,pa.Nomega,pa.Nz)

  end

  initialDistr

end
