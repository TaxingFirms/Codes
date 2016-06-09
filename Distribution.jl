
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
function stationarydist(E::Real,res::ResultsFP, pr::FirmProblem, p::Equilibrium, tau::Taxes, fp::FirmParam;  tol=10.0^-4.0)
  #Initiate Variables
  ii=zeros(Int,pr.Nomega*pr.Nz*pr.Nz);
  jj=zeros(Int,pr.Nomega*pr.Nz*pr.Nz);
  vv=zeros(pr.Nomega*pr.Nz*pr.Nz);
  #Construct the transition matrix
  k=1;
  for i_z=1:pr.Nz
    for i_omega = 1:pr.Nomega
      kprime= pr.kpolicygrid[i_omega,i_z];
      qprime= pr.qpolicygrid[i_omega,i_z];
      for i_zprime= 1:fp.Nz
        if !res.exitrule[i_omega, i_z,i_zprime]
          omegaprime = omegaprimefun(kprime,qprime,i_zprime,p,tau,fp)
          i_omegaprime = closestindex(omegaprime, pr.omega.step)
          if i_omegaprime<1 || i_omegaprime>pr.Nomega
            i_omega == pr.Nomega?
              i_omegaprime =pr.Nomega:
              error("omega' out of the grid ", " i_z = ", i_z, " i_omega = ", i_omega, " i_z' = ",i_zprime, " i_omega' = ",i_omegaprime)
          end
          ii[k]= find_dist_ind(i_omega, i_z, pr.Nomega);
          jj[k]= find_dist_ind(i_omegaprime, i_zprime, pr.Nomega);
          vv[k]= fp.ztrans[i_zprime,i_z];
          k+=1;
        end
      end
    end
  end
  #The initial vector was larger than needed because of exit. We trim it.
  rows=ii[1:k-1];
  cols=jj[1:k-1];
  vals=vv[1:k-1];
  transmat = sparse(rows,cols,vals,pr.Nz*pr.Nomega,pr.Nz*pr.Nomega);

  #Entrants
  ie=zeros(Int64,pr.Nz);
  ve=zeros(Float64,pr.Nz);
  k=1;
  for i_z in 1:fp.Nz
        ie[k]= find_dist_ind(1, i_z, pr.Nomega);
        ve[k]= E*fp.invariant_distr[i_z];
        k+=1;
  end

  entrants= sparsevec(ie,ve,pr.Nz*pr.Nomega);
  I=speye(Float64,pr.Nz*pr.Nomega); #Identity matrix

  distrguess = ones(pr.Nomega* pr.Nz);
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
  distr=reshape(distr_vectorized, (pr.Nomega,pr.Nz))

end

function stationarydist_iterate( E, res, pr, p, tau, fp; tol=eps() )
  distrguess = ones(pr.Nomega, pr.Nz);
  distrguess=  distrguess/sum(distrguess);

  distr=deepcopy(distrguess);
  dist=Inf;
  it=1; maxit = 10.0^5.0;
  while dist >tol && it<maxit
    println("it = ",it," dist = ", dist)
    distrprime = transitionrule(distr,E,res,pr,p,tau,fp)
    dist= norm(distrprime -distr);

    distr=deepcopy(distrprime);
    it+=1;
  end
  if it==maxit
    error("Distribution did not converge after maximum number of iterations")
  end

  distr
end



function transitionrule(distr::Matrix,E::Real,res::ResultsFP, pr::FirmProblem, p::Equilibrium, tau:: Taxes, fp::FirmParam)

  distrprime= zeros(pr.Nomega, pr.Nz);

  for i_z in 1:pr.Nz
    #First we consider entrants
    distrprime[1,i_z]+=E*fp.invariant_distr[i_z];
    #Non-exiting incumbents
    for i_omega in 1:pr.Nomega
      kprime=res.kprime[i_omega,i_z];
      qprime=res.qprime[i_omega,i_z];
      for i_zprime in 1:pr.Nz
        if !res.exitrule[i_omega,i_z,i_zprime]
          omegaprime = omegaprimefun(kprime,qprime,i_zprime,p,tau,fp)
          i_omegaprime = closestindex(omegaprime, pr.omega.step);

          #The block below checks that the index is within reasonable bounds
          if i_omegaprime<1 || i_omegaprime>pr.Nomega
            i_omega == pr.Nomega?
              i_omegaprime =pr.Nomega:
              error("omega' out of the grid ", "i_z = ", i_z, "i_omega = ", i_omega, "i_z' = ",i_zprime, "i_omega' = ",i_omegaprime)
          end
#=          if i_omegaprime<1 || i_omegaprime>pr.Nomega
            if i_omega==pr.Nomega || i_omegaprime < (pr.Nomega + 3)
              i_omegaprime =pr.Nomega
            elseif i_omega==1 || i_omegaprime > -3
              i_omegaprime =1;
            else
              error("omega' out of the grid ", "i_z = ", i_z, "i_omega = ", i_omega, "i_z' = ",i_zprime, "i_omega' = ",i_omegaprime)
            end
          end
=#
          distrprime[i_omegaprime,i_zprime] += fp.ztrans[i_zprime,i_z]*distr[i_omega,i_z];
        end

      end
    end
  end
distrprime
end

function distributionStupid(E::Real,res::ResultsFP, pr::FirmProblem, p::Equilibrium, tau::Taxes, fp::FirmParam;  tol=10.0^-4.0)

  initialDistr = zeros(Float64,pr.Nomega,pr.Nz)
  nextDistr = zeros(Float64,pr.Nomega,pr.Nz)

  initialDistr[1,:] = fp.invariant_distr
  diffValue = 10.0
  while diffValue > tol
    for i_z=1:pr.Nz, i_omega = 1:pr.Nomega
      kprime= pr.kpolicygrid[i_omega,i_z];
      qprime= pr.qpolicygrid[i_omega,i_z];
      for i_zprime= 1:fp.Nz
        if !res.exitrule[i_omega, i_z,i_zprime]
          omegaprime = omegaprimefun(kprime,qprime,i_zprime,p,tau,fp)
          i_omegaprime = closestindex(omegaprime, pr.omega.step)
          if i_omegaprime>pr.Nomega  && i_omega == pr.Nomega
            i_omegaprime = pr.Nomega
          elseif i_omegaprime < 1
              error("omega' out of the grid ", " i_z = ", i_z, " i_omega = ", i_omega, " i_z' = ",i_zprime, " i_omega' = ",i_omegaprime)
          end
          nextDistr[i_omegaprime,i_zprime] += fp.ztrans[i_zprime,i_z]*initialDistr[i_omega,i_z]
        end
      end # end i_z'

      if i_omega == 1
        nextDistr[i_omega,i_z] += E*fp.invariant_distr[i_z]
      end
    end # end i_omega, i_z

    diffValue    = maxabs(nextDistr - initialDistr)
    initialDistr = deepcopy(nextDistr)
    nextDistr    = zeros(Float64,pr.Nomega,pr.Nz)

  end

  initialDistr

end

###### CONSTANT EXOGENOUS EXIT
function stationarydistCE(E::Real, pr::FirmProblem, p::Equilibrium, tau::Taxes, fp::FirmParam;  tol=10.0^-4.0)
  #Initiate Variables
  ii=zeros(Int,pr.Nomega*pr.Nz*pr.Nz);
  jj=zeros(Int,pr.Nomega*pr.Nz*pr.Nz);
  vv=zeros(pr.Nomega*pr.Nz*pr.Nz);
  #Construct the transition matrix
  k=1;
  for i_z=1:pr.Nz
    for i_omega = 1:pr.Nomega
      kprime= pr.kpolicygrid[i_omega,i_z];
      qprime= pr.qpolicygrid[i_omega,i_z];
      for i_zprime= 1:fp.Nz
#        if !res.exitrule[i_omega, i_z,i_zprime]
          omegaprime = omegaprimefun(kprime,qprime,i_zprime,p,tau,fp)
          i_omegaprime = closestindex(omegaprime, pr.omega.step)
          if i_omegaprime<1
            i_omegaprime=1;
          elseif i_omegaprime>pr.Nomega
            i_omegaprime =pr.Nomega
          end
          ii[k]= find_dist_ind(i_omega, i_z, pr.Nomega);
          jj[k]= find_dist_ind(i_omegaprime, i_zprime, pr.Nomega);
          vv[k]= fp.ztrans[i_zprime,i_z]*(1-exitprob);
          k+=1;
#        end

      end
    end
  end
  #The initial vector was larger than needed because of exit. We trim it.

  transmat = sparse(ii,jj,vv,pr.Nz*pr.Nomega,pr.Nz*pr.Nomega);

  #Entrants
  ie=zeros(Int64,pr.Nz);
  ve=zeros(Float64,pr.Nz);
  k=1;
  for i_z in 1:fp.Nz
        ie[k]= find_dist_ind(1, i_z, pr.Nomega);
        ve[k]= E*fp.invariant_distr[i_z];
        k+=1;
  end

  entrants= sparsevec(ie,ve,pr.Nz*pr.Nomega);
  I=speye(Float64,pr.Nz*pr.Nomega); #Identity matrix

  distrguess = ones(pr.Nomega* pr.Nz);
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
  distr=reshape(distr_vectorized, (pr.Nomega,pr.Nz))

end
