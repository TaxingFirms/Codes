function aggregatesCE( exitprob::Float64, E::Real, distr::Matrix, res::ResultsFP, pr::FirmProblem, p::Equilibrium, tau::Taxes, fp::FirmParam)
  capital=0.0;
  debt=0.0;
  labor=0.0;
  gdp=0.0;
  corptax=0.0;
  netdistributions=0.0;
  liquidations=0.0;
  liquidationcosts=0.0;


  for i_z in 1:pr.Nz
    #First we consider entrants
#= These are all zero
    zprime =fp.zgrid[i_z];
    kprime=0;
    qprime=0;
    lprime = (zprime*fp.alphal*kprime^fp.alphak / p.w)^(1/(1-fp.alphal));

    capital += kprime*E*fp.invariant_distr[i_z];
    debt += qprime*E*fp.invariant_distr[i_z];
    labor += labor*E*fp.invariant_distr[i_z];
    gdp += (zprime*kprime^fp.alphak*lprime^fp.alphal -fp.f)*E*fp.invariant_distr[i_z];
    corptax += tau.c*(zprime*kprime^fp.alphak*lprime^fp.alphal - p.w*lprime - fp.delta*kprime - p.r*qprime )*E*fp.invariant_distr[i_z];
=#
    netdistributions+=res.distributions[1,i_z]*E*fp.invariant_distr[i_z]

    for i_omega in 1:pr.Nomega
      kprime=res.kprime[i_omega,i_z];
      qprime=res.qprime[i_omega,i_z];


      for i_zprime in 1:pr.Nz
        #Non-exiting incumbents
        #if !res.exitrule[i_omega,i_z,i_zprime] > 0
          zprime       = fp.zgrid[i_zprime];
          lprime       = (zprime*fp.alphal*kprime^fp.alphak / p.w)^(1/(1-fp.alphal));

          capital     += kprime*distr[i_omega,i_z]*fp.ztrans[i_zprime,i_z]*(1-exitprob);
          debt        += qprime*distr[i_omega,i_z]*fp.ztrans[i_zprime,i_z]*(1-exitprob);
          labor       += lprime*distr[i_omega,i_z]*fp.ztrans[i_zprime,i_z]*(1-exitprob);

          gdp         += (zprime*kprime^fp.alphak*lprime^fp.alphal -fp.f)*distr[i_omega,i_z]*fp.ztrans[i_zprime,i_z]*(1-exitprob);
          corptax     += tau.c*(zprime*kprime^fp.alphak*lprime^fp.alphal -p.w*lprime - fp.delta*kprime - p.r*qprime -fp.f)*distr[i_omega,i_z]*fp.ztrans[i_zprime,i_z]*(1-exitprob);

          omegaprime   = omegaprimefun(kprime,qprime,i_zprime,p,tau,fp)
          i_omegaprime = closestindex(omegaprime, pr.omega.step);
          #The block below checks that the index is within reasonable bounds
          if i_omegaprime<1 || i_omegaprime>pr.Nomega
            if i_omega==pr.Nomega || i_omegaprime < (pr.Nomega + 3)
              i_omegaprime =pr.Nomega
            elseif i_omega==1 || i_omegaprime > -3
              i_omegaprime =1;
            else
              error("omega' out of the grid ", "i_z = ", i_z, "i_omega = ", i_omega, "i_z' = ",i_zprime, "i_omega' = ",i_omegaprime)
            end
          end
          netdistributions+=res.distributions[i_omegaprime,i_zprime]*distr[i_omega,i_z]*fp.ztrans[i_zprime,i_z]*(1-exitprob);

        #Exiting incumbents
        #else
          liquidationcosts += (1-fp.kappa)*(1-fp.delta)*kprime*distr[i_omega,i_z]*fp.ztrans[i_zprime,i_z]*exitprob
          liquidations += (fp.kappa*(1-fp.delta)*kprime - (1+p.r)*qprime)*distr[i_omega,i_z]*fp.ztrans[i_zprime,i_z]*exitprob;
          capital += kprime*distr[i_omega,i_z]*fp.ztrans[i_zprime,i_z]*exitprob;
          debt += qprime*distr[i_omega,i_z]*fp.ztrans[i_zprime,i_z]*exitprob;
        #end
      end

    end

  end

  inctax = tau.i*p.r*debt;
  netdistributionscheck = sum(distr.*res.distributions);
  debtcheck = sum(distr.*res.qprime)
  capitalcheck= sum(distr.*res.kprime)

  if capital - capitalcheck >10.0^-4.0 ||
      debt - debtcheck >10.0^-4.0 ||
      netdistributions - netdistributionscheck >10.0^-4.0
    error("Consistency problem")
  end

  return capital, debt, labor, gdp, corptax, inctax, netdistributions, liquidations, liquidationcosts
end


###### CONSTANT EXOGENOUS EXIT
function stationarydistCE(E::Real, pr::FirmProblem, p::Equilibrium, tau::Taxes, fp::FirmParam;  tol=10.0^-4.0)
  #Initiate Variables
  ii=zeros(Int,pa.Nomega*pa.Nz*pa.Nz);
  jj=zeros(Int,pa.Nomega*pa.Nz*pa.Nz);
  vv=zeros(pa.Nomega*pa.Nz*pa.Nz);
  #Construct the transition matrix
  k=1;
  for i_z=1:pa.Nz
    for i_omega = 1:pa.Nomega
      kprime= pr.kpolicygrid[i_omega,i_z];
      qprime= pr.qpolicygrid[i_omega,i_z];
      for i_zprime= 1:fp.Nz
#        if !res.exitrule[i_omega, i_z,i_zprime]
          omegaprime = omegaprimefun(kprime,qprime,i_zprime,p,tau,fp)
          i_omegaprime = closestindex(omegaprime, pr.omega.step)
          if i_omegaprime<1
            i_omegaprime=1;
          elseif i_omegaprime>pa.Nomega
            i_omegaprime =pa.Nomega
          end
          ii[k]= find_dist_ind(i_omega, i_z, pa.Nomega);
          jj[k]= find_dist_ind(i_omegaprime, i_zprime, pa.Nomega);
          vv[k]= fp.ztrans[i_zprime,i_z]*(1-exitprob);
          k+=1;
#        end

      end
    end
  end
  #The initial vector was larger than needed because of exit. We trim it.

  transmat = sparse(ii,jj,vv,pa.Nz*pa.Nomega,pa.Nz*pa.Nomega);

  #Entrants
  ie=zeros(Int64,pa.Nz);
  ve=zeros(Float64,pa.Nz);
  k=1;
  for i_z in 1:fp.Nz
        ie[k]= find_dist_ind(1, i_z, pa.Nomega);
        ve[k]= E*fp.invariant_distr[i_z];
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
