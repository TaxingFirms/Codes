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
