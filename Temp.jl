function unit_entry(distr1::Matrix, res::ResultsFP, pr::FirmProblem, p::Equilibrium, tau::Taxes, fp::FirmParam)
  # This function uses the distribution for a UNITARY mass of entrant to compute
  # the aggregates needed to compute the actual mss of entrants

  labor_d1 bonds1 netdistributions1 liquidations1

  debt=0.0;
  labor=0.0;
  netdistributions=0.0;
  liquidations=0.0;

  for i_z in 1:pr.Nz
    #### First we consider entrants ###

    netdistributions+=res.distributions[1,i_z]*fp.invariant_distr[i_z]

    #### Next we consider incumbents ###
    for i_omega in 1:pr.Nomega
      kprime=res.kprime[i_omega,i_z];
      qprime=res.qprime[i_omega,i_z];


      for i_zprime in 1:pr.Nz
        #Non-exiting incumbents
        if !res.exitrule[i_omega,i_z,i_zprime]
          zprime       = fp.zgrid[i_zprime];
          lprime       = (zprime*fp.alphal*kprime^fp.alphak / p.w)^(1/(1-fp.alphal));

          debt        += qprime*distr1[i_omega,i_z]*fp.ztrans[i_zprime,i_z];
          labor       += lprime*distr1[i_omega,i_z]*fp.ztrans[i_zprime,i_z];

          omegaprime, i_omegaprime = predict_state(i_zprime, i_omega, i_z, pr, tau, fp)
          netdistributions+=res.distributions[i_omegaprime,i_zprime]*distr1[i_omega,i_z]*fp.ztrans[i_zprime,i_z];

        #Exiting incumbents
        else
          liquidations += (1-pr.taudtilde)*(fp.kappa*(1-fp.delta)*kprime - (1+p.r)*qprime)*distr1[i_omega,i_z]*fp.ztrans[i_zprime,i_z];
          debt += qprime*distr1[i_omega,i_z]*fp.ztrans[i_zprime,i_z];

        end
      end

    end

  end

  netdistributionscheck = sum(distr1.*res.distributions);
  debtcheck = sum(distr1.*res.qprime)

  if debt - debtcheck >10.0^-4.0 ||
      netdistributions - netdistributionscheck >10.0^-4.0
    error("Consistency problem")
  end


  return debt, labor, netdistributions, liquidations
end
