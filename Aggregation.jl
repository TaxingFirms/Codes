function unit_entry(distr1::Matrix, res::ResultsFP, pr::FirmProblem, p::Equilibrium, tau::Taxes, fp::FirmParam)
  # This function uses the distribution for a UNITARY mass of entrant to compute
  # the aggregates needed to compute the actual mss of entrants

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


################################################################################
################################################################################

function mass_of_entrants!( res::ResultsFP, pr::FirmProblem, p::Equilibrium, tau::Taxes, fp::FirmParam)
  #Computes the mass of entrants such that the labor market clears,
  distr1=stationarydist(1,res, pr, p, tau,fp);
  bonds1, labor_d1, netdistributions1, liquidations1 = unit_entry(distr1::Matrix, res::ResultsFP, pr::FirmProblem, p::Equilibrium, tau::Taxes, fp::FirmParam);
  p.E=( hp.H*(labor_d1 + p.w^(-1.0)*((1-tau.i)*p.r*bonds1+ netdistributions1 +liquidations1) ) )^-1.0;
end



function equilibrium_aggregates!( res::ResultsFP, pr::FirmProblem, p::Equilibrium, tau::Taxes, fp::FirmParam)

  capital1, bonds1, labor_d1, gdp1, corptax1, inctax1, netdistributions1, liquidations1, liquidationcosts1, m1 = aggregates(1, distr1, res, pr, p, tau, fp);
  investment1 = sum(distr1.*res.kprime) - (1-fp.delta)*capital1;
  grossdividends1=sum(distr1.*res.grossdividends);
  divtax1= tau.d*grossdividends1;
  financialcosts1= - sum(distr1.*res.financialcosts);
  G1 = divtax1 + corptax1 + inctax1;

  consumption = p.w/hp.H;
  p.consumption=consumption;
#  E= consumption/(gdp1 - G1 - financialcosts1 - investment1)
  p.distr= E*distr1;
  println("Mass = ", sum(E*distr1))
  p.E=E;
  p.netdistributions=E* netdistributions1;
  p.agginterests= E*(1-tau.i)*p.r*bonds1;
  p.collections=Taxes(E*divtax1,E*corptax1,E*inctax1,0);
  p.output= gdp1*E;

  labor_s = 1/hp.H - 1/p.w*E*((1-tau.i)*p.r*bonds1 + netdistributions1 +liquidations1);
  p.laborsupply=labor_s;

  p.welfare= log(consumption) - hp.H* labor_s;

  #Consisitency checks
  if abs(labor_s - E*labor_d1) > eps()
    error("labor market didn't clear")
  end

  if abs((gdp1 - G1 - financialcosts1 - investment1 -liquidationcosts1  )*E - consumption)> pr.omega.step
        println("goods market didn't clear")
  end

  # Moment update


end
