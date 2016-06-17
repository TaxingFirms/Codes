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

          omegaprime, i_omegaprime = predict_state(i_zprime, i_omega, i_z, p, pr, tau, fp)
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

  E=( hp.H*(labor_d1 + p.w^(-1.0)*((1-tau.i)*p.r*bonds1+ netdistributions1 +liquidations1) ) )^-1.0;
  p.distr = E*distr1;
  p.E = E;
end



################################################################################
################################################################################



function aggregates!(res::ResultsFP, pr::FirmProblem, p::Equilibrium, tau::Taxes, hp::HouseholdParam, fp::FirmParam)

  distr=p.distr;


  ########## Initialize accumulators ###########
  # Aggregates
  capital=0.0;
  debt=0.0;
  labor=0.0;
  gdp=0.0;
  corptax=0.0;
  liquidations=0.0;
  liquidationcosts=0.0;
  netdistributions=0.0;

  #Moments
  mean_inv_rate=0.0;
  var_inv_rate=0.0;
  mean_leverage=0.0;
  var_leverage=0.0;

  mean_dividends2k=0.0;
  var_dividends2k=0.0;
  mean_profits2k=0.0;

  var_profits2k=0.0;
  mean_eqis2k=0.0;
  freq_eqis2k=0.0;
  mean_tobinsq=0.0;

  cov_nw=0.0;

  # Auxiliary variables
  mean_inv_rate_shifted=0.0;
  mean_leverage_shifted=0.0;
  mean_dividends2k_shifted=0.0;
  mean_profits2k_shifted=0.0;
  mean_omega_shifted=0.0;
  mean_omegaprime_shifted=0.0;
  mass_incumbents=0.0;


  ##Define constants to avoid catastrophic cancellation
  # (Understanding this isn't important)
  mp_omega=convert(Int64,round(pr.Nomega/2));
  mp_z=convert(Int64,round(pr.Nz/2));
  mp_zprime=mp_z;
  vmp_omegaprime, mp_omegaprime = predict_state(mp_zprime, mp_omega, mp_z, p, pr, tau, fp);
  kprime = res.kprime[mp_omega,mp_z];
  qprime = res.qprime[mp_omega,mp_z];
  zprime = fp.zgrid[mp_zprime];
  lprime = (zprime*fp.alphal*kprime^fp.alphak / p.w)^(1/(1-fp.alphal));

  Kinv = ((res.kprime[mp_omegaprime,mp_zprime] - (1-fp.delta)*kprime)/kprime);
  Klev = qprime/kprime;
  Kdiv= res.grossdividends[mp_omegaprime,mp_zprime]/kprime;
  Kprof = (zprime*kprime^fp.alphak*lprime^fp.alphal - p.w*lprime - fp.f)/kprime;
  Komega = vmp_omegaprime;
  ##########################################



  ###### COMPUTE VALUES THAT DEPEND ON CURRENT K
  for i_z in 1:pr.Nz
    #### First we consider entrants ###
    netdistributions+=res.distributions[1,i_z]*p.E*fp.invariant_distr[i_z];
    #investment rate, leverage, etc are undefinied for entrants

    #### Next we consider incumbents ###
    for i_omega in 1:pr.Nomega
      kprime=res.kprime[i_omega,i_z];
      qprime=res.qprime[i_omega,i_z];
      omega = pr.omega.grid[i_omega];


      for i_zprime in 1:pr.Nz
        #Non-exiting incumbents
        if !res.exitrule[i_omega,i_z,i_zprime]
          zprime       = fp.zgrid[i_zprime];
          lprime       = (zprime*fp.alphal*kprime^fp.alphak / p.w)^(1/(1-fp.alphal));

          capital     += kprime*distr[i_omega,i_z]*fp.ztrans[i_zprime,i_z];
          debt        += qprime*distr[i_omega,i_z]*fp.ztrans[i_zprime,i_z];
          labor       += lprime*distr[i_omega,i_z]*fp.ztrans[i_zprime,i_z];

          gdp         += (zprime*kprime^fp.alphak*lprime^fp.alphal -fp.f)*distr[i_omega,i_z]*fp.ztrans[i_zprime,i_z];
          corptax     += tau.c*(zprime*kprime^fp.alphak*lprime^fp.alphal -p.w*lprime - fp.delta*kprime - p.r*qprime -fp.f)*distr[i_omega,i_z]*fp.ztrans[i_zprime,i_z];

          omegaprime, i_omegaprime = predict_state(i_zprime, i_omega, i_z, p, pr, tau, fp);
          netdistributions+=res.distributions[i_omegaprime,i_zprime]*distr[i_omega,i_z]*fp.ztrans[i_zprime,i_z];

          if kprime>0
            inv_rate = ((res.kprime[i_omegaprime,i_zprime] - (1-fp.delta)*kprime)/kprime);
            mean_inv_rate += inv_rate*distr[i_omega,i_z]*fp.ztrans[i_zprime,i_z];
            mean_inv_rate_shifted += (inv_rate- Kinv)*distr[i_omega,i_z]*fp.ztrans[i_zprime,i_z];
            var_inv_rate += (inv_rate-Kinv)^2*distr[i_omega,i_z]*fp.ztrans[i_zprime,i_z];

            leverage = qprime/kprime;
            mean_leverage += leverage*distr[i_omega,i_z]*fp.ztrans[i_zprime,i_z];
            mean_leverage_shifted += (leverage - Klev)*distr[i_omega,i_z]*fp.ztrans[i_zprime,i_z];
            var_leverage += (leverage - Klev)^2*distr[i_omega,i_z]*fp.ztrans[i_zprime,i_z];

            div2k= res.grossdividends[i_omegaprime,i_zprime]/kprime; #Before tax dividends
            mean_dividends2k += div2k*distr[i_omega,i_z]*fp.ztrans[i_zprime,i_z];
            mean_dividends2k_shifted += (div2k - Kdiv)*distr[i_omega,i_z]*fp.ztrans[i_zprime,i_z];
            var_dividends2k += (div2k - Kdiv)^2*distr[i_omega,i_z]*fp.ztrans[i_zprime,i_z];

            prof2k = (zprime*kprime^fp.alphak*lprime^fp.alphal -p.w*lprime -fp.f)/kprime;
            mean_profits2k += prof2k*distr[i_omega,i_z]*fp.ztrans[i_zprime,i_z];
            mean_profits2k_shifted += (prof2k - Kprof)*distr[i_omega,i_z]*fp.ztrans[i_zprime,i_z];
            var_profits2k +=  (prof2k - Kprof)^2*distr[i_omega,i_z]*fp.ztrans[i_zprime,i_z];

            eqis2k = res.grossequityis[i_omegaprime,i_zprime]/kprime;
            mean_eqis2k += eqis2k*distr[i_omega,i_z]*fp.ztrans[i_zprime,i_z];

            tobinsq = res.firmvalue[i_omegaprime,i_zprime]/kprime;
            mean_tobinsq += tobinsq*distr[i_omega,i_z]*fp.ztrans[i_zprime,i_z];

            cov_nw += (omegaprime - Komega)*(omega-Komega)*distr[i_omega,i_z]*fp.ztrans[i_zprime,i_z];
            mean_omega_shifted += (omega - Komega);
            mean_omegaprime_shifted += (omegaprime - Komega);
          end
          mass_incumbents += distr[i_omega,i_z]*fp.ztrans[i_zprime,i_z];

        #Exiting incumbents
        else
          liquidationcosts += (1-fp.kappa)*(1-fp.delta)*kprime*distr[i_omega,i_z]*fp.ztrans[i_zprime,i_z]
          liquidations += (1-pr.taudtilde)*(fp.kappa*(1-fp.delta)*kprime - (1+p.r)*qprime)*distr[i_omega,i_z]*fp.ztrans[i_zprime,i_z];
          capital += kprime*distr[i_omega,i_z]*fp.ztrans[i_zprime,i_z];
          debt += qprime*distr[i_omega,i_z]*fp.ztrans[i_zprime,i_z];
        end
      end

    end

  end
  #######################################




  ### COMPUTE VALUES THAT DO NOT DEPEND K
  investment = sum(distr.*res.kprime) - (1-fp.delta)*capital;
  grossdividends=sum(distr.*res.grossdividends);
  financialcosts= - sum(distr.*res.financialcosts);

  divtax= tau.d*grossdividends;
  inctax = tau.i*p.r*debt;

  G = divtax + corptax + inctax;
  ######################################



  ######## Consistency Checks #########
  netdistributionscheck = sum(distr.*res.distributions);
  debtcheck = sum(distr.*res.qprime)
  capitalcheck= sum(distr.*res.kprime)
  labor_s = 1/hp.H - 1/p.w*((1-tau.i)*p.r*debt + netdistributions +liquidations);
  consumption = p.w/hp.H;

  if capital - capitalcheck >10.0^-4.0 ||
      debt - debtcheck >10.0^-4.0 ||
      netdistributions - netdistributionscheck >10.0^-4.0
    error("Consistency problem")
  end

  if abs(labor_s - labor) > 10^-10.0
    error("labor market didn't clear")
  end

  if abs((gdp - G - financialcosts - investment -liquidationcosts  ) - consumption)> pr.omega.step
        println("goods market didn't clear")
  end

  ##############################



  ########## Moments ##########
  mean_inv_rate = mean_inv_rate / mass_incumbents;
  var_inv_rate= (var_inv_rate - mean_inv_rate_shifted^2/mass_incumbents )/ mass_incumbents;

  mean_leverage=mean_leverage/ mass_incumbents;
  var_leverage=(var_leverage - mean_leverage_shifted^2/mass_incumbents)/ mass_incumbents;

  mean_dividends2k= mean_dividends2k/ mass_incumbents;
  var_dividends2k=(var_dividends2k - mean_dividends2k_shifted^2/mass_incumbents)/ mass_incumbents;

  mean_profits2k=mean_profits2k/ mass_incumbents;
  var_profits2k=(var_profits2k - mean_profits2k_shifted^2/mass_incumbents )/ mass_incumbents;

  mean_eqis2k= mean_eqis2k/ mass_incumbents;

  freq_eqis2k=freq_eqis2k/ mass_incumbents;
  mean_tobinsq=mean_tobinsq/ mass_incumbents;

  cov_nw=(cov_nw - mean_omega_shifted*mean_omegaprime_shifted / mass_incumbents)/ mass_incumbents;

  p.m = Moments(mean_inv_rate, var_inv_rate, mean_leverage, var_leverage, mean_dividends2k, var_dividends2k, mean_profits2k, var_profits2k, mean_eqis2k, freq_eqis2k, mean_tobinsq, cov_nw)
  ##############################



  ######## Save values #########
  welfare = (log(consumption) - hp.H* labor_s)/(1-hp.beta);
  collections = Taxes(divtax,corptax,inctax,0);

 p.a=Aggregates(netdistributions, (1-tau.i)*p.r*debt, consumption, gdp, labor, welfare, collections, debt, capital, investment, grossdividends, financialcosts, G);


end
