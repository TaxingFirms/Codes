# Given the equilibirum wage and the equilibrium distribution consistent with
# such wage (both in eq::Equilibrium), this script computes the mass of entrants
# consistent with market clearing and saves it in eq.E. It also computes
# aggregates and moments and saves them in eq.a and eq.m


function unit_entry(distr1::Matrix, pr::FirmProblem, eq::Equilibrium, tau::Taxes, pa::Param)
  # This function uses the distribution for a UNITARY mass of entrant to compute
  # the aggregates needed to compute the actual mss of entrants
  taudtilde = 1-(1-tau.d)/(1-tau.g);
  debt=0.0;
  labor=0.0;
  netdistributions=0.0;
  liquidations=0.0;

  for i_z in 1:pa.Nz
    #### First we consider entrants ###

    netdistributions+=pr.distributions[1,i_z]*pa.invariant_distr[i_z]

    #### Next we consider incumbents ###
    for i_omega in 1:pa.Nomega
      kprime=pr.kpolicy[i_omega,i_z];
      qprime=pr.qpolicy[i_omega,i_z];


      for i_zprime in 1:pa.Nz
        #Non-exiting incumbents
        if !pr.exitrule[i_omega,i_z,i_zprime]
          zprime       = pa.zgrid[i_zprime];
          lprime       = (zprime*pa.alphal*kprime^pa.alphak / eq.w)^(1/(1-pa.alphal));

          debt        += qprime*distr1[i_omega,i_z]*pa.ztrans[i_zprime,i_z];
          labor       += lprime*distr1[i_omega,i_z]*pa.ztrans[i_zprime,i_z];

          omegaprime, i_omegaprime = predict_state(i_zprime, i_omega, i_z, pr, eq, tau, pa)
          netdistributions+=pr.distributions[i_omegaprime,i_zprime]*distr1[i_omega,i_z]*pa.ztrans[i_zprime,i_z];

        #Exiting incumbents
        else
          liquidations += (1-taudtilde)*(pa.kappa*(1-pa.delta)*kprime - (1+eq.r)*qprime)*distr1[i_omega,i_z]*pa.ztrans[i_zprime,i_z];
          debt += qprime*distr1[i_omega,i_z]*pa.ztrans[i_zprime,i_z];
        end
      end

    end

  end

  netdistributionscheck = sum(distr1.*pr.distributions);
  debtcheck = sum(distr1.*pr.qpolicy)

  if abs(debt - debtcheck) >10.0^-4.0 ||
      abs(netdistributions - netdistributionscheck) >10.0^-4.0
    error("Consistency problem")
  end


  return debt, labor, netdistributions, liquidations
end



################################################################################
################################################################################



function mass_of_entrants!( pr::FirmProblem, eq::Equilibrium, tau::Taxes, pa::Param, distribution::Function; verbose=true)
  #Computes the mass of entrants such that the labor market clears,
  verbose?distr1=distribution(1.0,pr, eq, tau,pa):distr1=distribution(1.0,pr, eq, tau,pa;verbose=false);
  bonds1, labor_d1, netdistributions1, liquidations1 = unit_entry(distr1, pr, eq, tau, pa);

  E=( pa.H*(labor_d1 + eq.w^(-1.0)*((1-tau.i)*eq.r*bonds1+ netdistributions1 +liquidations1) ) )^-1.0;
  eq.distr = E*distr1;
  eq.E = E;
end


function mass_of_entrantsGHH!( pr::FirmProblem, eq::Equilibrium, tau::Taxes, pa::Param, distribution::Function; verbose=true)
  #Computes the mass of entrants such that the labor market clears,
  verbose?distr1=distribution(1.0,pr, eq, tau,pa):distr1=distribution(1.0,pr, eq, tau,pa;verbose=false);
  bonds1, labor_d1, netdistributions1, liquidations1 = unit_entry(distr1, pr, eq, tau, pa);
  labor_s= (eq.w/pa.H)^pa.psi

  E= labor_s/labor_d1;
  eq.distr = E*distr1;
  eq.E = E;
end


################################################################################
################################################################################



function aggregates!(pr::FirmProblem, eq::Equilibrium, tau::Taxes, pa::Param; compute_moments=true)
# Computes aggregates and moments once the model is completely solved.
  distr=eq.distr;
  taudtilde = 1-(1-tau.d)/(1-tau.g);

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
  if compute_moments
    mean_inv_rate=0.0;
    var_inv_rate=0.0;
    mean_leverage=0.0;
    var_leverage=0.0;
    mean_dividends2k=0.0;
    var_dividends2k=0.0;
    mean_profits2k=0.0;
    var_profits2k=0.0;
    mean_eqis2k=0.0;
    mean_eqis=0.0;
    freq_eqis=0.0;
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
    mp_omega=convert(Int64,round(pa.Nomega/2));
    mp_z=convert(Int64,round(pa.Nz/2));
    mp_zprime=mp_z;
    vmp_omegaprime, mp_omegaprime = predict_state(mp_zprime, mp_omega, mp_z, pr, eq, tau, pa);
    kprime = pr.kpolicy[mp_omega,mp_z];
    qprime = pr.qpolicy[mp_omega,mp_z];
    zprime = pa.zgrid[mp_zprime];
    lprime = (zprime*pa.alphal*kprime^pa.alphak / eq.w)^(1/(1-pa.alphal));

    Kinv = ((pr.kpolicy[mp_omegaprime,mp_zprime] - (1-pa.delta)*kprime)/kprime);
    Klev = qprime/kprime;
    Kdiv= pr.grossdividends[mp_omegaprime,mp_zprime]/kprime;
    Kprof = (zprime*kprime^pa.alphak*lprime^pa.alphal - eq.w*lprime - pa.f)/kprime;
    Komega = vmp_omegaprime;
    ##########################################
  end



  ###### COMPUTE VALUES THAT DEPEND ON CURRENT K
  for i_z in 1:pa.Nz
    #### First we consider entrants ###
    netdistributions+=pr.distributions[1,i_z]*eq.E*pa.invariant_distr[i_z];
    kprime=0.0;
    qprime=0.0;
    lprime    = (zprime*pa.alphal*kprime^pa.alphak / wage)^(1/(1-pa.alphal));
    lprime_d += lprime*eq.E*pa.invariant_distr[i_z];
    gdp      += (zprime*kprime^pa.alphak*lprime^pa.alphal - pa.f)*eq.E*pa.invariant_distr[i_z];
    corptax  += tau.c*(zprime*kprime^pa.alphak*lprime^pa.alphal - wage*lprime - pa.delta*kprime - irate*qprime - pa.f)*eq.E*pa.invariant_distr[i_z];

    if compute_moments
      freq_eqis+=(1-pr.positivedistributions[1,i_z])*eq.E*pa.invariant_distr[i_z];
      mean_eqis+=pr.grossequityis*eq.E*pa.invariant_distr[i_z];
      #investment rate, leverage, etc are undefinied for entrants
    end

    #### Next we consider incumbents ###
    for i_omega in 1:pa.Nomega
      kprime=pr.kpolicy[i_omega,i_z];
      qprime=pr.qpolicy[i_omega,i_z];
      omega = pa.omega.grid[i_omega];


      for i_zprime in 1:pa.Nz
        #Non-exiting incumbents
        if !pr.exitrule[i_omega,i_z,i_zprime]
          zprime       = pa.zgrid[i_zprime];
          lprime       = (zprime*pa.alphal*kprime^pa.alphak / eq.w)^(1/(1-pa.alphal));

          capital     += kprime*distr[i_omega,i_z]*pa.ztrans[i_zprime,i_z];
          debt        += qprime*distr[i_omega,i_z]*pa.ztrans[i_zprime,i_z];
          labor       += lprime*distr[i_omega,i_z]*pa.ztrans[i_zprime,i_z];

          gdp         += (zprime*kprime^pa.alphak*lprime^pa.alphal - pa.f)*distr[i_omega,i_z]*pa.ztrans[i_zprime,i_z];
          corptax     += tau.c*(zprime*kprime^pa.alphak*lprime^pa.alphal -eq.w*lprime - pa.delta*kprime - eq.r*qprime - pa.f)*distr[i_omega,i_z]*pa.ztrans[i_zprime,i_z];

          omegaprime, i_omegaprime = predict_state(i_zprime, i_omega, i_z, pr, eq, tau, pa);
          netdistributions+=pr.distributions[i_omegaprime,i_zprime]*distr[i_omega,i_z]*pa.ztrans[i_zprime,i_z];

          if compute_moments
            freq_eqis+=(1-pr.positivedistributions[i_omegaprime,i_zprime])*distr[i_omega,i_z]*pa.ztrans[i_zprime,i_z];
            mean_eqis+=pr.grossequityis[i_omegaprime,i_zprime]*distr[i_omega,i_z]*pa.ztrans[i_zprime,i_z];

            if kprime>0
              inv_rate = ((pr.kpolicy[i_omegaprime,i_zprime] - (1-pa.delta)*kprime)/kprime);
              mean_inv_rate += inv_rate*distr[i_omega,i_z]*pa.ztrans[i_zprime,i_z];
              mean_inv_rate_shifted += (inv_rate- Kinv)*distr[i_omega,i_z]*pa.ztrans[i_zprime,i_z];
              var_inv_rate += (inv_rate-Kinv)^2*distr[i_omega,i_z]*pa.ztrans[i_zprime,i_z];

              leverage = qprime/kprime;
              mean_leverage += leverage*distr[i_omega,i_z]*pa.ztrans[i_zprime,i_z];
              mean_leverage_shifted += (leverage - Klev)*distr[i_omega,i_z]*pa.ztrans[i_zprime,i_z];
              var_leverage += (leverage - Klev)^2*distr[i_omega,i_z]*pa.ztrans[i_zprime,i_z];

              div2k= pr.grossdividends[i_omegaprime,i_zprime]/kprime; #Before tax dividends
              mean_dividends2k += div2k*distr[i_omega,i_z]*pa.ztrans[i_zprime,i_z];
              mean_dividends2k_shifted += (div2k - Kdiv)*distr[i_omega,i_z]*pa.ztrans[i_zprime,i_z];
              var_dividends2k += (div2k - Kdiv)^2*distr[i_omega,i_z]*pa.ztrans[i_zprime,i_z];

              prof2k = (zprime*kprime^pa.alphak*lprime^pa.alphal -eq.w*lprime -pa.f)/kprime;
              mean_profits2k += prof2k*distr[i_omega,i_z]*pa.ztrans[i_zprime,i_z];
              mean_profits2k_shifted += (prof2k - Kprof)*distr[i_omega,i_z]*pa.ztrans[i_zprime,i_z];
              var_profits2k +=  (prof2k - Kprof)^2*distr[i_omega,i_z]*pa.ztrans[i_zprime,i_z];

              eqis2k = pr.grossequityis[i_omegaprime,i_zprime]/kprime;
              mean_eqis2k += eqis2k*distr[i_omega,i_z]*pa.ztrans[i_zprime,i_z];

              tobinsq = pr.firmvaluegrid[i_omegaprime,i_zprime]/kprime;
              mean_tobinsq += tobinsq*distr[i_omega,i_z]*pa.ztrans[i_zprime,i_z];

              cov_nw += (omegaprime - Komega)*(omega-Komega)*distr[i_omega,i_z]*pa.ztrans[i_zprime,i_z];
              mean_omega_shifted += (omega - Komega);
              mean_omegaprime_shifted += (omegaprime - Komega);
            end
          end
          mass_incumbents += distr[i_omega,i_z]*pa.ztrans[i_zprime,i_z];

        #Exiting incumbents
        else
          liquidationcosts += (1-pa.kappa)*(1-pa.delta)*kprime*distr[i_omega,i_z]*pa.ztrans[i_zprime,i_z]
          liquidations += (1-taudtilde)*(pa.kappa*(1-pa.delta)*kprime - (1+eq.r)*qprime)*distr[i_omega,i_z]*pa.ztrans[i_zprime,i_z];
          capital += kprime*distr[i_omega,i_z]*pa.ztrans[i_zprime,i_z];
          debt += qprime*distr[i_omega,i_z]*pa.ztrans[i_zprime,i_z];
        end
      end

    end

  end
  #######################################




  ### COMPUTE VALUES THAT DO NOT DEPEND K
  investment = sum(distr.*pr.kpolicy) - (1-pa.delta)*capital;
  grossdividends=sum(distr.*pr.grossdividends);
  financialcosts= - sum(distr.*pr.financialcosts);

  divtax= tau.d*grossdividends;
  inctax = tau.i*eq.r*debt;

  G = divtax + corptax + inctax;
  ######################################



  ######## Consistency Checks #########
  netdistributionscheck = sum(distr.*pr.distributions);
  debtcheck = sum(distr.*pr.qpolicy)
  capitalcheck= sum(distr.*pr.kpolicy)
  labor_s = (eq.w/pa.H)^pa.psi;
  consumption = eq.w*labor_s +  (1-tau.i)*eq.r*debt+ netdistributions +liquidations ;

  if abs(capital - capitalcheck) >10.0^-4.0 ||
      abs(debt - debtcheck) >10.0^-4.0 ||
      abs(netdistributions - netdistributionscheck) >10.0^-4.0
    error("Consistency problem")
  end

  if abs(labor_s - labor) > 10^-10.0
    error("labor market didn't clear")
  end

  if abs((gdp - G - financialcosts - investment -liquidationcosts  ) - consumption)> pa.omega.step
        println("goods market didn't clear by: ", (gdp - G - financialcosts - investment -liquidationcosts  ) - consumption)
  end

  ##############################



  ########## Moments ##########
  if compute_moments
    mean_inv_rate = mean_inv_rate / mass_incumbents;
    var_inv_rate= (var_inv_rate - mean_inv_rate_shifted^2/mass_incumbents )/ mass_incumbents;

    mean_leverage=mean_leverage/ mass_incumbents;
    var_leverage=(var_leverage - mean_leverage_shifted^2/mass_incumbents)/ mass_incumbents;

    mean_dividends2k= mean_dividends2k/ mass_incumbents;
    var_dividends2k=(var_dividends2k - mean_dividends2k_shifted^2/mass_incumbents)/ mass_incumbents;

    mean_profits2k=mean_profits2k/ mass_incumbents;
    var_profits2k=(var_profits2k - mean_profits2k_shifted^2/mass_incumbents )/ mass_incumbents;

    mean_eqis2k= mean_eqis2k/ capital;

    freq_eqis=freq_eqis/ mass_incumbents;
    mean_tobinsq=mean_tobinsq/ mass_incumbents;

    cov_nw=(cov_nw - mean_omega_shifted*mean_omegaprime_shifted / mass_incumbents)/ mass_incumbents;

    eq.m = Moments(mean_inv_rate, var_inv_rate, mean_leverage, var_leverage, mean_dividends2k, var_dividends2k, mean_profits2k, var_profits2k, mean_eqis2k, freq_eqis, mean_tobinsq, cov_nw)
  end
  ##############################




  ######## Save values #########
  welfare = (log(consumption) - pa.H* labor_s)/(1-pa.beta);
  collections = Taxes(divtax,corptax,inctax,0);

  eq.a=Aggregates(netdistributions, (1-tau.i)*eq.r*debt, consumption, gdp, labor, welfare, collections, debt, capital, investment, grossdividends, financialcosts, G);


end
