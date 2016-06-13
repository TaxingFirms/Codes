
function aggregates(E::Real, distr::Matrix, res::ResultsFP, pr::FirmProblem, p::Equilibrium, tau::Taxes, fp::FirmParam)
  capital=0.0;
  debt=0.0;
  labor=0.0;
  gdp=0.0;
  corptax=0.0;
  netdistributions=0.0;
  liquidations=0.0;
  liquidationcosts=0.0;

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
  freq_equis2k=0.0;
  mean_tobinsq=0.0;

  cov_nw=0.0;

  # Auxiliary variables
  mean_inv_rate_shifted=0.0;
  mean_leverage_shifted=0.0;
  mean_dividends2k_shifted=0.0;
  mean_profits2k_shifted=0.0;
  mean_omega_shifted=0.0;
  mass_incumbents=0.0;


  ##Define constants to avoid catastrophic cancellation
  # (Understanding this isn't important)
  mp_omega=round(pr.Nomega/2);
  mp_z=round(pr.Nz/2);
  mp_zprime=mp_z;
  vmp_omegaprime, mp_omegaprime = predict_state(mp_zprime, mp_omega, mp_z, pr, tau, fp);
  kprime = res.kprime[mp_omega,mp_z];
  qprime = res.qprime[mp_omega,mp_z];
  zprime = fp.zgrid[mp_zprime];
  lprime = (zprime*fp.alphal*kprime^fp.alphak / p.w)^(1/(1-fp.alphal));

  Kinv = ((res.kprime[mp_omegaprime,mp_zprime] - (1-fp.delta)*kprime)/kprime);
  Klev = qprime/kprime;
  Kdiv= res.grossdividends[mp_omegaprime,mp_zprime]/kprime;
  Kprof = (zprime*kprime^fp.alphak*lprime^fp.alphal - p.w*lprime - fp.f)/kprime;
  Komega = vmp_omegaprime;
  #




  for i_z in 1:pr.Nz
    #### First we consider entrants ###


    #= These are all zero
    zprime =fp.zgrid[i_z];
    kprime=0;
    qprime=0;
    lprime = (zprime*fp.alphal*kprime^fp.alphak / p.w)^(1/(1-fp.alphal));

    capital += kprime*E*fp.invariant_distr[i_z];
    # inv_rate, leverage, payout_rate are undefined for entrantrants.
    debt += qprime*E*fp.invariant_distr[i_z];
    labor += labor*E*fp.invariant_distr[i_z];
    gdp += (zprime*kprime^fp.alphak*lprime^fp.alphal -fp.f)*E*fp.invariant_distr[i_z];
    corptax += tau.c*(zprime*kprime^fp.alphak*lprime^fp.alphal - p.w*lprime - fp.delta*kprime - p.r*qprime )*E*fp.invariant_distr[i_z];
    =#
    netdistributions+=res.distributions[1,i_z]*E*fp.invariant_distr[i_z]

    #investment rate, leverage, etc are undefinied for entrants

    #### Next we consider incumbents ###
    for i_omega in 1:pr.Nomega
      kprime=res.kprime[i_omega,i_z];
      qprime=res.qprime[i_omega,i_z];


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

          omegaprime, i_omegaprime = predict_state(i_zprime, i_omega, i_z, pr, tau, fp)

          inv_rate = ((res.kprime[i_omegaprime,i_zprime] - (1-fp.delta)*kprime)/kprime);
          mean_inv_rate += inv_rate*distr[i_omega,i_z]*fp.ztrans[i_zprime,i_z];
          mean_inv_rate_shifted += (inv_rate- Kinv)*distr[i_omega,i_z]*fp.ztrans[i_zprime,i_z];
          var_inv_rate += (inv_rate-Kinv)^2*distr[i_omega,i_z]*fp.ztrans[i_zprime,i_z];

          leverage = qprime/kprime;
          mean_leverage += leverage*distr[i_omega,i_z]*fp.ztrans[i_zprime,i_z];
          mean_leverage_shifted += (leverage - Klev)*distr[i_omega,i_z]*fp.ztrans[i_zprime,i_z];
          var_leverage += (leverage - Klev)^2*distr[i_omega,i_z]*fp.ztrans[i_zprime,i_z];

          netdistributions+=res.distributions[i_omegaprime,i_zprime]*distr[i_omega,i_z]*fp.ztrans[i_zprime,i_z];

          div2k= res.grossdividends[i_omegaprime,i_zprime]/kprime; #Before tax dividends
          mean_dividends2k += div2k*distr[i_omega,i_z]*fp.ztrans[i_zprime,i_z];
          mean_dividends2k_shifted += (div2k - Kdiv)*distr[i_omega,i_z]*fp.ztrans[i_zprime,i_z];
          var_dividends2k += (div2k - Kdiv)^2*distr[i_omega,i_z]*fp.ztrans[i_zprime,i_z];

          prof2k = (zprime*kprime^fp.alphak*lprime^fp.alphal -p.w*lprime -fp.f)/kprime;
          mean_profits2k += prof2k*distr[i_omega,i_z]*fp.ztrans[i_zprime,i_z];
          mean_profits2k_shifted += (prof2k - Kprof)*distr[i_omega,i_z]*fp.ztrans[i_zprime,i_z];
          var_profits2k +=  (prof2k - Kprof)^2*distr[i_omega,i_z]*fp.ztrans[i_zprime,i_z];

          eqis2k = res.grossequityis[i_omegaprime,i_zprime]/kprime;
          mean_eqis2k += equis2k*distr[i_omega,i_z]*fp.ztrans[i_zprime,i_z];

          tobinsq = res.firmvalue[i_omegaprime,i_zprime]/kprime;
          mean_tobinsq += tobinsq*distr[i_omega,i_z]*fp.ztrans[i_zprime,i_z];

          cov_nw += (omegaprime - Komega)*(omega-Komega)*distr[i_omega,i_z]*fp.ztrans[i_zprime,i_z];
          mean_omega_shifted += (omega - Komega);
          mean_omegaprime_shifted += (omegaprime - Komega);

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

  inctax = tau.i*p.r*debt;
  netdistributionscheck = sum(distr.*res.distributions);
  debtcheck = sum(distr.*res.qprime)
  capitalcheck= sum(distr.*res.kprime)

  if capital - capitalcheck >10.0^-4.0 ||
      debt - debtcheck >10.0^-4.0 ||
      netdistributions - netdistributionscheck >10.0^-4.0
    error("Consistency problem")
  end

  #Moments
  mean_inv_rate = mean_inv_rate / mass_incumbents;
  var_inv_rate= (var_inv_rate - mean_inv_rate_shifted^2/mass_incumbents )/ mass_incumbents;

  mean_leverage=mean_leverage/ mass_incumbents;
  var_leverage=(var_leverage - mean_leverage_shifted^2/mass_incumbents)/ mass_incumbents;

  mean_dividends2k= mean_dividends2k/ mass_incumbents;
  var_dividends2k=(var_dividends2k - mean_dividends2k_shifted^2/mass_incumbents)/ mass_incumbents;

  mean_profits2k=mean_profits2k/ mass_incumbents;
  var_profits2k=(var_profits2k - mean_profits2k_shifted^2/mass_incumbents )/ mass_incumbents;

  mean_eqis2k= mean_eqis2k/ mass_incumbents;

  freq_equis2k=freq_equis2k/ mass_incumbents;
  mean_tobinsq=mean_tobinsq/ mass_incumbents;

  cov_nw=(cov_nw - mean_omega_shifted*mean_omegaprime_shifted / mass_incumbents)/ mass_incumbents;


  return capital, debt, labor, gdp, corptax, inctax, netdistributions, liquidations, liquidationcosts
end



function equilibrium_aggregates!( res::ResultsFP, pr::FirmProblem, p::Equilibrium, tau::Taxes, fp::FirmParam)
  #Computes the mass of entrants such that the labor market clears,
  distr1=stationarydist(1,res, pr, p, tau,fp);
  capital1, bonds1, labor_d1, gdp1, corptax1, inctax1, netdistributions1, liquidations1, liquidationcosts1 = aggregates(1, distr1, res, pr, p, tau, fp);
  investment1 = sum(distr1.*res.kprime) - (1-fp.delta)*capital1;
  grossdividends1=sum(distr1.*res.grossdividends);
  divtax1= tau.d*grossdividends1;
  financialcosts1= - sum(distr1.*res.financialcosts);
  G1 = divtax1 + corptax1 + inctax1;

  consumption = p.w/hp.H;
  p.consumption=consumption;
#  E= consumption/(gdp1 - G1 - financialcosts1 - investment1)
  E=( hp.H*(labor_d1 + p.w^(-1.0)*((1-tau.i)*p.r*bonds1+ netdistributions1 +liquidations1) ) )^-1.0;
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
end



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
