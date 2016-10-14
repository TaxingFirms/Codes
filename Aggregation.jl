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
  liquidationtax = 0.0;

  for i_zprime in 1:pa.Nz
    #### First we consider entrants ###
      z0       = pa.zgrid[i_zprime];
      l0       = (z0*pa.alphal*pa.k0^pa.alphak / eq.w)^(1/(1-pa.alphal));

      debt        += 0.0*pa.invariant_distr[i_zprime];
      labor       += l0*pa.invariant_distr[i_zprime];
      omega0= omegaprimefun(pa.k0, 0.0, i_zprime, eq, tau, pa)
      i_omega0 = closestindex(omega0, pa.omega.step);
      netdistributions+=pr.distributions[i_omega0,i_zprime]*pa.invariant_distr[i_zprime];


    #### Next we consider incumbents ###
    for i_omega in 1:pa.Nomega
      for i_z in 1:pa.Nz
        kprime=pr.kpolicy[i_omega,i_z];
        qprime=pr.qpolicy[i_omega,i_z];

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
          liquidations += (1-tau.exit)*(pa.kappa*(1-pa.delta)*kprime - (1+eq.r)*qprime)*distr1[i_omega,i_z]*pa.ztrans[i_zprime,i_z];
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



function mass_of_entrants!( pr::FirmProblem, eq::Equilibrium, tau::Taxes, pa::Param, distr_routine::Function; verbose=true)
  #Computes the mass of entrants such that the labor market clears,
  verbose?distr1=distr_routine(1.0,pr, eq, tau,pa):distr1=distr_routine(1.0,pr, eq, tau,pa;verbose=false);
  bonds1, labor_d1, netdistributions1, liquidations1 = unit_entry(distr1, pr, eq, tau, pa);

  E=( pa.H*(labor_d1 + eq.w^(-1.0)*((1-tau.i)*eq.r*bonds1+ netdistributions1 +liquidations1) ) )^-1.0;
  eq.distr = E*distr1;
  eq.E = E;
end


function mass_of_entrantsGHH!( pr::FirmProblem, eq::Equilibrium, tau::Taxes, pa::Param, distr_routine::Function; verbose=true)
  #Computes the mass of entrants such that the labor market clears,
  verbose ? distr1=distr_routine(1.0,pr, eq, tau,pa): distr1=distr_routine(1.0,pr, eq, tau,pa;verbose=false);
  bonds1, labor_d1, netdistributions1, liquidations1 = unit_entry(distr1, pr, eq, tau, pa);
  labor_s= ((1-tau.l)*eq.w/pa.H)^pa.psi

  E= labor_s/labor_d1;
  eq.distr = E*distr1;
  eq.E = E;
end


################################################################################
################################################################################



function aggregates!(pr::FirmProblem, eq::Equilibrium, tau::Taxes, pa::Param; compute_moments::Bool=false, consistencychecks::Bool=true)
# Computes aggregates and moments once the model is completely solved.
  distr=eq.distr;
  taudtilde = 1-(1-tau.d)/(1-tau.g);
  kprimefun = map(x->CoordInterpGrid(pa.omega.grid,pr.kpolicy[:,x],BCnearest, InterpLinear),1:pa.Nz);
  qprimefun = map(x->CoordInterpGrid(pa.omega.grid,pr.qpolicy[:,x],BCnearest, InterpLinear),1:pa.Nz);

  ########## Initialize accumulators ###########
  # Aggregates
  capital=0.0;
  debt=0.0;
  labor=0.0;
  gdp=0.0;
  corptax=0.0;
  liquidations=0.0;
  liquidationcosts=0.0;
  liquidationtax=0.0;
  netdistributions=0.0;
  divs=0.0;


  ###### COMPUTE VALUES THAT DEPEND ON CURRENT K
  for i_zprime in 1:pa.Nz
    #### First we consider entrants ###
    omega0= omegaprimefun(pa.k0, 0.0, i_zprime, eq, tau, pa)
    omega0_ind = closestindex(omega0, pa.omega.step);

    netdistributions += pr.distributions[omega0_ind,i_zprime]*eq.E*pa.invariant_distr[i_zprime];
    divs += max(omega0 -pa.f - pr.kpolicy[omega0_ind,i_zprime] + pr.qpolicy[omega0_ind,i_zprime] ,0.0)*eq.E*pa.invariant_distr[i_zprime];
    k0=pa.k0;
    q0=0.0;
    z0 = pa.zgrid[i_zprime];
    l0    = (z0*pa.alphal*k0^pa.alphak / eq.w)^(1/(1-pa.alphal));
    labor    += l0*eq.E*pa.invariant_distr[i_zprime];
    gdp      += (z0*k0^pa.alphak*l0^pa.alphal - pa.f)*eq.E*pa.invariant_distr[i_zprime];
    corptax  += tau.c*max(z0*k0^pa.alphak*l0^pa.alphal - eq.w*l0 - pa.allowance*pa.delta*k0 - eq.r*q0 - pa.f,0.0)*eq.E*pa.invariant_distr[i_zprime];

    #### Next we consider incumbents ###
    for i_z in 1:pa.Nz

      for i_omega in 1:pa.Nomega
        kprime=pr.kpolicy[i_omega,i_z];
        qprime=pr.qpolicy[i_omega,i_z];
        omega = pa.omega.grid[i_omega];


        #Non-exiting incumbents
        if !pr.exitrule[i_omega,i_z,i_zprime]
          zprime     = pa.zgrid[i_zprime];
          lprime     = (zprime*pa.alphal*kprime^pa.alphak / eq.w)^(1/(1-pa.alphal));

          capital   += kprime*distr[i_omega,i_z]*pa.ztrans[i_zprime,i_z];
          debt      += qprime*distr[i_omega,i_z]*pa.ztrans[i_zprime,i_z];
          labor     += lprime*distr[i_omega,i_z]*pa.ztrans[i_zprime,i_z];

          gdp       += (zprime*kprime^pa.alphak*lprime^pa.alphal - pa.f)*distr[i_omega,i_z]*pa.ztrans[i_zprime,i_z];
          corptax   += tau.c*max(zprime*kprime^pa.alphak*lprime^pa.alphal -eq.w*lprime - pa.allowance*pa.delta*kprime - eq.r*qprime - pa.f,0.0)*distr[i_omega,i_z]*pa.ztrans[i_zprime,i_z];

          omegaprime, i_omegaprime = predict_state(i_zprime, i_omega, i_z, pr, eq, tau, pa);

          omegaprimecheck= zprime*kprime^pa.alphak*lprime^pa.alphal - eq.w*lprime+ (1-pa.delta)*kprime - (1+eq.r)*qprime - tau.c*max(zprime*kprime^pa.alphak*lprime^pa.alphal -eq.w*lprime - pa.allowance*pa.delta*kprime - eq.r*qprime - pa.f,0.0)

          divs += max(omegaprimecheck -pa.f - kprimefun[i_zprime][omegaprimecheck] + qprimefun[i_zprime][omegaprimecheck] , 0.0)*distr[i_omega,i_z]*pa.ztrans[i_zprime,i_z];
          divscheck = omegaprime - pa.f - pr.kpolicy[i_omegaprime,i_zprime] + pr.qpolicy[i_omegaprime,i_zprime] ;
          netdistributions+=pr.distributions[i_omegaprime,i_zprime]*distr[i_omega,i_z]*pa.ztrans[i_zprime,i_z];
          corptaxbase =  zprime*kprime^pa.alphak*lprime^pa.alphal -eq.w*lprime - pa.allowance*pa.delta*kprime - eq.r*qprime - pa.f;

        #Exiting incumbents
        else
          liquidationcosts += (1-pa.kappa)*(1-pa.delta)*kprime*distr[i_omega,i_z]*pa.ztrans[i_zprime,i_z]
          liquidations += (1-tau.exit)*(pa.kappa*(1-pa.delta)*kprime - (1+eq.r)*qprime)*distr[i_omega,i_z]*pa.ztrans[i_zprime,i_z];
          liquidationtax+= tau.exit*(pa.kappa*(1-pa.delta)*kprime - (1+eq.r)*qprime)*distr[i_omega,i_z]*pa.ztrans[i_zprime,i_z]
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

  labor_s = ((1-tau.l)*eq.w/pa.H)^pa.psi;
  deductionfactor=0.35;
  deduction_l = deductionfactor*tau.l*eq.w*labor_s;
  deduction_i = deductionfactor*tau.i*eq.r*debt;
  consumption = (1-tau.l)*eq.w*labor_s + (1-tau.i)*eq.r*debt + netdistributions + liquidations - eq.E*(pa.k0+pa.e) + deduction_l+ deduction_i;

  divtax= tau.d*grossdividends;
  inctax = tau.i*eq.r*debt - deduction_i;
  labtax = tau.l*eq.w*labor - deduction_l;

  G = divtax + corptax + inctax +labtax + liquidationtax;
  ######################################



  ######## Consistency Checks #########
  netdistributionscheck = sum(distr.*pr.distributions);
  debtcheck = sum(distr.*pr.qpolicy)
  capitalcheck= sum(distr.*pr.kpolicy)
  if consistencychecks
    if abs(capital - capitalcheck) >10.0^-4.0 ||
        abs(debt - debtcheck) >10.0^-4.0 ||
        abs(netdistributions - netdistributionscheck) >10.0^-4.0
        # || abs(divs - grossdividends) >10.0^-3.0
      error("Consistency problem")
    end

    if abs(labor_s - labor) > 10^-10.0
      error("labor market didn't clear")
    end

    if abs((gdp - G - financialcosts - investment -liquidationcosts  ) - consumption)>10.0^-2.0
          println("goods market didn't clear by: ", (gdp - G - financialcosts - investment -liquidationcosts  ) - consumption)
    end
  end


  ######## Save values #########
  logarg = consumption - (pa.H/(1+1/pa.psi))* labor_s^(1+1/pa.psi);
  logarg<0 && error("Negative Argument for log utility")

  welfare = pa.sigma==1.0 ?
    log( logarg ) /(1-pa.beta):
    (1/(1-pa.sigma)*( logarg )^(1-pa.sigma) ) /(1-pa.beta);
  collections = Taxes(divtax,corptax,inctax,0.0,labtax,liquidationtax);
  eq.a=Aggregates(netdistributions, (1-tau.i)*eq.r*debt, consumption, gdp, labor, welfare, collections, debt, capital, investment, grossdividends, financialcosts, G);


end


function aggregatesold!(pr::FirmProblem, eq::Equilibrium, tau::Taxes, pa::Param; compute_moments::Bool=false)
# Computes aggregates and moments once the model is completely solved.
  distr=eq.distr;
  taudtilde = 1-(1-tau.d)/(1-tau.g);
  kprimefun = map(x->CoordInterpGrid(pa.omega.grid,pr.kpolicy[:,x],BCnearest, InterpLinear),1:pa.Nz);
  qprimefun = map(x->CoordInterpGrid(pa.omega.grid,pr.qpolicy[:,x],BCnearest, InterpLinear),1:pa.Nz);

  ########## Initialize accumulators ###########
  # Aggregates
  capital=0.0;
  debt=0.0;
  labor=0.0;
  gdp=0.0;
  corptax=0.0;
  liquidations=0.0;
  liquidationcosts=0.0;
  liquidationtax=0.0;
  netdistributions=0.0;
  divs=0.0;


  #Moments
  mass_incumbents=0.0;
  if compute_moments
    mean_inv_rate=0.0;
    ss_inv_rate=0.0; #sum of squares
    mean_leverage=0.0;
    ss_leverage=0.0;
    mean_dividends2k=0.0;
    ss_dividends2k=0.0;
    mean_profits2k=0.0;
    ss_profits2k=0.0;
    mean_eqis=0.0;
    freq_eqis=0.0;
    mean_tobinsq=0.0;
    scov_profits2k=0.0;
    cov_nw = 0.0;

    # Auxiliary variables
    mean_inv_rate_shifted=0.0;
    mean_leverage_shifted=0.0;
    mean_dividends2k_shifted=0.0;
    mean_profits2k_shifted=0.0;
    mean_omega_shifted=0.0;
    mean_omegaprime_shifted=0.0;
    mean_profits2ksecond_shifted=0.0;


    ##Define constants to avoid catastrophic cancellation
    # (Understanding this isn't important)
    # Compute the value at midpoints
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
  for i_zprime in 1:pa.Nz
    #### First we consider entrants ###
    omega0= omegaprimefun(pa.k0, 0.0, i_zprime, eq, tau, pa)
    omega0_ind = closestindex(omega0, pa.omega.step);

    netdistributions += pr.distributions[omega0_ind,i_zprime]*eq.E*pa.invariant_distr[i_zprime];
    divs += max(omega0 -pa.f - pr.kpolicy[omega0_ind,i_zprime] + pr.qpolicy[omega0_ind,i_zprime] ,0.0)*eq.E*pa.invariant_distr[i_zprime];
    k0=pa.k0;
    q0=0.0;
    z0 = pa.zgrid[i_zprime];
    l0    = (z0*pa.alphal*k0^pa.alphak / eq.w)^(1/(1-pa.alphal));
    labor    += l0*eq.E*pa.invariant_distr[i_zprime];
    gdp      += (z0*k0^pa.alphak*l0^pa.alphal - pa.f)*eq.E*pa.invariant_distr[i_zprime];
    corptax  += tau.c*max(z0*k0^pa.alphak*l0^pa.alphal - eq.w*l0 - pa.allowance*pa.delta*k0 - eq.r*q0 - pa.f,0.0)*eq.E*pa.invariant_distr[i_zprime];

    if compute_moments
      freq_eqis+=(1-pr.positivedistributions[omega0_ind,i_zprime])*eq.E*pa.invariant_distr[i_zprime];
      mean_eqis+= - pr.grossequityis[omega0_ind,i_zprime]*eq.E*pa.invariant_distr[i_zprime];

      #investment rate, leverage, etc are undefinied for entrants
    end

    #### Next we consider incumbents ###
    for i_z in 1:pa.Nz

      for i_omega in 1:pa.Nomega
        kprime=pr.kpolicy[i_omega,i_z];
        qprime=pr.qpolicy[i_omega,i_z];
        omega = pa.omega.grid[i_omega];


        #Non-exiting incumbents
        if !pr.exitrule[i_omega,i_z,i_zprime]
          zprime     = pa.zgrid[i_zprime];
          lprime     = (zprime*pa.alphal*kprime^pa.alphak / eq.w)^(1/(1-pa.alphal));

          capital   += kprime*distr[i_omega,i_z]*pa.ztrans[i_zprime,i_z];
          debt      += qprime*distr[i_omega,i_z]*pa.ztrans[i_zprime,i_z];
          labor     += lprime*distr[i_omega,i_z]*pa.ztrans[i_zprime,i_z];

          gdp       += (zprime*kprime^pa.alphak*lprime^pa.alphal - pa.f)*distr[i_omega,i_z]*pa.ztrans[i_zprime,i_z];
          corptax   += tau.c*max(zprime*kprime^pa.alphak*lprime^pa.alphal -eq.w*lprime - pa.allowance*pa.delta*kprime - eq.r*qprime - pa.f,0.0)*distr[i_omega,i_z]*pa.ztrans[i_zprime,i_z];

          omegaprime, i_omegaprime = predict_state(i_zprime, i_omega, i_z, pr, eq, tau, pa);

          omegaprimecheck= zprime*kprime^pa.alphak*lprime^pa.alphal - eq.w*lprime+ (1-pa.delta)*kprime - (1+eq.r)*qprime - tau.c*max(zprime*kprime^pa.alphak*lprime^pa.alphal -eq.w*lprime - pa.allowance*pa.delta*kprime - eq.r*qprime - pa.f,0.0)

          divs += max(omegaprimecheck -pa.f - kprimefun[i_zprime][omegaprimecheck] + qprimefun[i_zprime][omegaprimecheck] , 0.0)*distr[i_omega,i_z]*pa.ztrans[i_zprime,i_z];
          divscheck = omegaprime - pa.f - pr.kpolicy[i_omegaprime,i_zprime] + pr.qpolicy[i_omegaprime,i_zprime] ;
          netdistributions+=pr.distributions[i_omegaprime,i_zprime]*distr[i_omega,i_z]*pa.ztrans[i_zprime,i_z];
          corptaxbase =  zprime*kprime^pa.alphak*lprime^pa.alphal -eq.w*lprime - pa.allowance*pa.delta*kprime - eq.r*qprime - pa.f;


          if compute_moments
            freq_eqis+=(1-pr.positivedistributions[i_omegaprime,i_zprime])*distr[i_omega,i_z]*pa.ztrans[i_zprime,i_z];
            mean_eqis+= - pr.grossequityis[i_omegaprime,i_zprime]*distr[i_omega,i_z]*pa.ztrans[i_zprime,i_z];
            ksecond = pr.kpolicy[i_omegaprime,i_zprime];

            if kprime>0
              inv_rate = ((ksecond - (1-pa.delta)*kprime)/kprime);
              mean_inv_rate += inv_rate*distr[i_omega,i_z]*pa.ztrans[i_zprime,i_z];
              mean_inv_rate_shifted += (inv_rate- Kinv)*distr[i_omega,i_z]*pa.ztrans[i_zprime,i_z];
              ss_inv_rate += (inv_rate-Kinv)^2*distr[i_omega,i_z]*pa.ztrans[i_zprime,i_z];

              leverage = qprime/kprime;
              mean_leverage += leverage*distr[i_omega,i_z]*pa.ztrans[i_zprime,i_z];
              mean_leverage_shifted += (leverage - Klev)*distr[i_omega,i_z]*pa.ztrans[i_zprime,i_z];
              ss_leverage += (leverage - Klev)^2*distr[i_omega,i_z]*pa.ztrans[i_zprime,i_z];

              div2k= pr.grossdividends[i_omegaprime,i_zprime]/kprime; #Before tax dividends
              mean_dividends2k += div2k*distr[i_omega,i_z]*pa.ztrans[i_zprime,i_z];
              mean_dividends2k_shifted += (div2k - Kdiv)*distr[i_omega,i_z]*pa.ztrans[i_zprime,i_z];
              ss_dividends2k += (div2k - Kdiv)^2*distr[i_omega,i_z]*pa.ztrans[i_zprime,i_z];

              prof2k = (zprime*kprime^pa.alphak*lprime^pa.alphal -eq.w*lprime -pa.f)/kprime;
              mean_profits2k += prof2k*distr[i_omega,i_z]*pa.ztrans[i_zprime,i_z];
              mean_profits2k_shifted += (prof2k - Kprof)*distr[i_omega,i_z]*pa.ztrans[i_zprime,i_z];
              ss_profits2k +=  (prof2k - Kprof)^2*distr[i_omega,i_z]*pa.ztrans[i_zprime,i_z];

              tobinsq = pr.firmvaluegrid[i_omegaprime,i_zprime]/kprime;
              mean_tobinsq += tobinsq*distr[i_omega,i_z]*pa.ztrans[i_zprime,i_z];

              cov_nw += (omegaprime - Komega)*(omega-Komega)*distr[i_omega,i_z]*pa.ztrans[i_zprime,i_z];
              mean_omega_shifted += (omega - Komega);
              mean_omegaprime_shifted += (omegaprime - Komega);

              for i_zsecond in 1:pa.Nz
                if !pr.exitrule[i_omegaprime,i_zprime,i_zsecond] && ksecond >0
                  zsecond = pa.zgrid[i_zsecond];
                  lsecond = (zsecond*pa.alphal*ksecond^pa.alphak / eq.w)^(1/(1-pa.alphal));
                  prof2ksecond = (zsecond*ksecond^pa.alphak*lsecond^pa.alphal - eq.w*lsecond -pa.f)/ksecond;

                  mean_profits2ksecond_shifted += (prof2ksecond - Kprof)*distr[i_omega,i_z]*pa.ztrans[i_zprime,i_z]*pa.ztrans[i_zsecond,i_zprime];
                  scov_profits2k += (prof2k - Kprof)*(prof2ksecond - Kprof)*distr[i_omega,i_z]*pa.ztrans[i_zprime,i_z]*pa.ztrans[i_zsecond,i_zprime];
                end
              end

            end
          end
          mass_incumbents += distr[i_omega,i_z]*pa.ztrans[i_zprime,i_z];

        #Exiting incumbents
        else
          liquidationcosts += (1-pa.kappa)*(1-pa.delta)*kprime*distr[i_omega,i_z]*pa.ztrans[i_zprime,i_z]
          liquidations += (1-taudtilde)*(pa.kappa*(1-pa.delta)*kprime - (1+eq.r)*qprime)*distr[i_omega,i_z]*pa.ztrans[i_zprime,i_z];
          liquidationtax+= taudtilde*(pa.kappa*(1-pa.delta)*kprime - (1+eq.r)*qprime)*distr[i_omega,i_z]*pa.ztrans[i_zprime,i_z]
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
  labtax = tau.l*eq.w*labor;

  G = divtax + corptax + inctax +labtax + liquidationtax;
  ######################################



  ######## Consistency Checks #########
  netdistributionscheck = sum(distr.*pr.distributions);
  debtcheck = sum(distr.*pr.qpolicy)
  capitalcheck= sum(distr.*pr.kpolicy)
  labor_s = ((1-tau.l)*eq.w/pa.H)^pa.psi;
  consumption = (1-tau.l)*eq.w*labor_s +  (1-tau.i)*eq.r*debt+ netdistributions +liquidations - eq.E*(pa.k0+pa.e);

  if abs(capital - capitalcheck) >10.0^-4.0 ||
      abs(debt - debtcheck) >10.0^-4.0 ||
      abs(netdistributions - netdistributionscheck) >10.0^-4.0
      # || abs(divs - grossdividends) >10.0^-3.0
    error("Consistency problem")
  end

  if abs(labor_s - labor) > 10^-10.0
    error("labor market didn't clear")
  end

  if abs((gdp - G - financialcosts - investment -liquidationcosts  ) - consumption)>10.0^-2.0
        println("goods market didn't clear by: ", (gdp - G - financialcosts - investment -liquidationcosts  ) - consumption)
  end

  ##############################


  ########## Moments ##########
  if compute_moments
    mean_inv_rate = mean_inv_rate / mass_incumbents;
    var_inv_rate=  (ss_inv_rate - mean_inv_rate_shifted^2/mass_incumbents )/ mass_incumbents;
    sd_inv_rate = sqrt(var_inv_rate);

    mean_leverage=mean_leverage/ mass_incumbents;
    var_leverage=(ss_leverage - mean_leverage_shifted^2/mass_incumbents)/ mass_incumbents;
    sd_leverage= sqrt(var_leverage);

    mean_dividends2k= mean_dividends2k/ mass_incumbents;
    var_dividends2k=(ss_dividends2k - mean_dividends2k_shifted^2/mass_incumbents)/ mass_incumbents;
    sd_dividends2k = sqrt(var_dividends2k);

    mean_profits2k=mean_profits2k/ mass_incumbents;
    var_profits2k=(ss_profits2k - mean_profits2k_shifted^2/mass_incumbents )/ mass_incumbents;
    sd_profits2k= sqrt(var_profits2k);

    mean_eqis2k= mean_eqis/ capital;
    freq_eqis=freq_eqis/ mass_incumbents;
    mean_tobinsq=mean_tobinsq/ mass_incumbents;

    cov_nw=(cov_nw - mean_omega_shifted*mean_omegaprime_shifted / mass_incumbents)/ mass_incumbents;
    cov_profits2k= (scov_profits2k - mean_profits2k_shifted*mean_profits2ksecond_shifted / mass_incumbents)/mass_incumbents;
    correl_profits2k = cov_profits2k/var_profits2k;

    turnover=eq.E/sum(eq.distr);
    labor = eq.a.labor;

    eq.m = Moments(mean_inv_rate, sd_inv_rate, mean_leverage, sd_leverage, mean_dividends2k, sd_dividends2k, mean_profits2k, sd_profits2k, mean_eqis2k, freq_eqis, mean_tobinsq, correl_profits2k,turnover,labor)
  end
  ##############################




  ######## Save values #########
  welfare = pa.sigma==1 ?
    (log(consumption - (pa.H/(1+pa.psi))* labor_s^(1+pa.psi) ) ) /(1-pa.beta):
    (1/(1-pa.sigma)*(consumption - (pa.H/1+pa.psi)* labor_s^(1+pa.psi) )^(1-pa.sigma) ) /(1-pa.beta);
  collections = Taxes(divtax,corptax,inctax,0.0,labtax);

  eq.a=Aggregates(netdistributions, (1-tau.i)*eq.r*debt, consumption, gdp, labor, welfare, collections, debt, capital, investment, grossdividends, financialcosts, G);


end
