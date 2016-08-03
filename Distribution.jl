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
function stationarydist(E::Real, pr::FirmProblem, eq::Equilibrium, tau::Taxes, pa::Param;  tol=eps(), verbose=true)
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
          omegaprime, i_omegaprime = predict_state(i_zprime, i_omega, i_z, pr, eq, tau, pa);
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
        omega0= omegaprimefun(pa.k0, 0.0, i_z, eq, tau, pa)
        omega0_ind = closestindex(omega0, pa.omega.step);
        ie[k]= find_dist_ind(omega0_ind, i_z, pa.Nomega);
        ve[k]= E*pa.invariant_distr[i_z];
        k+=1;
  end

  entrants= sparsevec(ie,ve,pa.Nz*pa.Nomega);
  I=speye(Float64,pa.Nz*pa.Nomega); #Identity matrix

  distrguess = ones(pa.Nomega* pa.Nz);
  distrguess=  distrguess/sum(distrguess);

  distr_vectorized=deepcopy(distrguess);
  dist=Inf;
  it=1; maxit = 10.0^4.0;
  while dist >tol && it<maxit
    verbose && println("it = ",it," dist = ", dist)
    #println(it)
    distrprime= transmat'*distr_vectorized +entrants;
    dist= norm(distrprime -distr_vectorized);

    distr_vectorized=deepcopy(distrprime);
      distr=reshape(distr_vectorized, (pa.Nomega,pa.Nz))
      probexit = sum(distr.*pr.exitprobability)
      if probexit ==0
        error("No firm wants to exit")
      end
    it+=1;
  end
  if it==maxit
    error("Distribution did not converge after maximum number of iterations. Max(exitrule) = ", maximum(pr.exitrule), " Min(exitrule) = ", minimum(pr.exitrule), "Determinant of transition matrix = ", det(transmat'))
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
  it=1; maxit = 10.0^3.0;
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
    omega0= omegaprimefun(pa.k0, 0.0, i_z, eq, tau, pa)
    omega0_ind = closestindex(omega0, pa.omega.step);
    distrprime[omega0_ind,i_z]+=E*pa.invariant_distr[i_z];
    #Non-exiting incumbents
    for i_omega in 1:pa.Nomega
      kprime=pr.kpolicy[i_omega,i_z];
      qprime=pr.qpolicy[i_omega,i_z];
      for i_zprime in 1:pa.Nz
        if !pr.exitrule[i_omega,i_z,i_zprime]
          omegaprime, i_omegaprime = predict_state(i_zprime, i_omega, i_z, pr, eq, tau, pa);
          distrprime[i_omegaprime,i_zprime] += pa.ztrans[i_zprime,i_z]*distr[i_omega,i_z];
        end

      end
    end
  end
distrprime
end

function distributionStupid(E::Real,pr::FirmProblem, eq::Equilibrium, tau::Taxes, pa::Param;  tol=eps())

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
          omegaprime, i_omegaprime = predict_state(i_zprime, i_omega, i_z, pr, eq, tau, pa);
          nextDistr[i_omegaprime,i_zprime] += pa.ztrans[i_zprime,i_z]*initialDistr[i_omega,i_z]
        end
      end # end i_z'
      omega0= omegaprimefun(pa.k0, 0.0, i_z, eq, tau, pa)
      omega0_ind = closestindex(omega0, pa.omega.step);
      if i_omega == omega0_ind;
        nextDistr[i_omega,i_z] += E*pa.invariant_distr[i_z]
      end
    end # end i_omega, i_z

    diffValue    = maxabs(nextDistr - initialDistr)
    initialDistr = deepcopy(nextDistr)
    nextDistr    = zeros(Float64,pa.Nomega,pa.Nz)

  end

  initialDistr

end

function computeMomentsCutoff(E::Real,pr::FirmProblem, eq::Equilibrium, tau::Taxes, pa::Param; cutoffCapital::Float64=0.0, toPrint::Bool=true)

  indexCutoff = 1
  for i in 1:pa.Nomega
    if (pa.omega.grid[i]-pa.omega.lb)/(pa.omega.ub-pa.omega.lb) >= cutoffCapital
      indexCutoff = i
      break
    end
  end

  # Now, since taking averages, we need a factor to weight the distribution
  # so that it adds up to one, taking into account the ones we are cutting off
  # dividing by this makes the mass add up to one, but maintains relative proportions
  # of firms in each state
  massCorrection  = 0.0
  mass2           =0.0
  masslosses      =0.0
  for i_omega in indexCutoff:pa.Nomega, i_z in 1:pa.Nz
    firstKPrime = pr.kpolicy[i_omega,i_z]; # since we drop obs with k= 0, adjust for that as well.
    massCorrection += firstKPrime > 0.0 ? eq.distr[i_omega,i_z]:0.0;
    mass2 += eq.distr[i_omega,i_z];
    for i_zprime in 1:pa.Nz
      masslosses += profits(pa.zgrid[i_zprime],firstKPrime,eq,pa) < 0.0 && firstKPrime > 0.0 ? eq.distr[i_omega,i_z]*pa.ztrans[i_zprime,i_z]:0.0;
    end
  end

  # Compute moments

  mean_inv_rate     = 0.0
  var_inv_rate      = 0.0
  mean_leverage     = 0.0
  var_leverage      = 0.0
  mean_dividends2k  = 0.0
  var_dividends2k   = 0.0
  mean_profits2k    = 0.0
  var_profits2k     = 0.0
  mean_eqis         = 0.0
  freq_equis        = 0.0
  mean_tobinsq      = 0.0
  autocov_profits2k = 0.0
  losses2k          = 0.0
  capital= 0.0

  # First the means
  for i_omega in indexCutoff:pa.Nomega, i_z in 1:pa.Nz

    firstKPrime = pr.kpolicy[i_omega,i_z]
    firstQPrime = pr.qpolicy[i_omega,i_z]

    capital+= firstKPrime*eq.distr[i_omega,i_z]/mass2;

    # assuming leverage is debt/capital, make sure previous capital was not zero
    mean_leverage += firstKPrime > 0.0 ? (firstQPrime/firstKPrime)*eq.distr[i_omega,i_z]/massCorrection : 0.0
    # if dividends below zero actually they are distributions
    mean_eqis      += !pr.positivedistributions[i_omega,i_z] ? -pr.grossequityis[i_omega,i_z]*eq.distr[i_omega,i_z]/mass2 : 0.0
    freq_equis     += !pr.positivedistributions[i_omega,i_z] ? eq.distr[i_omega,i_z]/mass2  : 0.0

    for i_zprime in 1:pa.Nz
      #INCUMBENTS
#      if !pr.exitrule[i_omega, i_z,i_zprime]
        omegaprime, i_omegaprime = predict_state(i_zprime, i_omega, i_z, pr, eq, tau, pa);
        secondKPrime = pr.kpolicy[i_omegaprime,i_zprime]

        # Dividend to capital ratio

        mean_dividends2k += pr.positivedistributions[i_omegaprime,i_zprime] && firstKPrime>0.0 ?
          (pr.grossdividends[i_omegaprime,i_zprime]/firstKPrime)*pa.ztrans[i_zprime,i_z]*eq.distr[i_omega,i_z]/massCorrection : 0.0

        # Investment, make sure previous capital was not zero
        mean_inv_rate += firstKPrime > 0.0 ?
          ((secondKPrime-(1-pa.delta)*firstKPrime)/firstKPrime)*pa.ztrans[i_zprime,i_z]*eq.distr[i_omega,i_z]/massCorrection :
          0.0

        # Profits
        mean_profits2k +=firstKPrime > 0.0 ?
         (profits(pa.zgrid[i_zprime],firstKPrime,eq,pa)/firstKPrime)*pa.ztrans[i_zprime,i_z]*eq.distr[i_omega,i_z]/massCorrection :0.0

        #Losses
        losses2k +=profits(pa.zgrid[i_zprime],firstKPrime,eq,pa) < 0.0 && firstKPrime > 0.0 ?
         (profits(pa.zgrid[i_zprime],firstKPrime,eq,pa)/firstKPrime)*pa.ztrans[i_zprime,i_z]*eq.distr[i_omega,i_z]/masslosses:0.0

        # Tobin's Q
        mean_tobinsq += firstKPrime > 0.0 ?
          ((pr.firmvaluegrid[i_omegaprime,i_zprime]+ firstQPrime )/firstKPrime)*pa.ztrans[i_zprime,i_z]*eq.distr[i_omega,i_z]/massCorrection:0.0;
#      end

    end
  end

  # Now the variances
  for i_omega in indexCutoff:pa.Nomega, i_z in 1:pa.Nz
    firstKPrime  = pr.kpolicy[i_omega,i_z]
    firstQPrime  = pr.qpolicy[i_omega,i_z]
    var_leverage += firstKPrime > 0.0 ?
      ((firstQPrime/firstKPrime-mean_leverage)^2)*eq.distr[i_omega,i_z]/massCorrection : 0.0

    for i_zprime in 1:pa.Nz
      omegaprime, i_omegaprime = predict_state(i_zprime, i_omega, i_z, pr, eq, tau, pa);
      secondKPrime = pr.kpolicy[i_omegaprime,i_zprime]
      secondQPrime = pr.qpolicy[i_omegaprime,i_zprime]

      # Investment, make sure previous capital was not zero
      var_inv_rate += firstKPrime > 0.0 ?
        ((((secondKPrime-(1-pa.delta)*firstKPrime)/firstKPrime)-mean_inv_rate)^2)*pa.ztrans[i_zprime,i_z]*eq.distr[i_omega,i_z]/massCorrection :
        0.0

      # Dividends
      var_dividends2k += pr.positivedistributions[i_omegaprime,i_zprime] && firstKPrime>0.0 ?
          ((pr.grossdividends[i_omegaprime,i_zprime]/firstKPrime  - mean_dividends2k)^2.0)*pa.ztrans[i_zprime,i_z]*eq.distr[i_omega,i_z]/massCorrection : 0.0

      # Profits
      currentProfits = profits(pa.zgrid[i_zprime],firstKPrime,eq,pa)/firstKPrime;
      var_profits2k += firstKPrime > 0.0 ?
        ((currentProfits-mean_profits2k)^2.0)*pa.ztrans[i_zprime,i_z]*eq.distr[i_omega,i_z]/massCorrection:0.0;

      for i_zprime2 in 1:pa.Nz
        currentProfits2 =profits(pa.zgrid[i_zprime2],secondKPrime,eq,pa)/secondKPrime;
        autocov_profits2k += (firstKPrime >0.0 && secondKPrime > 0.0) ?
          ((currentProfits2-mean_profits2k)*(currentProfits-mean_profits2k))*pa.ztrans[i_zprime2,i_zprime]*pa.ztrans[i_zprime,i_z]*eq.distr[i_omega,i_z]/massCorrection:0.0;
      end
    end
  end

  autocov_profits2k = autocov_profits2k/var_profits2k;
  turnover=eq.E/sum(eq.distr);
  labor = eq.a.laborsupply;

  resultingMoments = Moments(mean_inv_rate,sqrt(var_inv_rate),mean_leverage,sqrt(var_leverage),mean_dividends2k,sqrt(var_dividends2k),
    mean_profits2k,sqrt(var_profits2k),mean_eqis/capital,freq_equis,mean_tobinsq,autocov_profits2k)

  if toPrint
#    using DataFrames
    namesMoments = ["Mean Investment","Mean Leverage","Mean Profits","SD Profits",
    "Mean Equity Issuance","Frequency of Equity Issuance","Autocovariance Profits",
    "Turnover","Time At Work","SD Leverage","Mean Dividends","SD Dividends","SD Investment",
    "Means Tobins Q"]

    valueMoments = [mean_inv_rate,mean_leverage,mean_profits2k,sqrt(var_profits2k),
    mean_eqis/capital,freq_equis,autocov_profits2k,turnover,labor,sqrt(var_leverage),
    sqrt(var_inv_rate),mean_dividends2k,sqrt(var_dividends2k),mean_tobinsq]
    println(DataFrame(names=namesMoments,aggVals=valueMoments))
  end

  resultingMoments
end
