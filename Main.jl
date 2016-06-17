## Performance tips
# @time
#first index -> inner loop
# PREALLOCATE OUTPUTS
# Avoid unnecessary arrays
# Use @inbounds
# define type asap


# Define types. This will allow to be more organized when passing parameters to functions.
immutable GridObject
  ub:: Real  #upper bound
  lb:: Real  #lower bound
  step:: Real #Distance between grid points
  N:: Int    #Number of gridpoints
  grid:: AbstractArray #Grid
end

immutable FirmParam
  alphak::Real #Capital share of output
  alphal::Real #Labor share of output
  f::Real
  lambda0::Real
  lambda1::Real
  delta::Real # Capital depreciation
  theta::Real # Collateral
  kappa::Real #Liquidation cost
  e::Real # Entry cost
  collateral_factor::Real #theta*(1-delta)
  leverageratio::Real # 1/(1-theta*(1-delta)), leverage at colateral and no divindend
  zgrid::Array{Float64,1}
  ztrans::Array{Float64,2}
  invariant_distr::Array #Invariant distribution
  Nz::Int64
  Nk::Int64
  Nq::Int64
  Nomega::Int64
end

immutable HouseholdParam
  beta::Real #Discount rate
  sigma::Real  #Risk aversion/ ies
  H::Real  # Labor supply parameter
  psi::Real #Labor supply
end

type Taxes
  d::Real
  c::Real
  i::Real
  g::Real
end


type FirmProblem
  #input
  betatilde::Real
  taudtilde::Real
  omega::GridObject
  firmvalueguess::Matrix
  kprime::GridObject
  qprime::GridObject
  discounted_interest::Real # betatilde*(1+(1-tau.c)*r)
  #output
  firmvaluegrid::Matrix
  kpolicygrid:: Matrix
  qpolicygrid::Matrix
  Nk::Int64
  Nq::Int64
  Nz::Int64
  Nomega::Int64
  InterpolationGrid::Array{CoordInterpGrid,1}
end

type ResultsFP
  firmvalue :: Array{Float64,2}
  interpolation::Array{CoordInterpGrid,1}
  kprime :: Array{Float64,2}
  qprime :: Array{Float64,2}
  distributions :: Array{Float64,2}
  financialcosts :: Array{Float64,2}
  grossdividends :: Array{Float64,2}
  grossequityis :: Array{Float64,2}
  exitprobability :: Array{Float64,2}
  exitrule:: Array{Bool,3}
  positivedistributions :: Array{Bool,2}
end

type Moments
  mean_inv_rate::Float64
  var_inv_rate::Float64
  mean_leverage::Float64
  var_leverage::Float64

  mean_dividends2k::Float64
  var_dividends2k::Float64
  mean_profits2k::Float64

  var_profits2k::Float64
  mean_eqis2k::Float64
  freq_equis2k::Float64
  mean_tobinsq::Float64

  cov_nw::Float64
end

type Aggregates
  netdistributions::Float64
  agginterests::Float64
  consumption::Float64
  output::Float64
  laborsupply::Float64
  welfare::Float64
  collections::Taxes
  bonds::Float64
  capital::Float64
  investment::Float64
  grossdividends::Float64
  financialcosts::Float64
  G::Float64
end

type Equilibrium
  r::Real #interest rate
  w::Real #wage
  #Results
  distr::Array{Float64,2}
  E::Float64
  a::Aggregates
  m::Moments
end

###########################################################################
# 0.PARAMETER DEFINITION

function init_hhparameters(bbeta=0.98,ssigma=1.0,psi=1)
  #Guess H such that labor supply in deterministic steady sate =1
  ddelta=0.14; aalphak=0.3; aalphal = 0.65;
  K= aalphak/((bbeta^(-1.0) -1 + ddelta));
  s= ddelta*K; #Savings
  c=1-s;
  H=aalphal/c;
  HouseholdParam(bbeta, ssigma, H, psi);
end

# Initialize firm parameters
function init_firmparameters(hp;aalphak=0.3, aalphal = 0.65, ff=0.0145, llambda0= 0.08, llambda1= 0.028, ddelta= 0.14, ttheta=0.4, kappa=1, e=0.01,rhoz= 0.76, ssigmaz= 0.0352, Nz::Int=9,  Nk::Int=80, Nq::Int=40, Nomega::Int=100)
  mc = tauchen(Nz,rhoz,ssigmaz); # Process of firm productivity z
  logshocks = mc.state_values;
  shocks=exp(logshocks);
  trans = mc.p;
  invariant=trans^100;
  invariant_dist=invariant[1,:];
  #aalpha= aalphak/(1-aalphal);
  #A = ( (aalpha)/(hp.beta^-1.0 - 1.0 + ddelta)  )^(-aalpha);
  A = ( (aalphak)/(hp.beta^-1.0 - 1.0 + ddelta)  )^(-aalphak);
  zgrid = A*shocks
  ztrans=trans';
  FirmParam(aalphak, aalphal, ff, llambda0, llambda1, ddelta, ttheta, kappa, e,ttheta*(1-ddelta), 1/(1-ttheta*(1-ddelta)),zgrid, ztrans,invariant_dist,Nz,Nk,Nq,Nomega);
end

#Initialize taxes
function init_taxes(;ttaud=0.15, ttauc = 0.35, ttaui = 0.3, ttaug = 0.15)
  Taxes(ttaud, ttauc, ttaui, ttaug);
  end

#Guess Prices
function guess_prices(tau,fp,hp)
  r=(hp.beta^(-1.0) -1)/(1-tau.i);
  w=0.69 #fp.alphal;

  #Initiate Results
  distr= Array(Float64,(fp.Nomega,fp.Nz));
  E=convert(Float64,NaN);
  netdistributions = convert(Float64,NaN);
  agginterests = convert(Float64,NaN);
  consumption = convert(Float64,NaN);
  output = convert(Float64,NaN);
  laborsupply= convert(Float64,NaN);
  welfare = convert(Float64,NaN);
  collections=Taxes(NaN,NaN,NaN,NaN);

  aggregates= Aggregates(NaN,NaN,NaN,NaN,NaN,NaN,collections,NaN,NaN,NaN,NaN,NaN,NaN)
  moments= Moments(NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN)

  Equilibrium(r,w,distr,E,aggregates,moments);
end

#Compute omega prime
function omegaprimefun(kprime::Real ,qprime::Real ,i_zprime::Int ,p::Equilibrium ,tau:: Taxes ,fp::FirmParam)
  zprime=fp.zgrid[i_zprime];
  lprime= (zprime*fp.alphal*kprime^fp.alphak / p.w)^(1/(1-fp.alphal));

  return (1-tau.c)*(zprime*kprime^fp.alphak*lprime^fp.alphal -p.w*lprime - fp.delta*kprime - p.r*qprime) + kprime - qprime +tau.c*fp.f
end

function grossdistributions(omega::Real,kprime::Real,qprime::Real,fp::FirmParam)
  omega - kprime + qprime -fp.f
end

#Predict future state
function predict_state(i_zprime::Int64, i_omega::Int64, i_z::Int64, pr::FirmProblem, tau:: Taxes, fp::FirmParam)
  kprime = pr.kpolicygrid[i_omega,i_z];
  qprime = pr.qpolicygrid[i_omega,i_z];
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
  return omegaprime, i_omegaprime
end
