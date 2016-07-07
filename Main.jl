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

type Param
  ##Household parameters
  beta::Real #Discount rate
  sigma::Real  #Risk aversion/ ies
  H::Real  # Labor supply parameter
  psi::Real #Labor supply

  ## Firm parameters
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
  invariant_distr::Array{Float64,1} #Invariant distribution
  Nz::Int64
  Nk::Int64
  Nq::Int64
  Nomega::Int64

  omega::GridObject
  kprime::GridObject
  qprime::GridObject
end

type Taxes
  d::Real
  c::Real
  i::Real
  g::Real
end


type FirmProblem
  #constants
  betatilde::Real
  taudtilde::Real
  discounted_interest::Real # betatilde*(1+(1-tau.c)*r)

  #input
  firmvalueguess::Matrix

  #output
  firmvaluegrid::Matrix
  kpolicy:: Matrix
  qpolicy::Matrix

  InterpolationGrid::Array{CoordInterpGrid,1}

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

# Initialize parameters
function init_parameters(;bbeta=0.98,ssigma=1.0,psi=1,aalphak::Float64=0.3, aalphal::Float64 = 0.65, ff::Float64=0.0145, llambda0::Float64= 0.08, llambda1::Float64= 0.028, ddelta::Float64= 0.14, ttheta::Float64=0.45, kappa::Float64=1.0, e::Float64=0.00,rhoz::Float64= 0.76, ssigmaz::Float64= 0.0352, Nz::Int64=9, Nk::Int64=80, Nq::Int64=40, Nomega::Int64=100, A::Float64=0.76)
  #Guess H such that labor supply in deterministic steady sate =1
  K= aalphak/((bbeta^(-1.0) -1 + ddelta));
  s= ddelta*K; #Savings
  c=1-s;
  H=aalphal/c;

  mc = tauchen(Nz,rhoz,ssigmaz); # Process of firm productivity z
  logshocks = mc.state_values;
  shocks=exp(logshocks);
  trans = mc.p;
  invariant=trans^100;
  invariant_dist=collect(invariant[1,:]);
  zgrid = A*shocks;
  ztrans=trans';

  #construct grid for capital
  wage=aalphal;
  auxconst= (1/wage)^(aalphal/(1-aalphal))*(aalphal^(aalphal/(1-aalphal)) - aalphal^(1/(1-aalphal)));
  ggamma = zgrid.^(1/(1-aalphal)).*auxconst;
  maxexpgamma =ztrans[:,end]'*ggamma;
  kmax = (  ((aalphak/(1-aalphal))*maxexpgamma[1] )/(bbeta^-1.0 -1 +ddelta )  )^((1.0 - aalphal)/(1.0 - aalphak - aalphal));
  kub=1.05*kmax; 
  klb = 0.0;
  kstep = (kub-klb)/(Nk-1);
  kgrid = 0:kstep:kub;
  kprime=GridObject(kub, 0, kstep, Nk, kgrid);

  #construct grid for debt
  qub = ttheta*(1-ddelta)*kub ;
    #lower bound is such that there is enough cash to finance max investment when a max(z) shock follows two min(z) shocks
    zprime=zgrid[1];
    minexpgamma =ztrans[:,1]'*ggamma;
    k1 = (  ((aalphak/(1-aalphal))*minexpgamma[1] )/(bbeta^-1.0 -1 +ddelta )  )^((1.0 - aalphal)/(1.0 - aalphak - aalphal));
    lprime= (zprime*aalphal*k1^aalphak / wage)^(1/(1-aalphal));
  qlb= -(kmax - 0.4*(zprime*k1^aalphak*lprime^aalphal -wage*lprime - ddelta*k1 -ff) +k1)/(bbeta^-1.0 -1);
  qstep = (qub-qlb)/(Nq-1);
  qgrid = qlb:qstep:qub;
  qprime=GridObject(qub,qlb,qstep,Nq,qgrid);

  #grid for net worth
    zprime=zgrid[end];
    lprime= (zprime*aalphal*kub^aalphak / wage)^(1/(1-aalphal));
  omegaub = zprime*kub^aalphak*lprime^aalphal -wage*lprime + (1-ddelta)*kub -ff;
  omegalb=0;
  omegastep=(omegaub-omegalb)/(Nomega-1);
  omegagrid = omegalb:omegastep:omegaub; #collect(linspace(0,1.2*kmax,Nomega));
  Nomega=length(omegagrid)
  omega=GridObject(omegaub, omegalb, omegastep,Nomega, omegagrid);

  Param(bbeta, ssigma, H, psi, aalphak, aalphal, ff, llambda0, llambda1, ddelta, ttheta, kappa, e,ttheta*(1-ddelta), 1/(1-ttheta*(1-ddelta)),zgrid, ztrans,invariant_dist,Nz,Nk,Nq,Nomega, omega, kprime, qprime);
end

#Initialize taxes
function init_taxes(;ttaud::Float64 =0.15, ttauc::Float64 = 0.35, ttaui::Float64 = 0.3, ttaug::Float64 = 0.15)
  Taxes(ttaud, ttauc, ttaui, ttaug);
  end

#Guess Prices
function init_equilibirium(wguess::Float64,tau::Taxes,pa::Param)
  r=(pa.beta^(-1.0) -1)/(1-tau.i);
  w=wguess; #pa.alphal;

  #Initiate Results
  distr= Array(Float64,(pa.Nomega,pa.Nz));
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


#Initialize firm problem
function init_firmproblem(eq::Equilibrium, tau::Taxes, pa::Param ; guessvalue = false, firmvalueguess::Matrix= ones(pa.Nomega,pa.Nz))

  #Nk, Nq, Nz, Nomega = pa.Nk, pa.Nq, pa.Nz, pa.Nomega

  betatilde = (1.0 + (1-tau.i)/(1-tau.g)*eq.r )^(-1);
  taudtilde = 1-(1-tau.d)/(1-tau.g);

  #guess firm value
  if !guessvalue
    firmvalueguess = repmat(pa.omega.grid,1,pa.Nz);
  end
  firmvaluegrid  = copy(firmvalueguess);
  kpolicy    = similar(firmvaluegrid);
  qpolicy    = similar(firmvaluegrid);

  #Preallocate results
  distributions = Array(Float64,(pa.Nomega,pa.Nz));
  financialcosts = zeros(Float64,(pa.Nomega,pa.Nz));
  grossdividends = zeros(Float64,(pa.Nomega,pa.Nz));
  grossequityis = zeros(Float64,(pa.Nomega,pa.Nz));
  exitprobability = Array(Float64,(pa.Nomega,pa.Nz));
  exitrule = falses(pa.Nomega,pa.Nz,pa.Nz);
  positivedistributions = falses(pa.Nomega,pa.Nz);

  #Initialize the firm problem object
  FirmProblem( betatilde, taudtilde, betatilde*(1+(1-tau.c)*eq.r), firmvalueguess, firmvaluegrid, kpolicy, qpolicy,
   map(x->CoordInterpGrid(pa.omega.grid,firmvalueguess[:,x],BCnearest, InterpLinear),1:pa.Nz),
   distributions, financialcosts, grossdividends, grossequityis, exitprobability, exitrule, positivedistributions);
end

#Compute omega prime
function omegaprimefun(kprime::Real, qprime::Real, i_zprime::Int, eq::Equilibrium, tau:: Taxes, pa::Param)
  zprime=pa.zgrid[i_zprime];
  lprime= (zprime*pa.alphal*kprime^pa.alphak / eq.w)^(1/(1-pa.alphal));

  return (1-tau.c)*(zprime*kprime^pa.alphak*lprime^pa.alphal -eq.w*lprime - pa.delta*kprime - eq.r*qprime) + kprime - qprime +tau.c*pa.f
end

function grossdistributions(omega::Real,kprime::Real,qprime::Real,pa::Param)
  omega - kprime + qprime -pa.f
end

#Predict future state
function predict_state(i_zprime::Int64, i_omega::Int64, i_z::Int64, pr::FirmProblem, eq::Equilibrium, tau:: Taxes, pa::Param)
  kprime = pr.kpolicy[i_omega,i_z];
  qprime = pr.qpolicy[i_omega,i_z];
  omegaprime   = omegaprimefun(kprime,qprime,i_zprime,eq,tau,pa);
  i_omegaprime = closestindex(omegaprime, pa.omega.step);
  #The block below checks that the index is within reasonable bounds
  if i_omegaprime<1 || i_omegaprime>pa.Nomega
    if i_omega==pr.Nomega || i_omegaprime < (pa.Nomega + 3)
      i_omegaprime =pa.Nomega
    elseif i_omega==1 || i_omegaprime > -3
      i_omegaprime =1;
    else
      error("omega' out of the grid ", "i_z = ", i_z, "i_omega = ", i_omega, "i_z' = ",i_zprime, "i_omega' = ",i_omegaprime)
    end
  end
  return omegaprime, i_omegaprime
end
