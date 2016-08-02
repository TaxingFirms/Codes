function simulation(S::Int64, T::Int64,pr::FirmProblem,pa::Param; seed::Int64 =1234)
  #Interpolate policy functions
  kprimefun = map(x->CoordInterpGrid(pa.omega.grid,pr.kpolicy[:,x],BCnearest, InterpLinear),1:pa.Nz);
  qprimefun = map(x->CoordInterpGrid(pa.omega.grid,pr.qpolicy[:,x],BCnearest, InterpLinear),1:pa.Nz);
  grossdividendfun= map(x->CoordInterpGrid(pa.omega.grid,pr.grossdividends[:,x],BCnearest, InterpLinear),1:pa.Nz);
  grossequitfun = map(x->CoordInterpGrid(pa.omega.grid,pr.grossequityis[:,x],BCnearest, InterpLinear),1:pa.Nz);
  z_history_ind = generate_shocks(seed,S, T, pa);

  # Declare variables
  capital = Array(Float64,(T,S));
  debt = Array(Float64,(T,S));
  dividends = Array(Float64,(T,S));
  networth = Array(Float64,(T,S));

  for i_s=1:S
    networth[1,i_s]=omegaprimefun(pa.k0, 0.0, z_history_ind[1,i_s], eq, tau, pa);
    capital[1,i_s]= pa.k0;
    debt[1,i_s]=0.0;
    dividends= grossdividendfun[networth[1,i_s],z_history_ind[1,i_s]];

    for i_t=1:T
    end
  end

end





function generate_shocks(seed::Int64,S::Int64, T::Int64, pa::Param)
  #Set random seed
  srand(seed);
  #Generate S simulation of T draw, uniform in [0,1)
  random_uniform = rand(Float64,T,S);
  z_history_ind = Array(Int64,(T,S));
  #Transform first random uniform into draw from ergodic distribution
    #Generate thresholds = normcdf(midpoints)
    thresholds_invariant, thresholds_conditional = generate_thresholds(pa);
    #Map uniform random to point in grid of shocks

    for i_s =1:S
      z_history_ind[1,i_s] = inverse_cumulative(thresholds_invariant,random_uniform[1,i_s], pa.Nz)
      for i_t=2:T
        z_ind = z_history_ind[i_t-1,i_s];
        z_history_ind[i_t,i_s] = inverse_cumulative(thresholds_conditional[:,z_ind],random_uniform[i_t,i_s], pa.Nz);
      end
    end
    return z_history_ind
end


function inverse_cumulative(thresholds::Array{Float64,1},random::Float64, Nz)
  ind_shock=10;
  for i_z = 1:Nz
    if thresholds[i_z] <= random < thresholds[i_z+1]
      ind_shock = i_z;
    end
  end
  ind_shock
end

function generate_thresholds(pa::Param)
  # Generate thresholds to invert uniform draws
  # 1. Invariant distributuon
  mc = tauchen(pa.Nz,pa.rhoz,pa.sigmaz);
  thresholds_invariant= Array(Float64,(pa.Nz+1,));
  thresholds_invariant[1]=0.0
  for i=1:(pa.Nz-1)
    midpoint= (mc.state_values[i] + mc.state_values[i+1])/2;
    thresholds_invariant[i+1] = normcdf(0,sqrt(pa.sigmaz^2.0/(1-pa.rhoz^2.0)), midpoint); # The variance of the invariant distribution is ssigmaz^2.0/(1-rhoz^2.0)
  end
  thresholds_invariant[pa.Nz+1]=1.0;
  # 2. Conditional on current z
  thresholds_conditional= Array(Float64,(pa.Nz+1,pa.Nz));
  for i_z=1:pa.Nz
    thresholds_conditional[1,i_z]=0.0;
    for i=1:(pa.Nz-1)
      midpoint= (mc.state_values[i] + mc.state_values[i+1])/2;
      thresholds_conditional[i+1,i_z] = normcdf(pa.rhoz*mc.state_values[i_z], pa.sigmaz, midpoint); # The variance of the invariant distribution is ssigmaz^2.0/(1-rhoz^2.0)
    end
    thresholds_conditional[pa.Nz+1,i_z]=1.0;
  end
return thresholds_invariant, thresholds_conditional
end
