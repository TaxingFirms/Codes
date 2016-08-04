function magnitudes(tau::Taxes, pa::Param)

  optimalcaps = Array(Float64,(pa.Nz,));

  #In steady state, r is given by,
  r=(pa.beta^(-1.0) -1.0)/(1.0-tau.i);
  w=pa.alphal;

  gamma = pa.alphak/(1.0-pa.alphal);
  ppsi =  (1/w)^(aalphal/(1-aalphal))*(aalphal^(aalphal/(1-aalphal)) - aalphal^(1/(1-aalphal)))


  betatilde = (1.0 + (1.0-tau.i)/(1.0-tau.g)*r )^(-1.0);
  userscost = (betatilde^-1.0-1.0)/(1.0-tau.c)+pa.delta;
  for i_z=1:pa.Nz
    # For an unconstrained firm
    expected_z= ztrans[:,i_z]'*(zgrid.^(1.0/(1.0-pa.alphal)) )* ppsi;
    optimalcaps[i_z] = (gamma*expected_z[1]/userscost)^(1.0/(1.0-gamma));
  end


  #conditional on issuing equity, users's cost becomes.

  userscost2 = ((1.0+pa.lambda1)/betatilde -1.0)/(1.0-tau.c)+pa.delta;
  capital_equity = Array(Float64,(pa.Nz,));

  for i_z=1:pa.Nz
    expected_z= ztrans[:,i_z]'*(zgrid.^(1.0/(1.0-pa.alphal)) )* (1.0-pa.alphal);
    capital_equity[i_z] = (gamma*expected_z[1]/userscost2)^(1.0/(1.0-gamma));
  end

  #Profits at unconstrained capital given same shock
  profits_unc = Array(Float64,(pa.Nz,));
  for i_z = 1:pa.Nz
    labor= (pa.zgrid[i_z]*pa.alphal*optimalcaps[i_z]^pa.alphak / w)^(1/(1-pa.alphal));
    profits_unc[i_z]=pa.zgrid[i_z]*optimalcaps[i_z]^pa.alphak*labor^pa.alphal -w*labor - pa.f;
  end

  #Profits at equity issued capital
  profits_equity = Array(Float64,(pa.Nz,));
  for i_z = 1:pa.Nz
    labor= (pa.zgrid[i_z]*pa.alphal*optimalcaps[i_z]^pa.alphak / w)^(1/(1-pa.alphal));
    profits_equity[i_z]=pa.zgrid[i_z]*optimalcaps[i_z]^pa.alphak*labor^pa.alphal -w*labor - pa.f;
  end



  optimalcaps, capital_equity
end
