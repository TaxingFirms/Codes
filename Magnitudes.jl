function magnitudes(tau::Taxes, pa::Param; r::Float64 = (pa.beta^(-1.0) -1.0)/(1.0-tau.i), w::Float64=pa.alphal)


  #In steady state, r is given by,
  # r= (pa.beta^(-1.0) -1.0)/(1.0-tau.i);  w=pa.alphal;

  gamma = pa.alphak/(1.0-pa.alphal);
  ppsi =  w^(-pa.alphal/(1-pa.alphal))*(pa.alphal^(pa.alphal/(1-pa.alphal)) - pa.alphal^(1/(1-pa.alphal)));
  betatilde = (1.0 + (1.0-tau.i)/(1.0-tau.g)*r )^(-1.0);

  userscost = (betatilde^-1.0-1.0)/(1.0-tau.c)+pa.delta;
  capital_unc = Array(Float64,(pa.Nz,));

  for i_z=1:pa.Nz
    # For an unconstrained firm
    expected_z= pa.ztrans[:,i_z]'*(pa.zgrid.^(1.0/(1.0-pa.alphal)) )* ppsi;
    capital_unc[i_z] = (gamma*expected_z[1]/userscost)^(1.0/(1.0-gamma));
  end

  #conditional on issuing equity, users's cost becomes.

  userscost2 = ((1.0+pa.lambda1)/betatilde -1.0)/(1.0-tau.c)+pa.delta;
  capital_equity = Array(Float64,(pa.Nz,));

  for i_z=1:pa.Nz
    expected_z= pa.ztrans[:,i_z]'*(pa.zgrid.^(1.0/(1.0-pa.alphal)) )* ppsi;
    capital_equity[i_z] = (gamma*expected_z[1]/userscost2)^(1.0/(1.0-gamma));
  end

  #Profits at unconstrained capital given same shock
  profits_unc = Array(Float64,(pa.Nz,));
  for i_z = 1:pa.Nz
    labor= (pa.zgrid[i_z]*pa.alphal*capital_unc[i_z]^pa.alphak / w)^(1/(1-pa.alphal));
    profits_unc[i_z]=(1-tau.c)*(pa.zgrid[i_z]*capital_unc[i_z]^pa.alphak*labor^pa.alphal -w*labor - pa.delta*capital_unc[i_z] - pa.f);
  end

  #Profits at equity issued capital
  profits_equity = Array(Float64,(pa.Nz,));
  for i_z = 1:pa.Nz
    labor= (pa.zgrid[i_z]*pa.alphal*capital_equity[i_z]^pa.alphak / w)^(1/(1-pa.alphal));
    profits_equity[i_z]=(1-tau.c)*(pa.zgrid[i_z]*capital_equity[i_z]^pa.alphak*labor^pa.alphal -w*labor - pa.delta*capital_equity[i_z] - pa.f);
  end
  capital_unc, capital_equity, profits_unc, profits_equity
end
