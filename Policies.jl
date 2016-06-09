

#copies the optimal policies in a new object after the VFI is over, just to organize thing a bit.
#also allocates memory for get_policies
function copy_opt_policies(pr::FirmProblem)
  firmvalue= pr.firmvaluegrid;
  interpolation = pr.InterpolationGrid;
  kprime = pr.kpolicygrid;
  qprime = pr.qpolicygrid;

  distributions = Array(Float64,(pr.Nomega,pr.Nz));
  financialcosts = zeros(Float64,(pr.Nomega,pr.Nz));
  grossdividends = zeros(Float64,(pr.Nomega,pr.Nz));
  grossequityis = zeros(Float64,(pr.Nomega,pr.Nz));
  exitprobability = Array(Float64,(pr.Nomega,pr.Nz));
  exitrule = falses(pr.Nomega,pr.Nz,pr.Nz);
  positivedistributions = falses(pr.Nomega,pr.Nz);

  ResultsFP(firmvalue, interpolation,kprime,qprime, distributions, financialcosts, grossdividends, grossequityis, exitprobability, exitrule, positivedistributions)
end

#Gets the exit rules, distributions and other quantities of interest

function getpolicies!(res::ResultsFP, pr::FirmProblem, p::Equilibrium, tau::Taxes, fp::FirmParam)
  for (i_z,z) in enumerate(fp.zgrid)
    for (i_omega, omega) in enumerate(pr.omega.grid)
      kprime= res.kprime[i_omega, i_z];
      qprime= res.qprime[i_omega, i_z];

      if grossdistributions(omega,kprime,qprime,fp) >=0
        res.positivedistributions[i_omega, i_z]= true;
        res.distributions[i_omega, i_z]= (1-tau.d)*grossdistributions(omega,kprime,qprime,fp);
        res.grossdividends[i_omega, i_z]= grossdistributions(omega,kprime,qprime,fp);
      else
        res.distributions[i_omega, i_z]= (1+fp.lambda1)*grossdistributions(omega,kprime,qprime,fp) - fp.lambda0;
        res.financialcosts[i_omega, i_z]= fp.lambda1*grossdistributions(omega,kprime,qprime,fp) -fp.lambda0;
        res.grossequityis[i_omega, i_z]= grossdistributions(omega,kprime,qprime,fp);
      end



      prexit=0.0;
       #This just avoids computing this constant again and again while computing omegaprime
      exitvalue = (1-pr.taudtilde)*(fp.kappa*(1-fp.delta)*kprime - (1+p.r)*qprime);
      for (i_zprime,zprime) in enumerate(fp.zgrid)
        lprime=(fp.alphal*zprime*(kprime^fp.alphak)/p.w)^(1/(1-fp.alphal));
        omegaprime = omegaprimefun(kprime,qprime,i_zprime,p,tau,fp);
        contvalue = firmvaluefunction(omegaprime,i_zprime,pr);
        if exitvalue>=contvalue #if indiferent, firms choose to exit
          res.exitrule[i_omega, i_z,i_zprime]=true;
          prexit+=fp.ztrans[i_zprime,i_z];
        end
      end
      res.exitprobability[i_omega, i_z]=prexit;
    end
  end
end



