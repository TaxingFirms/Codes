function plotpolicies(pr::FirmProblem, pa::Param)
  for i=1:pa.Nz
  figure()
  d= plot(pa.omega.grid, pr.distributions[:,i] )
  k= plot(pa.omega.grid, pr.kpolicy[:,i] )
  q= plot(pa.omega.grid, pr.qpolicy[:,i] )
    xlabel("Net worth")
    title("Policy functions at z=%f")
    legend("dkq", loc="best")
  end
end

function plotsimulations(pr::FirmProblem, pa::Param)
  capital, debt, networth, dividends, investment, z_history_ind = simulation(50,50,pr,pa);
  figure()
  plot(capital)
  figure()
  plot(debt)
  figure()
  plot(dividends)
end

function plotvaluefcn(pr::FirmProblem,pa::Param)
  figure()
  for i=1:pa.Nz
   plot(pa.omega.grid, pr.firmvaluegrid[:,i] )
    xlabel("Net worth")
    title("Value function")
  end
end
