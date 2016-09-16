
function plotsimulations(pr::FirmProblem, pa::Param)
  capital, debt, networth, dividends, investment, z_history_ind = simulation(50,50,pr,pa);
  figure()
  plot(capital)
  figure()
  plot(debt)
  figure()
  plot(dividends)
end


function plotpolicies(pr::FirmProblem,pa::Param)
  for j=1:pa.Nz
    figure()
    d= plot(pa.omega.grid, pr.distributions[:,j] )
    k= plot(pa.omega.grid, pr.kpolicy[:,j] )
    q= plot(pa.omega.grid, pr.qpolicy[:,j] )
      xlabel("Net worth")
    tit = string("Policy functions (at z= ",j,")")
      title(tit)
      legend("dkq", loc="best")
  end
end


function plotkpolicy(pr::FirmProblem,pa::Param)
  for j=1:pa.Nz
    figure()
    k= plot(pa.omega.grid, pr.kpolicy[:,j] )
      xlabel("Net worth")
      ylabel("""k' """)
    tit = string("Capital (at z= ",j,")")
      title(tit)
  end
end


function plotkpolicyshift(pr::FirmProblem,pr1::FirmProblem,pa::Param)
  for j=1:pa.Nz
    figure()
    k= plot(pa.omega.grid, pr.kpolicy[:,j])
    k1= plot(pa.omega.grid, pr1.kpolicy[:,j])
      xlabel("Net worth")
      ylabel("""k' """)
    tit = string("Capital (at z= ",j,")")
      title(tit)
  end
end
