function comparativestatics( eq::Equilibrium, pr:: FirmProblem, tau::Taxes , pa::Param; tauvec::Array{Float64,1}=[0.15, 0.20, 0.25, 0.30, 0.35, 0.40])
  # INPUT: initial equilibrium
  # OUTPUT: matrix of tax basis
  Ntau,=size(tauvec);
  compstats = Array(Economy,(Ntau,4));

  eqvoid = init_equilibirium(0.0,tau,pa);
  prvoid  = init_firmproblem(pa);

  #DIVIDENDS
  for j=1:Ntau
    tau0 = Taxes(tauvec[j], tau.c,tau.i,tau.g,tau.l)
    if (1-tau0.c)*(1-tau0.g)>(1-tau0.i)
      compstats[j,1]= Economy(prvoid, eqvoid, tau0, 0.0)
    else
      pr0, eq0 = SolveSteadyState(tau0,pa);
      compstats[j,1]= Economy(pr0, eq0, tau0, 0.0);
    end
  end
  save("compstats.jld","compstats",compstats,"pa",pa);

  #CORPORATE
  for j=1:Ntau
    tau0 = Taxes(tau.d, tauvec[j],tau.i,tau.g,tau.l)
    if (1-tau0.c)*(1-tau0.g)>(1-tau0.i)
      compstats[j,2]= Economy(prvoid, eqvoid, tau0, 0.0)
    else
      pr0, eq0 = SolveSteadyState(tau0,pa)
      compstats[j,2]= Economy(pr0, eq0, tau0, 0.0)
    end
  end
  save("compstats.jld","compstats",compstats,"pa",pa);

  #INTEREST
  for j=1:Ntau
    tau0 = Taxes(tau.d, tau.c,tauvec[j],tau.g,tau.l)
    if (1-tau0.c)*(1-tau0.g)>(1-tau0.i)
      compstats[j,3]= Economy(prvoid, eqvoid, tau0, 0.0)
    else
      pr0, eq0 = SolveSteadyState(tau0,pa)
      compstats[j,3]= Economy(pr0, eq0, tau0, 0.0)
    end
  end
  save("compstats.jld","compstats",compstats,"pa",pa);

  #CAPGAINS
  for j=1:Ntau
    println(j,j,j,j,j,j,j,j,j,j,j,j,j,j,j,j,j,j,j,j,j,j,j,j,j,j,j,j,j,j,j,j)
    tau0 = Taxes(tau.d, tau.c,tau.i,tauvec[j],tau.l)
    if (1-tau0.c)*(1-tau0.g)>(1-tau0.i)
      compstats[j,4]= Economy(prvoid, eqvoid, tau0, 0.0)
    else
      pr0, eq0 = SolveSteadyState(tau0,pa)
      compstats[j,4]= Economy(pr0, eq0, tau0, 0.0)
    end
  end
  save("compstats.jld","compstats",compstats,"pa",pa);
end


function drawtaxbases()
  compstats,pa=load("compstats.jld", "compstats","pa");
  Ntau,~=size(compstats)

  taxbases2div=Array(Float64,(Ntau,4));
  taxbases2corp=Array(Float64,(Ntau,4));
  taxbases2int=Array(Float64,(Ntau,4));
  taxbases2cg=Array(Float64,(Ntau,4));

  tauvec=[0.15, 0.20, 0.25, 0.30, 0.35, 0.40]
  for i=1:Ntau
    taxbases2div[i,1] = compstats[i,1].eq.a.collections.d/compstats[i,1].tau.d
    taxbases2div[i,2] = compstats[i,1].eq.a.collections.c/compstats[i,1].tau.c
    taxbases2div[i,3] = compstats[i,1].eq.a.collections.i/compstats[i,1].tau.i
    taxbases2div[i,4] = compstats[i,1].eq.a.collections.l/compstats[i,1].tau.l
  end

  for i=1:Ntau
    taxbases2corp[i,1] = compstats[i,2].eq.a.collections.d/compstats[i,2].tau.d
    taxbases2corp[i,2] = compstats[i,2].eq.a.collections.c/compstats[i,2].tau.c
    taxbases2corp[i,3] = compstats[i,2].eq.a.collections.i/compstats[i,2].tau.i
    taxbases2corp[i,4] = compstats[i,2].eq.a.collections.l/compstats[i,2].tau.l
  end

  for i=1:Ntau
    taxbases2int[i,1] = compstats[i,3].eq.a.collections.d/compstats[i,3].tau.d
    taxbases2int[i,2] = compstats[i,3].eq.a.collections.c/compstats[i,3].tau.c
    taxbases2int[i,3] = compstats[i,3].eq.a.collections.i/compstats[i,3].tau.i
    taxbases2int[i,4] = compstats[i,3].eq.a.collections.l/compstats[i,3].tau.l
  end

  for i=1:1
    taxbases2cg[i,1] = compstats[i,4].eq.a.collections.d/compstats[i,4].tau.d
    taxbases2cg[i,2] = compstats[i,4].eq.a.collections.c/compstats[i,4].tau.c
    taxbases2cg[i,3] = compstats[i,4].eq.a.collections.i/compstats[i,4].tau.i
    taxbases2cg[i,4] = compstats[i,4].eq.a.collections.l/compstats[i,4].tau.l
  end

  figure()
  plot(tauvec, taxbases2div)
  xlabel("Dividend Tax")
  ylabel("Tax Base")
  legend(["dividends","corporare","interest", "labor"], loc="best")

  figure()
  plot(tauvec, taxbases2corp)
  xlabel("Corporate Tax")
  ylabel("Tax Base")
  legend(["dividends","corporare","interest", "labor"], loc="best")

  figure()
  plot(tauvec, taxbases2int)
  xlabel("Interest Tax")
  ylabel("Tax Base")
  legend(["dividends","corporare","interest", "labor"], loc="best")
end

##################################################
##################################################
##################################################

function drawtaxcollections()
  compstats,pa=load("compstats.jld", "compstats","pa");
  Ntau,~=size(compstats)

  taxclctn2div=Array(Float64,(Ntau,4));
  taxclctn2corp=Array(Float64,(Ntau,4));
  taxclctn2int=Array(Float64,(Ntau,4));
  taxclctn2cg=Array(Float64,(Ntau,4));

  tauvec=[0.15, 0.20, 0.25, 0.30, 0.35, 0.40]
  for i=1:Ntau
    taxclctn2div[i,1] = compstats[i,1].eq.a.collections.d
    taxclctn2div[i,2] = compstats[i,1].eq.a.collections.c
    taxclctn2div[i,3] = compstats[i,1].eq.a.collections.i
    taxclctn2div[i,4] = compstats[i,1].eq.a.collections.l
  end

  for i=1:Ntau
    taxclctn2corp[i,1] = compstats[i,2].eq.a.collections.d
    taxclctn2corp[i,2] = compstats[i,2].eq.a.collections.c
    taxclctn2corp[i,3] = compstats[i,2].eq.a.collections.i
    taxclctn2corp[i,4] = compstats[i,2].eq.a.collections.l
  end

  for i=1:Ntau
    taxclctn2int[i,1] = compstats[i,3].eq.a.collections.d
    taxclctn2int[i,2] = compstats[i,3].eq.a.collections.c
    taxclctn2int[i,3] = compstats[i,3].eq.a.collections.i
    taxclctn2int[i,4] = compstats[i,3].eq.a.collections.l
  end

#  for i=1:1
#    taxclctn2cg[i,1] = compstats[i,4].eq.a.collections.d
#    taxclctn2cg[i,2] = compstats[i,4].eq.a.collections.c
#    taxclctn2cg[i,3] = compstats[i,4].eq.a.collections.i
#    taxclctn2cg[i,4] = compstats[i,4].eq.a.collections.l
#  end

  figure()
  plot(tauvec, taxclctn2div)
  xlabel("Dividend Tax")
  ylabel("Tax Collections")
  legend(["dividends","corporare","interest", "labor"], loc="best")

  figure()
  plot(tauvec, taxclctn2corp)
  xlabel("Corporate Tax")
  ylabel("Tax Collections")
  legend(["dividends","corporare","interest", "labor"], loc="best")

  figure()
  plot(tauvec, taxclctn2int)
  xlabel("Interest Tax")
  ylabel("Tax Collections")
  legend(["dividends","corporare","interest", "labor"], loc="center")
end
