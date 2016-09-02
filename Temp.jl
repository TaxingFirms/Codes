rpr,req,rtau = taxreform2(0.1, eq, tau, pa; tol=10.0^-2.0,update=0.7, maxroutine=maximizationfast);


using JLD
type Economy
  pr::FirmProblem
  eq::Equilibrium
  tau::Taxes
  cev::Float64
end
ref,pa=load("Reform2Benchmark.jld", "ref","pa");
refhe,pahe=load("Reform2HighElasticity.jld", "ref","pa");
refla,pala= load("Reform2LowAllowance.jld", "ref","pa");


for j=1:6
  println(ref[j].eq.a.consumption)
end


cevarray= Array(Float64,(6,));
for j=1:6
  cevarray[j]= (ref[j].eq.a.consumption - ref[1].eq.a.consumption)/ref[1].eq.a.consumption - pa.H/(ref[1].eq.a.consumption*(1+1/pa.psi))*( (ref[j].eq.w*(1-ref[j].tau.l)/pa.H)^(1+pa.psi) - (ref[1].eq.w*(1-ref[1].tau.l)/pa.H)^(1+pa.psi) );
  println(cevarray[j])
end

i=7
figure()
d= plot(pa.omega.grid, pr0.distributions[:,i] )
k= plot(pa.omega.grid, pr0.kpolicy[:,i] )
  xlabel("Net worth")
  title("Policy functions at z=%f")
  legend("dkq", loc="best")

legend("dkq", loc="best")
d2= plot(pa.omega.grid, pr1.distributions[:,i] )
k2= plot(pa.omega.grid, pr1.kpolicy[:,i] )
  xlabel("Net worth")
  title("Policy functions at z=%f")
  legend("dkq", loc="best")


  i=4
  figure()
  d= plot(pa.omega.grid, pr0.distributions[:,i] )
  k= plot(pa.omega.grid, pr0.kpolicy[:,i] )
    xlabel("Net worth")
    title("Policy functions at z=%f")
    legend("dkq", loc="best")

  legend("dkq", loc="best")
  d2= plot(pa.omega.grid, pr1.distributions[:,i] )
  k2= plot(pa.omega.grid, pr1.kpolicy[:,i] )
    xlabel("Net worth")
    title("Policy functions at z=%f")
    legend("dkq", loc="best")
