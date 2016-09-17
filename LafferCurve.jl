
immutable EquilArgument
  tau::Taxes
  pa::Param
  wguess::Float64
end

type AllInfo
  pr::FirmProblem
  eq::Equilibrium
end

function compute_equilibrium(arg::EquilArgument)
  pr,eq= SolveSteadyState(arg.tau,arg.pa;wguess=arg.wguess, VFIfunction= firmVFI! ,VFItol=10.0^-3.0, displayit0=false, displayw = false);
end

function laffer_tauc(taucvec::FloatRange{Float64},filename::ASCIIString)
  pr,eq,tau,pa=load("ModelResults.jld", "pr","eq","tau","pa");
  benchmarkG=eq.a.G;

  Ntau=size(taucvec)[1]
  collections = Array(Float64,(Ntau,));
  collections2gdp = Array(Float64,(Ntau,));

  allinfo = Array(AllInfo,(Ntau,));

  for j=1:Ntau
    println(j,j,j,j,j,j,j,j,j,j,j,j,j,j,j)
    pr,eq= SolveSteadyState(Taxes(tau.d,taucvec[j],tau.i,tau.g,tau.l),pa;wguess=eq.w, VFItol=10.0^-3.0, displayit0=false, displayw = false);
    collections[j]= eq.a.G/benchmarkG
    collections2gdp[j]= eq.a.G/eq.a.output
    allinfo[j]= AllInfo(pr,eq)
    println(eq.a.G/benchmarkG)
  end

  save(filename,"collections",collections,"collections2gdp", collections2gdp,"allinfo",allinfo)
end


function laffer_tauc_parallel(taucvec::FloatRange{Float64},filename::ASCIIString,benchmark::ASCIIString)

  Ntau=size(taucvec)[1]
  arguments = Array(EquilArgument,(Ntau,));

  pr,eq,tau,pa=load(benchmark, "pr","eq","tau","pa");

  for j=1:Ntau
    tauhat=Taxes(tau.d,taucvec[j],tau.i,tau.g,tau.l)
    arguments[j]= EquilArgument(tauhat,pa,eq.w)
  end

  bigvector=pmap(compute_equilibrium,arguments)

  benchmarkG=eq.a.G;
  collections = Array(Float64,(Ntau,));
  collections2gdp = Array(Float64,(Ntau,));
  allinfo = Array(AllInfo,(Ntau,));

  for j=1:Ntau
    pr,eq=bigvector[j]
    collections[j]= eq.a.G/benchmarkG
    collections2gdp[j]= eq.a.G/eq.a.output
    allinfo[j]= AllInfo(pr,eq)
    println(eq.a.G/benchmarkG)
  end


  save(filename,"collections",collections,"collections2gdp", collections2gdp,"allinfo",allinfo)
end


function laffer_taue_parallel(tauevec::FloatRange{Float64},filename::ASCIIString,benchmark::ASCIIString)

  Ntau=size(tauevec)[1]
  arguments = Array(EquilArgument,(Ntau,));

  pr,eq,tau,pa=load(benchmark, "pr","eq","tau","pa");

  for j=1:Ntau
    tauhat=Taxes(tauevec[j],tau.c,tau.i,tauevec[j],tau.l)
    arguments[j]= EquilArgument(tauhat,pa,eq.w)
  end

  bigvector=pmap(compute_equilibrium,arguments)

  benchmarkG=eq.a.G;
  collections = Array(Float64,(Ntau,));
  collections2gdp = Array(Float64,(Ntau,));
  allinfo = Array(AllInfo,(Ntau,));

  for j=1:Ntau
    pr,eq=bigvector[j]
    collections[j]= eq.a.G/benchmarkG
    collections2gdp[j]= eq.a.G/eq.a.output
    allinfo[j]= AllInfo(pr,eq)
    println(eq.a.G/benchmarkG)
  end


  save(filename,"collections",collections,"collections2gdp", collections2gdp,"allinfo",allinfo)
end


function laffer_taui_parallel(tauivec::FloatRange{Float64},filename::ASCIIString,benchmark::ASCIIString)

  Ntau=size(tauivec)[1]
  arguments = Array(EquilArgument,(Ntau,));

  pr,eq,tau,pa=load(benchmark, "pr","eq","tau","pa");

  for j=1:Ntau
    tauhat=Taxes(tau.d,tau.c,tauivec[j],tau.g,tau.l)
    arguments[j]= EquilArgument(tauhat,pa,eq.w)
  end

  bigvector=pmap(compute_equilibrium,arguments)

  benchmarkG=eq.a.G;
  collections = Array(Float64,(Ntau,));
  collections2gdp = Array(Float64,(Ntau,));
  allinfo = Array(AllInfo,(Ntau,));

  for j=1:Ntau
    pr,eq=bigvector[j]
    collections[j]= eq.a.G/benchmarkG
    collections2gdp[j]= eq.a.G/eq.a.output
    allinfo[j]= AllInfo(pr,eq)
    println(eq.a.G/benchmarkG)
  end


  save(filename,"collections",collections,"collections2gdp", collections2gdp,"allinfo",allinfo)
end


function laffer_taud_parallel(taudvec::FloatRange{Float64},filename::ASCIIString,benchmark::ASCIIString)

  Ntau=size(taudvec)[1]
  arguments = Array(EquilArgument,(Ntau,));

  pr,eq,tau,pa=load(benchmark, "pr","eq","tau","pa");

  for j=1:Ntau
    tauhat=Taxes(taudvec[j],tau.c,tau.i,tau.g,tau.l)
    arguments[j]= EquilArgument(tauhat,pa,eq.w)
  end

  bigvector=pmap(compute_equilibrium,arguments)

  benchmarkG=eq.a.G;
  collections = Array(Float64,(Ntau,));
  collections2gdp = Array(Float64,(Ntau,));
  allinfo = Array(AllInfo,(Ntau,));

  for j=1:Ntau
    pr,eq=bigvector[j]
    collections[j]= eq.a.G/benchmarkG
    collections2gdp[j]= eq.a.G/eq.a.output
    allinfo[j]= AllInfo(pr,eq)
    println(eq.a.G/benchmarkG)
  end


  save(filename,"collections",collections,"collections2gdp", collections2gdp,"allinfo",allinfo)
end
