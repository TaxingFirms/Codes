function labormarket!(w::Float64,E::Float64, exitrule::Array{Bool,3}, res::ResultsFP, pr::FirmProblem, p::Equilibrium, tau::Taxes, fp::FirmParam)
#This function canges the wage and updats the value an policy finctions in pr and res
p.w=w;
println("Wage = ", w, " = ", p.w)
firmVFIParallelOmega!(pr,p,tau,fp);
res=copy_opt_policies(pr);
getpolicies!(res,pr,p,tau,fp);  #r is updated exctracts policies

distr=stationarydist(E,res, pr, p, tau,fp);

capital, bonds, labor_d, gdp, corptax, inctax, netdistributions, liquidations, liquidationcosts = aggregates(E, distr, res, pr, p, tau, fp);

labor_s = 1/hp.H - 1/p.w*((1-tau.i)*p.r*bonds + netdistributions +liquidations);

return labor_d - labor_s


end


function clear_labormarket!(E::Float64, exitrule::Array{Bool,3}, res::ResultsFP, pr::FirmProblem, p::Equilibrium, tau::Taxes, fp::FirmParam; tol=10^-4.0)
# It updates the wages until labormarket return 0
excessdemand =labormarket!(p.w::Float64,E::Float64, exitrule::Array{Bool,3}, res::ResultsFP, pr::FirmProblem, p::Equilibrium, tau::Taxes, fp::FirmParam);
println("Excess of labor demand = ", excessdemand )
excessdemand>0?
  direction = +1:
  direction = -1;

f(w)=labormarket!(w::Float64,E::Float64, exitrule::Array{Bool,3}, res::ResultsFP, pr::FirmProblem, p::Equilibrium, tau::Taxes, fp::FirmParam);


#Look for a change of sign
center= p.w;
radius=10^-2.0;
flag = false;
newexcessdemand=excessdemand;
while !flag
  radius*=2.0;
  newcenter=center + direction*radius;
  println("New wage = ", newcenter)
  newexcessdemand = f(newcenter);
  println("Excess of labor demand = ", newexcessdemand, "Wage = ", newcenter)


  if excessdemand*newexcessdemand < 0
    flag=true;
  else
    center=newcenter;
    excessdemand=newexcessdemand;
  end
end

fzero(f,center,newcenter)
end


###############################
## Constant Exit Probability ##
###############################


function labormarketCE!(w::Float64,E::Float64, exitprob::Float64, res::ResultsFP, pr::FirmProblem, p::Equilibrium, tau::Taxes, fp::FirmParam)
  #This function canges the wage and updats the value an policy finctions in pr and res
   p.w = w;
  # println("Wage = ", w)
  firmVFIParallelOmegaCE!(exitprob,pr,p,tau,fp);
  res = copy_opt_policies(pr);
  getpolicies!(res,pr,p,tau,fp);  #r is updated exctracts policies

  distr = stationarydistCE(E, pr, p, tau, fp);

  capital, bonds, labor_d, gdp, corptax, inctax, netdistributions, liquidations, liquidationcosts = aggregatesCE(exitprob, E, distr, res, pr, p, tau, fp)

  hp = init_hhparameters();

  labor_s = 1/hp.H - 1/p.w*((1-tau.i)*p.r*bonds + netdistributions +liquidations);


  #Save variables
  grossdividends=sum(distr.*res.grossdividends)
  divtax= tau.d*grossdividends;
  consumption = p.w/hp.H;
  p.consumption=consumption;
  p.distr= distr;
  p.netdistributions=netdistributions;
  p.agginterests= (1-tau.i)*p.r*bonds;
  p.collections=Taxes(divtax,corptax,inctax,0);
  p.output= gdp;
  p.laborsupply=labor_s;
  p.welfare= log(consumption) - hp.H* labor_s;


  return labor_d - labor_s

end


# function clear_labormarketCE!(E::Float64, exitprob::Float64, res::ResultsFP, pr::FirmProblem, p::Equilibrium, tau::Taxes, fp::FirmParam; tol=10^-4.0)
# # It updates the wages until labormarket return 0
#   excessdemand =labormarketCE!(p.w,E, exitprob, res, pr, p, tau, fp);
#   println("Excess of labor demand = ", excessdemand )
#   excessdemand>0?
#     direction = +1:
#     direction = -1;

#   f(x)=labormarket!(x,E, exitprob, res, pr, p, tau, fp);


#   #Look for a change of sign
#   center= p.w;
#   radius=10^-2.0;
#   flag = false;
#   newexcessdemand=excessdemand;
#   while !flag
#     radius*=2.0;
#     newcenter=center + direction*radius;
#     println("New wage = ", newcenter)
#     newexcessdemand = f(newcenter);
#     println("Excess of labor demand = ", newexcessdemand, "Wage = ", newcenter)


#     if excessdemand*newexcessdemand < 0
#       flag=true;
#     else
#       center=newcenter;
#       excessdemand=newexcessdemand;
#     end
#   end

#   fzero(f,center,newcenter)
# end




function clear_labormarketCE!(E::Float64, exitprob::Float64, res::ResultsFP, pr::FirmProblem, p::Equilibrium, tau::Taxes, fp::FirmParam; tol=10^-4.0,maxit = 1000)

  topWage = 0.9
  minWage = 0.3
  currentWage = (topWage+minWage)/2.0

  excessdemand =labormarketCE!(currentWage,E, exitprob, res, pr, p, tau, fp);

  count = 1

  while abs(excessdemand) > tol || count > maxit

    println("Excess of labor demand = ", excessdemand, " at wage ", currentWage )
    if excessdemand < 0.0
      topWage = currentWage
    else
      minWage = currentWage
    end

    currentWage = (topWage+minWage)/2.0
    println("New wage ", currentWage)
    count += 1

    excessdemand =labormarketCE!(currentWage,E, exitprob, res, pr, p, tau, fp);

  end

  p.w = currentWage

end
