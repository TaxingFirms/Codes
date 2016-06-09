## Find zero
function myfzero(f, a::Float64, b::Float64 , fa::Float64, fb::Float64; xtol::Real=0.0, xtolrel::Real=0.0, verbose::Bool=false)

    if (xtol < 0.0) | (xtolrel < 0.0)
        error("Tolerances must be non-negative")
    end

    x0, y0 = a, fa
    x2, y2 = b, fb

    y0 == 0 && return a
    y2 == 0.0 && return b
    if sign(y0) * sign(y2) > 0
        error("[a,b] is not a bracket")
    end

    x1 = x0 - y0/(y2 - y0)*(x2 - x0)
    y1 = f(x1)


    while x0 < x1 && x1 < x2
        if sign(y0) == sign(y1)
            x0, x2 = x1, x2
            y0, y2 = y1, y2
        else
            x0, x2 = x0, x1
            y0, y2 = y0, y1
        end

        x1 = secant(f,x0, x2)
        y1 = f(x1)

        sign(y1) == 0 && return x1
        abs(x2 - x0) <= max(xtol, xtolrel*abs(x1)) && return(x1)

        verbose && println(@sprintf("xi =%18.15f,\t f(xi) = %18.15f", x1, f(x1)))
    end

    return abs(y0) < abs(y2) ? x0 : x2
end

# take a secant step, if the resulting guess is very close to a or b,
function secant(f::Function, a::Real, b::Real)
    c = a - f(a)/(f(b) - f(a))*(b - a)
    tol = eps(1.0)*5
    if c <= a + abs(a)*tol || c >= b - abs(b)*tol
        return a + (b - a)/2
    end
    return c
end



sameshock = zeros(Float64,(pr.Nomega,pr.Nz));

for i_z in 1:pr.Nz, i_omega in 1:pr.Nomega
  kprime=res.kprime[i_omega,i_z];
  qprime=res.qprime[i_omega,i_z];
  if !res.exitrule[i_omega,i_z,i_z] > 0
    omegaprime = omegaprimefun(kprime,qprime,i_z,p,tau,fp)
    i_omegaprime = closestindex(omegaprime, pr.omega.step);

    #The block below checks that the index is within reasonable bounds
    if i_omegaprime<1 || i_omegaprime>pr.Nomega
      if i_omega==pr.Nomega || i_omegaprime < (pr.Nomega + 3)
        i_omegaprime =pr.Nomega
      elseif i_omega==1 || i_omegaprime > -3
        i_omegaprime =1;
      else
        error("omega' out of the grid ", "i_z = ", i_z, "i_omega = ", i_omega, "i_z' = ",i_zprime, "i_omega' = ",i_omegaprime)
      end
    end
    sameshock[i_omega,i_z] = omegaprime;
  end
end


function aggregateValue(policy::Array{Float64,2},res::ResultsFP, pr::FirmProblem, p::Prices, tau::Taxes, fp::FirmParam)

  value = 0.0
  for i_omega in 1:pr.Nomega, i_z in 1:pr.Nz
    if distribuition[i_omega,i_z] > 0
      value += policy[i_omega,i_z]*distribuition*exitDecision
    end
  end
  value
end



type State
  i_omega::Int64
  omega::Real
  i_z::Int64
  p::Prices
  pr::FirmProblem
  tau::Taxes
  fp::FirmParam
end

function firmVFI!(statespace::Array{State,2}, tol=10.0^-3, maxit=500)
  dist=Inf;
  dif=similar(pr.firmvalueguess);
  it=1;
  while dist > tol
    println("it=", it);
    firmbellman!(statespace)
    dif = pr.firmvalueguess - pr.firmvaluegrid;
    dist= norm(dif, Inf);
    println("it=", it, "   dist=", dist);
    pr.firmvalueguess = deepcopy(pr.firmvaluegrid)
    pr.InterpolationGrid = map(x->CoordInterpGrid(pr.omega.grid,pr.firmvaluegrid[:,x],BCnearest, InterpLinear),1:pr.Nz)
    it>= maxit ? error("maximum number of iterations reached"):it+=1;
  end
end



function init_state(p::Prices, pr::FirmProblem, tau::Taxes, fp::FirmParam)
  N=length(pr.omega.grid)
  input= Array(State,(N,pr.Nz));
  for j= 1:pr.Nz
    for i = 1:N
      input[i,j]= State(i,pr.omega.grid[i],j , p, pr, tau, fp);
    end
  end
  input
end


function maximize_aux(state::State)
  v=maximizationstep(state.omega, state.i_z, state.p, state.pr, state.tau, state.fp,false);
  pr.firmvaluegrid[state.i_omega,state.i_z] = deepcopy(v);
end


# Firm Problem
function firmbellman!(input::Array{State,2})
  pmap(maximize_aux,input)
  #value= pmap(maximize_aux,input)
end


type ResultsFP
  firmvalue :: Array{Float64,2}
  interpolation::Array{CoordInterpGrid,1}
  kprime :: Array{Float64,2}
  qprime :: Array{Float64,2}
  capital :: Array{Float64,2}
  labor :: Array{Float64,2}
  output :: Array{Float64,2}
  debt :: Array{Float64,2}
  distributions :: Array{Float64,2}
  financialcosts :: Array{Float64,2}
  exitprobability :: Array{Float64,2}
  exitrule:: Array{Bool,3}
  expostexitrule :: Array{Float64,2}
  positivedistributions :: Array{Bool,2}
end

function copy_opt_policies(pr::FirmProblem)
  firmvalue= pr.firmvaluegrid;
  interpolation = pr.InterpolationGrid;
  kprime = pr.kpolicygrid;
  qprime = pr.qpolicygrid;

  distributions = Array(Float64,(pr.Nomega,pr.Nz));
  financialcosts = zeros(Float64,(pr.Nomega,pr.Nz));
  capital = Array(Float64,(pr.Nomega,pr.Nz));
  labor = Array(Float64,(pr.Nomega,pr.Nz));
  output = Array(Float64,(pr.Nomega,pr.Nz));
  debt = Array(Float64,(pr.Nomega,pr.Nz));
  exitprobability = Array(Float64,(pr.Nomega,pr.Nz));
  exitrule = falses(pr.Nomega,pr.Nz,pr.Nz);
  positivedistributions = falses(pr.Nomega,pr.Nz);

  ResultsFP(firmvalue, interpolation,kprime,qprime,capital, labor, output,debt, distributions, financialcosts ,exitprobability, exitrule, positivedistributions)
end

#Gets the exit rules, distributions and other quantities of interest

function getpolicies!(res::ResultsFP, pr::FirmProblem, p::Prices, tau::Taxes, fp::FirmParam)
  for (i_z,z) in enumerate(fp.zgrid)
    for (i_omega, omega) in enumerate(pr.omega.grid)
      kprime= res.kprime[i_omega, i_z];
      qprime= res.qprime[i_omega, i_z];
      if grossdistributions(omega,kprime,qprime,fp) >=0
        res.positivedistributions[i_omega, i_z]= true;
        res.distributions[i_omega, i_z]= (1-tau.d)*grossdistributions(omega,kprime,qprime,fp);
      else
        res.distributions[i_omega, i_z]= (1+fp.lambda1)*grossdistributions(omega,kprime,qprime,fp) - fp.lambda0;
        res.financialcosts[i_omega, i_z]= fp.lambda1*grossdistributions(omega,kprime,qprime,fp) -fp.lambda0;
      end

      prexit=0.0;
       #This just avoids computing this constant again and again while computing omegaprime
      exitvalue = fp.theta*(1-fp.delta)*kprime - (1+p.r)*qprime;
      for (i_zprime,zprime) in enumerate(fp.zgrid)
        lprime=(fp.alphal*zprime*(kprime^fp.alphak)/p.w)^(1/(1-fp.alphal));
        omegaprime = omegaprimefun(kprime,qprime,i_zprime,p,tau,fp);
        contvalue = firmvaluefunction(omegaprime,i_zprime,pr);

        i_omegaprime = closestindex(omegaprime, pr.omega.step);
        #The block below checks that the index is within reasonable bounds
        if i_omegaprime<1 || i_omegaprime>pr.Nomega
          if i_omega==pr.Nomega || i_omegaprime < (pr.Nomega + 3)
            i_omegaprime =pr.Nomega
          elseif i_omega==1 || i_omegaprime > -3
            i_omegaprime =1;
          else
            error("omega' out of the grid ", "i_z = ", i_z, "i_omega = ", i_omega, "i_z' = ",i_zprime, "i_omega' = ",i_omegaprime)
          end
        end

        res.labor[i_omegaprime, i_zprime] = lprime;
        res.capital[i_omegaprime, i_zprime] = kprime;
        res.output[i_omegaprime, i_zprime] = zprime*(kprime^fp.alphak)*(lprime^fp.alphal) -fp.f;
        res.debt[i_omegaprime, i_zprime] = qprime;

        if exitvalue>=contvalue #if indiferent, firms choose to exit
          res.exitrule[i_omega, i_z,i_zprime]=true;
          prexit+=fp.ztrans[i_zprime,i_z];
          res.expostexitrule[i_omegaprime, i_zprime] =true;
        end
      end
      res.exitprobability[i_omega, i_z]=prexit;
    end
  end
end

