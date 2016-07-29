# This script has the functions to compute a wage consistent with free entry
# It calls the value function iteration until the value function in
# pr::FirmProblem is zero for entrants. It saves the wage in eq.w

function free_entry!(eq::Equilibrium, pr::FirmProblem, tau:: Taxes, pa::Param,  VFIfunction::Function, maxroutine::Function; tol = 10^-7.0, verbose=true)
  f(x) = expvalentry!(x,pr,eq,tau,pa,VFIfunction,maxroutine,verbose);

  expvalentry= compute_expvalentry(pr,pa);
  verbose && println(" Expected Value Entrants = ",expvalentry, " w = ",eq.w);

  #The folowing block is meant to speed up bisection a little.
  center= eq.w;
  radius=10.0^-2.0;
  flag = false;
  while !flag
    radius*=1.5;
    expvalentry>0?
      newcenter = center + radius:
      newcenter = max(center - radius,eps());

    global newvalentry
    newvalentry = f(newcenter);

    if newvalentry*expvalentry < 0
      flag=true;
    else
      center=newcenter;
      expvalentry= newvalentry;
    end
    verbose && println(" flag is ", flag)
  end

  expvalentry>0?
    myfzero(f,center, expvalentry, center+radius, newvalentry; xtol= tol ):
    myfzero(f, max(center - radius,eps()), newvalentry, center, expvalentry; xtol= tol );
end

function expvalentry!(w::Real,pr::FirmProblem, eq::Equilibrium, tau:: Taxes, pa::Param, VFIfunction::Function, maxroutine::Function, verbose::Bool )
  #Computes the interest rate consistent with free entry.
  eq.w=w;

  VFIfunction(pr,eq,tau,pa; maximizationroutine=maxroutine, verbose = verbose);
  expvalentry=compute_expvalentry(pr,pa);
  verbose && println("w= ",w, " expvalentry = ", expvalentry);

  return expvalentry
end

function compute_expvalentry(pr::FirmProblem, pa::Param)
  #Computes the expected value just once, without updating prices

  expvalentry=0;
  for i_z = 1:pa.Nz
    omega0= omegaprimefun(pa.k0, 0.0, i_z, eq, tau, pa);
    omega0_ind = closestindex(omega0, pa.omega.step);
    expvalentry += pr.firmvaluegrid[1,i_z]*pa.invariant_distr[i_z] - ( pa.e  +  pa.k0);
  end
  expvalentry
end


function _middle(x::Float64, y::Float64)
  # Use the usual float rules for combining non-finite numbers
  if !isfinite(x) || !isfinite(y)
    return x + y
  end

  # Always return 0.0 when inputs have opposite sign
  if sign(x) != sign(y) && x != 0.0 && y != 0.0
    return 0.0
  end

  negate = x < 0.0 || y < 0.0

  xint = reinterpret(UInt64, abs(x))
  yint = reinterpret(UInt64, abs(y))
  unsigned = reinterpret(Float64, (xint + yint) >> 1)

  return negate ? -unsigned : unsigned
end


function myfzero(f, x0::Float64, y0::Float64, x2::Float64, y2::Float64; xtol::Real=0.0, xtolrel::Real=0.0, verbose::Bool=false)

    if (xtol < 0.0) | (xtolrel < 0.0)
        error("Tolerances must be non-negative")
    end

    #x0, y0 = a, f(a)
    #x2, y2 = b, f(b)

    if sign(y0) * sign(y2) > 0
        error("[a,b] is not a bracket")
    end

    x1 = _middle(x0, x2)
    y1 = f(x1)


    while x0 < x1 && x1 < x2
        if sign(y0) == sign(y1)
            x0, x2 = x1, x2
            y0, y2 = y1, y2
        else
            x0, x2 = x0, x1
            y0, y2 = y0, y1
        end

        x1 = _middle(x0, x2)
        y1 = f(x1)

        sign(y1) == 0 && return x1
        abs(x2 - x0) <= max(xtol, xtolrel*abs(x1)) && return(x1)

        verbose && println(@sprintf("xi =%18.15f,\t f(xi) = %18.15f", x1, f(x1)))
    end

    return abs(y0) < abs(y2) ? x0 : x2
end
