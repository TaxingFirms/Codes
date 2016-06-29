# This script has the functions to compute a wage consistent with free entry
# It calls the value function iteration until the value function in
# pr::FirmProblem is zero for entrants. It saves the wage in eq.w

function free_entry!(eq::Equilibrium, pr::FirmProblem, tau:: Taxes, pa::Param, VFIfunction::Function; xtol = .00001)
  f(x) = expvalentry!(x,pr,eq,tau,pa,VFIfunction);

  expvalentry= compute_expvalentry(pr,eq,tau,pa);
  println(" Expected Value Entrants = ",expvalentry, " w = ",eq.w);

  #The folowing block is meant to speed up bisection a little.
  center= eq.w;
  radius=10^-2.0;
  flag = false;
  while !flag
    radius*=2.0;
    expvalentry>0?
      newcenter = center + radius:
      newcenter = center - radius;

    global newvalentry
    newvalentry = f(newcenter);

    if newvalentry*expvalentry < 0
      flag=true;
    else
      center=newcenter;
      expvalentry= newvalentry;
    end
    println("center = ", center, " radius = ", radius, " flag is ", flag)
  end

  expvalentry>0?
    fzero(f,center,center+radius; xtol=10^-3.0 ):
    fzero(f,max(center-radius,0); xtol=10^-3.0 )
end

function expvalentry!(w::Real,pr::FirmProblem, eq::Equilibrium, tau:: Taxes, pa::Param, VFIfunction::Function)
  #Computes the interest rate consistent with free entry.
  eq.w=w;
  println("w= ",w);

  VFIfunction(pr,eq,tau,pa);
  expvalentry=compute_expvalentry(pr,eq,tau,pa);
  println("w= ",w, " expvalentry = ", expvalentry);

  return expvalentry
end

function compute_expvalentry(pr::FirmProblem, eq::Equilibrium, tau:: Taxes, pa::Param)
  #Computes the expected value just once, without updating prices

  expvalentry=0;
  for i_z = 1:pa.Nz
    expvalentry = pr.firmvaluegrid[1,i_z]*pa.invariant_distr[i_z] - pa.e;
  end
  expvalentry
end
