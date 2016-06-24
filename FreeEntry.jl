function free_entry!(pr::FirmProblem, p::Equilibrium, tau:: Taxes, fp::FirmParam, hp::HouseholdParam; xtol = .00001)
  f(x) = expvalentry!(x,pr,p,tau,fp,hp);

  expvalentry= compute_expvalentry(pr,p,tau,fp,hp);
  println(" Expected Value Entrants = ",expvalentry, " w = ",p.w);

  #The folowing block is meant to speed up bisection a little.
  center= p.w;
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

function expvalentry!(w::Real,pr::FirmProblem, p::Equilibrium, tau:: Taxes, fp::FirmParam, hp::HouseholdParam)
  #Computes the interest rate consistent with free entry.
  p.w=w;
  println("w= ",w);
  pr  = init_firmproblem(p,tau,fp,hp); #When the prices change, optimal size changes and we need to update the grids.

  firmVFIParallelOmega!(pr,p,tau,fp);
  expvalentry=compute_expvalentry(pr,p,tau,fp,hp);
  println("w= ",w, " expvalentry = ", expvalentry);

  return expvalentry
end

function compute_expvalentry(pr::FirmProblem, p::Equilibrium, tau:: Taxes, fp::FirmParam, hp::HouseholdParam)
  #Computes the expected value just once, without updating prices

  expvalentry=0;
  for i_z = 1:pr.Nz
    expvalentry = pr.firmvaluegrid[1,i_z]*fp.invariant_distr[i_z] - fp.e;
  end
  expvalentry
end
