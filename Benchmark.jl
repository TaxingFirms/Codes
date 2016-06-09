using QuantEcon:tauchen
using Grid:CoordInterpGrid, BCnan, BCnearest, InterpLinear
using Roots:fzero
using JLD
using PyPlot

include("Main.jl")
include("Firms.jl")
include("Policies.jl")
include("Distribution.jl")
include("FreeEntry.jl")
include("Aggregation.jl")
include("TaxReforms.jl")


pr,tau,fp,res,p=load("/Users/danielwillsr/Dropbox/1FirmTaxation/SimpleDiscreteTime/ModelResults.jld", "pr", "tau", "fp", "res","p")
hp = init_hhparameters();
fp  = init_firmparameters(hp);


omegagrid=pr.omega.grid
omegaprime =  sameshocktrans(pr,tau,fp)

plot(omegagrid, res.distributions[:,7])
plot(omegagrid, res.kprime[:,7])
plot(omegagrid, res.qprime[:,7])
plot(omegagrid, omegaprime[:,7])

ind= closestindex(omegaprime[end,7],pr.omega.step)
omegagrid[ind]

function sameshocktrans(pr,tau,fp)
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

  return sameshock
end
