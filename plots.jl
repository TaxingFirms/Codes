#@everywhere

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


pr,tau,fp,res,p=load("/Users/danielwillsr/Dropbox/1FirmTaxation/SimpleDiscreteTime/ModelResults.jld", "pr","tau","fp","res","p");
#pr,tau,fp,res,p=load("/Users/gcamilo/Dropbox/research/1FirmTaxation/SimpleDiscreteTime/ModelResults.jld", "pr","tau","fp","res","p");

omegagrid=pa.omega.grid;

regions=zeros(size(res.distributions));

for i_omega in 1:fp.Nomega, i_z in 1:fp.Nz
  nd = (pr.omega.grid[i_omega] - res.kprime[i_omega,i_z] + res.qprime[i_omega,i_z])
  if nd > pr.omega.step #Positive dividends
    regions[i_omega,i_z]=1;
  elseif nd < -pr.omega.step  #Negative dividends
   regions[i_omega,i_z]=-1;
  end
end



figure()
d= plot(pa.omega.grid, pr.distributions[:,1] )
k= plot(pa.omega.grid, pr.kpolicy[:,1] )
q= plot(pa.omega.grid, pr.qpolicy[:,1] )
  xlabel("Net worth")
  title("Policy functions (at z=7)")
  legend("dkq", loc="best")



figure()
subplot(121)
imshow(res.exitprobability,aspect="auto",extent=(fp.zgrid[1],fp.zgrid[end],pr.omega.ub,0))
    title("Exit Probability")
    xlabel("Productivity")
    ylabel("Net worth")

subplot(122)
imshow(regions,aspect="auto",extent=(fp.zgrid[1],fp.zgrid[end],pr.omega.ub,0))
    title("Distribution Regions")
    xlabel("Productivity")
    ylabel("Net worth")


#######
# Distribution of firms by net worth
#######

# Simple
oneDimension = squeeze(sum(p.distr,2),2)
figure()
plot(omegagrid,p.distr[:,2],"r--",linewidth=2.0,label="z low")
xlim(0,2)
plot(omegagrid,p.distr[:,5],"g+-",linewidth=2.0,label="z mid")
plot(omegagrid,p.distr[:,8],"b",linewidth=2.0,label="z high")
legend()
xlabel("Firm Net Worth")
ylabel("Firm Mass")


# Different Taxes
# prLow,tauLow,fpLow,resLow,pLow=load("/Users/danielwillsr/Dropbox/1FirmTaxation/SimpleDiscreteTime/Reform025.jld", "pr","tau","fp","res","p");
prLow,tauLow,fpLow,resLow,pLow=load("/Users/gcamilo/Dropbox/research/1FirmTaxation/SimpleDiscreteTime/Reform025.jld", "pr","tau","fp","res","p");
# prHigh,tauHigh,fpHigh,resHigh,pHigh=load("/Users/danielwillsr/Dropbox/1FirmTaxation/SimpleDiscreteTime/Reform045.jld", "pr","tau","fp","res","p");
prHigh,tauHigh,fpHigh,resHigh,pHigh=load("/Users/gcamilo/Dropbox/research/1FirmTaxation/SimpleDiscreteTime/Reform045.jld", "pr","tau","fp","res","p");
oneDimensionLow = squeeze(sum(pLow.distr,2),2)
oneDimensionHigh = squeeze(sum(pHigh.distr,2),2)

# PDF

figure()
plot(collect(prLow.omega.grid),oneDimensionLow/sum(oneDimensionLow),"g+-",linewidth=2.0,label="Taxes = 25%")
plot(omegagrid,oneDimension/sum(oneDimension),"r--",linewidth=2.0,label="Taxes = 35%")
xlim(0,2)
plot(collect(prHigh.omega.grid),oneDimensionHigh/sum(oneDimensionHigh),"b",linewidth=2.0,label="Taxes = 45%")
legend()
xlabel("Firm Net Worth")
ylabel("Firm Mass")

# CDF

figure()
plot(collect(prLow.omega.grid),cumsum(oneDimensionLow/sum(oneDimensionLow)),"g+-",linewidth=2.0,label="Taxes = 25%")
plot(omegagrid,cumsum(oneDimension/sum(oneDimension)),"r--",linewidth=2.0,label="Taxes = 35%")
xlim(0,2)
plot(collect(prHigh.omega.grid),cumsum(oneDimensionHigh/sum(oneDimensionHigh)),"b",linewidth=2.0,label="Taxes = 45%")
legend()
xlabel("Firm Net Worth")
ylabel("Firm Mass")
