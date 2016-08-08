
taxesb=[0.3 0.25 0.2 0.15 0.1 0.05 0.0];
taxesa=[0.4 0.45 0.5];
~,Nrb=size(taxesb);
~,Nra=size(taxesa);

rexceptions=Array(Any,(Nrb+Nra+1,));

type Economy
  pr::FirmProblem
  eq::Equilibrium
  tau::Taxes
end

ref=Array(Economy,(Nrb+Nra+1,));
ref[1]=Economy(pr,eq,tau);

for j=1:Nrb
#  try
    rpr,req,rtau = taxreform2(taxesb[j], ref[j].eq, ref[j].tau, pa; tol=5.0^-4.0,update=0.7);
    ref[j+1]=Economy(rpr,req,rtau);
#  catch rexceptions[j]
#    println("Error computing reform ",j)
#  end
end

save("CounterfactualBoth.jld","ref",ref,"pa",pa);
