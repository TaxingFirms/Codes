include("running.jl")
pr,eq,tau,pa=load("ModelResults.jld", "pr","eq","tau","pa");

include("TaxMaximization.jl")





govexp=eq.a.G;


Ntau=3;
dimtau=2;

#1.1 Create tax vectors: first component is taui, second is tauc
s=SobolSeq(dimtau); skip(s,Ntau);
sobolspace = Array(Float64,(Ntau,dimtau));
for j=1:Ntau
  sobolspace[j,:]= next(s);
end
#1.2 Create algorithm to move through the space continuously
taxvec = Array(Taxes,(Ntau,));
for i=1:Ntau
  taxvec[i]= Taxes(sobolspace[1],0.0,sobolspace[2],tau.g,tau.l)
end
#2. Initialize vector of objects (tau, eq, pr, pa) where eq,pr,pa are initialized but empty
welfarevec = Array(Float64,(Ntau,));


for i=1:Ntau
  @time welfarevec[i] =close_gov_tauc!(govexp,taxvec[i], pa; update= 0.9, verbose = true, tol =10.0^-3.0,  wguess = eq.w, outsideparallel=false)
end
