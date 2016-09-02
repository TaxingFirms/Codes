type Economy
  pr::FirmProblem
  eq::Equilibrium
  tau::Taxes
  cev::Float64
end

#taxesb=[0.33 0.31 0.29 0.27 0.25 0.23 0.21 0.19 0.17 0.15 0.13 0.11 0.09 0.07 0.05 0.03 0.01];
function runreforms(filename::ASCIIString, pr::FirmProblem, eq::Equilibrium, tau::Taxes, pa::Param)
  taxesb=[0.3 0.25 0.2 0.15 0.1 0.05 0.0];
  taxesa=[0.4 0.45 0.5];
  ~,Nrb=size(taxesb);
  ~,Nra=size(taxesa);



  ref=Array(Economy,(Nrb+Nra+1,));
  ref[1]=Economy(pr,eq,tau,0.0);

  for j=1:Nrb
    if isnan(ref[j].tau.d)
      break
    end
      rpr,req,rtau = taxreform2(taxesb[j], ref[j].eq, ref[j].tau, pa; tol=10.0^-2.0,update=0.7,momentsprint=true);
      cev= (req.a.consumption - eq.a.consumption)/eq.a.consumption - pa.H/(eq.a.consumption*(1+1/pa.psi))*( (req.w/pa.H)^(1+pa.psi) - (eq.w/pa.H)^(1+pa.psi) );
      println("cev = ",cev)
      ref[j+1]=Economy(rpr,req,rtau,cev);
      save(filename,"ref",ref,"pa",pa);
  end


  rpr,req,rtau = taxreform2(taxesa[1], ref[1].eq, ref[1].tau, pa; tol=10.0^-2.0,update=0.7,momentsprint=true);
  ref[Nra+1]=Economy(rpr,req,rtau);

  for j=2:Nra
      rpr,req,rtau = taxreform2(taxesa[j], ref[j-1].eq, ref[j-1].tau, pa; tol=10.0^-2.0,update=0.7,momentsprint=true);
      ref[Nra+j]=Economy(rpr,req,rtau);
  end

  save(filename,"ref",ref,"pa",pa);
end

function runreform3(filename::ASCIIString, pr::FirmProblem, eq::Equilibrium, tau::Taxes, pa::Param)
  taxesb=[0.3 0.25 0.2 0.15 0.1 0.05 0.0];
  taxesa=[0.4 0.45 0.5];
  ~,Nrb=size(taxesb);
  ~,Nra=size(taxesa);



  ref=Array(Economy,(Nrb+Nra+1,));
  ref[1]=Economy(pr,eq,tau,0.0);

  for j=1:Nrb
    if isnan(ref[j].tau.d)
      break
    end
      rpr,req,rtau = taxreform3(taxesb[j], ref[j].eq, ref[j].tau, pa; tol=10.0^-2.0,update=0.7,momentsprint=true);
      cev= (req.a.consumption - eq.a.consumption)/eq.a.consumption - pa.H/(eq.a.consumption*(1+1/pa.psi))*( (req.w/pa.H)^(1+pa.psi) - (eq.w/pa.H)^(1+pa.psi) );
      println("cev = ",cev)
      ref[j+1]=Economy(rpr,req,rtau,cev);
      save(filename,"ref",ref,"pa",pa);
  end


  rpr,req,rtau = taxreform3(taxesa[1], ref[1].eq, ref[1].tau, pa; tol=10.0^-2.0,update=0.7,momentsprint=true);
  ref[Nra+1]=Economy(rpr,req,rtau);

  for j=2:Nra
      rpr,req,rtau = taxreform2(taxesa[j], ref[j-1].eq, ref[j-1].tau, pa; tol=10.0^-2.0,update=0.7,momentsprint=true);
      ref[Nra+j]=Economy(rpr,req,rtau);
  end

  save(filename,"ref",ref,"pa",pa);
end
