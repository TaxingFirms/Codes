

type Economy
  pr::FirmProblem
  eq::Equilibrium
  tau::Taxes
  cev::Float64
end


function Reform2Vector(filename::ASCIIString, vector::Array, pr::FirmProblem, eq::Equilibrium, tau::Taxes, pa::Param)
    #Loops through several taxes and saves a file named "filename" with the results
    #It is recommended that first tax starts close to benchmark and moves monotonically
    ~,Nv=size(vector);

    ref=Array(Economy,(Nv+1,));
    ref[1]=Economy(pr,eq,tau,0.0);

    for j=1:Nv
        initialguess = copy(ref[j].pr.firmvaluegrid)
        wguess = ref[j].eq.w
        rpr,req,rtau = taxreform2(taxesb[j], ref[j].eq, ref[j].tau, pa; tol=10.0^-2.0,update=0.7,momentsprint=true,firmvalueguess=initialguess);
        cev= (req.a.consumption - eq.a.consumption)/eq.a.consumption - pa.H/(eq.a.consumption*(1+pa.psi))*( (req.w*(1-rtau.l)/pa.H)^(1+pa.psi) - (eq.w*(1-tau.l)/pa.H)^(1+pa.psi) );
        ref[j+1]=Economy(rpr,req,rtau,cev);
        save(filename,"ref",ref,"pa",pa);
    end
    ref
end

