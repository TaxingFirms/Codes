
function Reform1Vector(filename::ASCIIString, taxvector::Array, pr::FirmProblem, eq::Equilibrium, tau::Taxes, pa::Param; bctol::Float64=5.0*10.0^-3.0, update::Float64=0.9)
    #Loops through several taxes and saves a file named "filename" with the results
    #It is recommended that first tax starts close to benchmark and moves monotonically
    ~,Nv=size(taxvector);

    ref=Array(Economy,(Nv+1,));
    ref[1]=Economy(pr,eq,tau,0.0);

    for j=1:Nv
        #initialguess = copy(ref[j].pr.firmvaluegrid)
        rpr,req,rtau = taxreform1(taxvector[j], ref[1].eq.a.G, ref[j].pr, ref[j].eq, ref[j].tau, pa; tol=bctol,update=update,momentsprint=false);
        cev= (req.a.consumption - eq.a.consumption)/eq.a.consumption - pa.H/(eq.a.consumption*(1+pa.psi))*( (req.w*(1-rtau.l)/pa.H)^(1+pa.psi) - (eq.w*(1-tau.l)/pa.H)^(1+pa.psi) );
        ref[j+1]=Economy(rpr,req,rtau,cev);
        save(filename,"ref",ref,"pa",pa);
    end
    ref
end



function Reform2Vector(filename::ASCIIString, taxvector::Array, pr::FirmProblem, eq::Equilibrium, tau::Taxes, pa::Param; bctol::Float64=5.0*10.0^-3.0, update::Float64=0.9)
    #Loops through several taxes and saves a file named "filename" with the results
    #It is recommended that first tax starts close to benchmark and moves monotonically
    ~,Nv=size(taxvector);

    ref=Array(Economy,(Nv+1,));
    ref[1]=Economy(pr,eq,tau,0.0);

    for j=1:Nv
        #initialguess = copy(ref[j].pr.firmvaluegrid)
        rpr,req,rtau = taxreform2(taxvector[j],ref[1].eq.a.G, ref[j].pr, ref[j].eq, ref[j].tau, pa; tol=bctol,update=update,momentsprint=false);
        cev= (req.a.consumption - eq.a.consumption)/eq.a.consumption - pa.H/(eq.a.consumption*(1+pa.psi))*( (req.w*(1-rtau.l)/pa.H)^(1+pa.psi) - (eq.w*(1-tau.l)/pa.H)^(1+pa.psi) );
        ref[j+1]=Economy(rpr,req,rtau,cev);
        save(filename,"ref",ref,"pa",pa);
    end
    ref
end


function Reform2TauIVector(filename::ASCIIString, taxvector::Array, pr::FirmProblem, eq::Equilibrium, tau::Taxes, pa::Param; bctol::Float64=5.0*10.0^-3.0, update::Float64=0.9)
    #Loops through several taxes and saves a file named "filename" with the results
    #It is recommended that first tax starts close to benchmark and moves monotonically
    ~,Nv=size(taxvector);

    ref=Array(Economy,(Nv+1,));
    ref[1]=Economy(pr,eq,tau,0.0);

    for j=1:Nv
        #initialguess = copy(ref[j].pr.firmvaluegrid)
        rpr,req,rtau = taxreform2_taui(taxvector[j],ref[1].eq.a.G, ref[j].pr, ref[j].eq, ref[j].tau, pa; tol=bctol,update=update,momentsprint=false);
        cev= (req.a.consumption - eq.a.consumption)/eq.a.consumption - pa.H/(eq.a.consumption*(1+pa.psi))*( (req.w*(1-rtau.l)/pa.H)^(1+pa.psi) - (eq.w*(1-tau.l)/pa.H)^(1+pa.psi) );
        ref[j+1]=Economy(rpr,req,rtau,cev);
        save(filename,"ref",ref,"pa",pa);
    end
    ref
end


function Reform3Vector(filename::ASCIIString, taxvector::Array, pr::FirmProblem, eq::Equilibrium, tau::Taxes, pa::Param; bctol::Float64=5.0*10.0^-3.0, update::Float64=0.9)
    #Loops through several taxes and saves a file named "filename" with the results
    #It is recommended that first tax starts close to benchmark and moves monotonically
    ~,Nv=size(taxvector);

    ref=Array(Economy,(Nv+1,));
    ref[1]=Economy(pr,eq,tau,0.0);

    for j=1:Nv
        #initialguess = copy(ref[j].pr.firmvaluegrid)
        rpr,req,rtau = taxreform3(taxvector[j], ref[1].eq.a.G, ref[j].pr, ref[j].eq, ref[j].tau, pa; tol=bctol,update=update,momentsprint=false);
        cev= (req.a.consumption - eq.a.consumption)/eq.a.consumption - pa.H/(eq.a.consumption*(1+pa.psi))*( (req.w*(1-rtau.l)/pa.H)^(1+pa.psi) - (eq.w*(1-tau.l)/pa.H)^(1+pa.psi) );
        ref[j+1]=Economy(rpr,req,rtau,cev);
        save(filename,"ref",ref,"pa",pa);
    end
    ref
end


function Reform4Vector(filename::ASCIIString, taxvector::Array, pr::FirmProblem, eq::Equilibrium, tau::Taxes, pa::Param; bctol::Float64=5.0*10.0^-3.0, update::Float64=0.9)
    #Loops through several taxes and saves a file named "filename" with the results
    #It is recommended that first tax starts close to benchmark and moves monotonically
    ~,Nv=size(taxvector);

    ref=Array(Economy,(Nv+1,));
    ref[1]=Economy(pr,eq,tau,0.0);

    for j=1:Nv
        #initialguess = copy(ref[j].pr.firmvaluegrid)
        rpr,req,rtau = taxreform4(taxvector[j], ref[1].eq.a.G, ref[j].pr, ref[j].eq, ref[j].tau, pa; tol=bctol,update=update,momentsprint=false);
        cev= (req.a.consumption - eq.a.consumption)/eq.a.consumption - pa.H/(eq.a.consumption*(1+pa.psi))*( (req.w*(1-rtau.l)/pa.H)^(1+pa.psi) - (eq.w*(1-tau.l)/pa.H)^(1+pa.psi) );
        ref[j+1]=Economy(rpr,req,rtau,cev);
        save(filename,"ref",ref,"pa",pa);
    end
    ref
end



function Reform5Vector(filename::ASCIIString, taxvector::Array, pr::FirmProblem, eq::Equilibrium, tau::Taxes, pa::Param; bctol::Float64=5.0*10.0^-3.0, update::Float64=0.9)
    #Loops through several taxes and saves a file named "filename" with the results
    #It is recommended that first tax starts close to benchmark and moves monotonically
    ~,Nv=size(taxvector);

    ref=Array(Economy,(Nv+1,));
    ref[1]=Economy(pr,eq,tau,0.0);
    println(@sprintf(" Dividend Base \t LabTax Collts \t  Consumption \t Mass Entrants \t     Wage    \t    Welfare  \t    TFP   "))
    println(@sprintf("%9.4f \t %9.4f \t %9.4f \t %9.4f \t %9.4f  \t %9.4f \t %9.4f ", eq.a.collections.d/rtau.d, eq.a.collections.l, eq.a.consumption, eq.E, eq.w, eq.a.welfare , eq.a.output/(eq.a.capital^pa.alphak*eq.a.laborsupply^pa.alphal)))


    for j=1:Nv
        #initialguess = copy(ref[j].pr.firmvaluegrid)
        rpr,req,rtau = taxreform5(taxvector[j], ref[1].eq.a.G, ref[j].pr, ref[j].eq, ref[j].tau, pa; tol=bctol,update=update,momentsprint=false);
        cev= (req.a.consumption - eq.a.consumption)/eq.a.consumption - pa.H/(eq.a.consumption*(1+pa.psi))*( (req.w*(1-rtau.l)/pa.H)^(1+pa.psi) - (eq.w*(1-tau.l)/pa.H)^(1+pa.psi) );
        ref[j+1]=Economy(rpr,req,rtau,cev);
        save(filename,"ref",ref,"pa",pa);

        println(@sprintf(" Dividend Base \t LabTax Collts \t  Consumption \t Mass Entrants \t     Wage    \t    Welfare  \t    TFP   "))
        println(@sprintf("%9.4f \t %9.4f \t %9.4f \t %9.4f \t %9.4f  \t %9.4f \t %9.4f ", req.a.collections.d/rtau.d, req.a.collections.l, req.a.consumption, req.E, req.w, req.a.welfare , req.a.output/(req.a.capital^pa.alphak*req.a.laborsupply^pa.alphal)))
    end
    ref
end



function Reform6Vector(filename::ASCIIString, taxvector::Array, pr::FirmProblem, eq::Equilibrium, tau::Taxes, pa::Param; bctol::Float64=5.0*10.0^-3.0, update::Float64=0.9)
    #Loops through several taxes and saves a file named "filename" with the results
    #It is recommended that first tax starts close to benchmark and moves monotonically
    ~,Nv=size(taxvector);

    ref=Array(Economy,(Nv+1,));
    ref[1]=Economy(pr,eq,tau,0.0);

    for j=1:Nv
        #initialguess = copy(ref[j].pr.firmvaluegrid)
        rpr,req,rtau = taxreform6(taxvector[j], ref[1].eq.a.G, ref[j].pr, ref[j].eq, ref[j].tau, pa; tol=bctol,update=update,momentsprint=true);
        cev= (req.a.consumption - eq.a.consumption)/eq.a.consumption - pa.H/(eq.a.consumption*(1+pa.psi))*( (req.w*(1-rtau.l)/pa.H)^(1+pa.psi) - (eq.w*(1-tau.l)/pa.H)^(1+pa.psi) );
        ref[j+1]=Economy(rpr,req,rtau,cev);
        save(filename,"ref",ref,"pa",pa);
    end
    ref
end

function ReformVector(taxreform::Function,filename::ASCIIString, taxvector::Array, pr::FirmProblem, eq::Equilibrium, tau::Taxes, pa::Param; bctol::Float64=5.0*10.0^-3.0, update::Float64=0.9)
  ~,Nv=size(taxvector);

  ref=Array(Economy,(Nv+1,));
  ref[1]=Economy(pr,eq,tau,0.0);

  for j=1:Nv
      #initialguess = copy(ref[j].pr.firmvaluegrid)
      rpr,req,rtau = taxreform(taxvector[j], ref[1].eq.a.G, ref[j].pr, ref[j].eq, ref[j].tau, pa; tol=bctol,update=update,momentsprint=true);
      cev= (req.a.consumption - eq.a.consumption)/eq.a.consumption - pa.H/(eq.a.consumption*(1+pa.psi))*( (req.w*(1-rtau.l)/pa.H)^(1+pa.psi) - (eq.w*(1-tau.l)/pa.H)^(1+pa.psi) );
      ref[j+1]=Economy(rpr,req,rtau,cev);
      save(filename,"ref",ref,"pa",pa);
  end
  ref
end
