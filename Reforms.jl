
function Reform1Vector(filename::ASCIIString, taxvector::Array, pr::FirmProblem, eq::Equilibrium, tau::Taxes, pa::Param)
    #Loops through several taxes and saves a file named "filename" with the results
    #It is recommended that first tax starts close to benchmark and moves monotonically
    ~,Nv=size(taxvector);

    ref=Array(Economy,(Nv+1,));
    ref[1]=Economy(pr,eq,tau,0.0);

    for j=1:Nv
        #initialguess = copy(ref[j].pr.firmvaluegrid)
        rpr,req,rtau = taxreform1(taxvector[j], ref[1].eq.a.G, ref[j].eq, ref[j].tau, pa; tol=5.0*10.0^-3.0,update=0.98,momentsprint=false);
        cev= (req.a.consumption - eq.a.consumption)/eq.a.consumption - pa.H/(eq.a.consumption*(1+pa.psi))*( (req.w*(1-rtau.l)/pa.H)^(1+pa.psi) - (eq.w*(1-tau.l)/pa.H)^(1+pa.psi) );
        ref[j+1]=Economy(rpr,req,rtau,cev);
        save(filename,"ref",ref,"pa",pa);
    end
    ref
end



function Reform2Vector(filename::ASCIIString, taxvector::Array, pr::FirmProblem, eq::Equilibrium, tau::Taxes, pa::Param)
    #Loops through several taxes and saves a file named "filename" with the results
    #It is recommended that first tax starts close to benchmark and moves monotonically
    ~,Nv=size(taxvector);

    ref=Array(Economy,(Nv+1,));
    ref[1]=Economy(pr,eq,tau,0.0);

    for j=1:Nv
        #initialguess = copy(ref[j].pr.firmvaluegrid)
        rpr,req,rtau = taxreform2(taxvector[j],ref[1].eq.a.G, ref[j].eq, ref[j].tau, pa; tol=5.0*10.0^-3.0,update=0.98,momentsprint=false);
        cev= (req.a.consumption - eq.a.consumption)/eq.a.consumption - pa.H/(eq.a.consumption*(1+pa.psi))*( (req.w*(1-rtau.l)/pa.H)^(1+pa.psi) - (eq.w*(1-tau.l)/pa.H)^(1+pa.psi) );
        ref[j+1]=Economy(rpr,req,rtau,cev);
        save(filename,"ref",ref,"pa",pa);
    end
    ref
end


function Reform2TauIVector(filename::ASCIIString, taxvector::Array, pr::FirmProblem, eq::Equilibrium, tau::Taxes, pa::Param)
    #Loops through several taxes and saves a file named "filename" with the results
    #It is recommended that first tax starts close to benchmark and moves monotonically
    ~,Nv=size(taxvector);

    ref=Array(Economy,(Nv+1,));
    ref[1]=Economy(pr,eq,tau,0.0);

    for j=1:Nv
        #initialguess = copy(ref[j].pr.firmvaluegrid)
        rpr,req,rtau = taxreform2_taui(taxvector[j],ref[1].eq.a.G, ref[j].eq, ref[j].tau, pa; tol=5.0*10.0^-3.0,update=0.98,momentsprint=false);
        cev= (req.a.consumption - eq.a.consumption)/eq.a.consumption - pa.H/(eq.a.consumption*(1+pa.psi))*( (req.w*(1-rtau.l)/pa.H)^(1+pa.psi) - (eq.w*(1-tau.l)/pa.H)^(1+pa.psi) );
        ref[j+1]=Economy(rpr,req,rtau,cev);
        save(filename,"ref",ref,"pa",pa);
    end
    ref
end


function Reform3Vector(filename::ASCIIString, taxvector::Array, pr::FirmProblem, eq::Equilibrium, tau::Taxes, pa::Param)
    #Loops through several taxes and saves a file named "filename" with the results
    #It is recommended that first tax starts close to benchmark and moves monotonically
    ~,Nv=size(taxvector);

    ref=Array(Economy,(Nv+1,));
    ref[1]=Economy(pr,eq,tau,0.0);

    for j=1:Nv
        #initialguess = copy(ref[j].pr.firmvaluegrid)
        rpr,req,rtau = taxreform3(taxvector[j], ref[1].eq.a.G, ref[j].eq, ref[j].tau, pa; tol=5.0*10.0^-3.0,update=0.98,momentsprint=false);
        cev= (req.a.consumption - eq.a.consumption)/eq.a.consumption - pa.H/(eq.a.consumption*(1+pa.psi))*( (req.w*(1-rtau.l)/pa.H)^(1+pa.psi) - (eq.w*(1-tau.l)/pa.H)^(1+pa.psi) );
        ref[j+1]=Economy(rpr,req,rtau,cev);
        save(filename,"ref",ref,"pa",pa);
    end
    ref
end



function Reform4Vector(filename::ASCIIString, vector::Array, pr::FirmProblem, eq::Equilibrium, tau::Taxes, pa::Param)
    #Loops through several taxes and saves a file named "filename" with the results
    #It is recommended that first tax starts close to benchmark and moves monotonically
    ~,Nv=size(vector);

    ref=Array(Economy,(Nv+1,));
    ref[1]=Economy(pr,eq,tau,0.0);

    for j=1:Nv
        #initialguess = copy(ref[j].pr.firmvaluegrid)
        rpr,req,rtau = taxreform4(taxesb[j], ref[j].eq, ref[j].tau, pa; tol=10.0^-3.0,update=0.9,momentsprint=true);
        cev= (req.a.consumption - eq.a.consumption)/eq.a.consumption - pa.H/(eq.a.consumption*(1+pa.psi))*( (req.w*(1-rtau.l)/pa.H)^(1+pa.psi) - (eq.w*(1-tau.l)/pa.H)^(1+pa.psi) );
        ref[j+1]=Economy(rpr,req,rtau,cev);
        save(filename,"ref",ref,"pa",pa);
    end
    ref
end


function Reform5Vector(filename::ASCIIString, vector::Array, pr::FirmProblem, eq::Equilibrium, tau::Taxes, pa::Param)
    #Loops through several taxes and saves a file named "filename" with the results
    #It is recommended that first tax starts close to benchmark and moves monotonically
    ~,Nv=size(vector);

    ref=Array(Economy,(Nv+1,));
    ref[1]=Economy(pr,eq,tau,0.0);

    for j=1:Nv
        #initialguess = copy(ref[j].pr.firmvaluegrid)
        rpr,req,rtau = taxreform5(taxesb[j], ref[j].eq, ref[j].tau, pa; tol=10.0^-3.0,update=0.7,momentsprint=true);
        cev= (req.a.consumption - eq.a.consumption)/eq.a.consumption - pa.H/(eq.a.consumption*(1+pa.psi))*( (req.w*(1-rtau.l)/pa.H)^(1+pa.psi) - (eq.w*(1-tau.l)/pa.H)^(1+pa.psi) );
        ref[j+1]=Economy(rpr,req,rtau,cev);
        save(filename,"ref",ref,"pa",pa);
    end
    ref
end
