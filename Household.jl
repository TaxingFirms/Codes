type HouseholdProblem
  c::Float64
  l::Float64
  a::Float64
  mu::Float64
end

function householdProblem!(hpr::HouseholdProblem, aprime::Float64, mgutilityprime::Float64, r::Float64, w::Float64, tau::Taxes, pa::Param)
  labor= (w/pa.H)^pa.psi;
  dsc = (1+ (1-tau.i)*r)*pa.beta;
  hpr.a= (mgutilityprime/dsc)^(-1/pa.sigma) - w*labor + aprime + (pa.H/(1+1/pa.psi))*labor^(1+pa.psi);
  hpr.c= assets +  w*labor - aprime/dsc;
  hpr.mu= (consumption - (pa.H/(1+1/pa.psi))*labor^(1+pa.psi))^(-pa.sigma);
  hpr.l=labor;
end
