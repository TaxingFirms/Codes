# 0. PARAMETERS
# 0.1 Household Parameters -Conesa, Kitao, Krueger-
ssigma = 1.0;       # Risk aversion utility parameter
bbeta = 0.96;    # Discount factor
#The borrowing limit is implicitely defined by the lower bound of the grid.

# 0.2 Firm Parameters -Li, Whited, Wu-
alpha_k = 0.31; # Production function F = zK^alphaK * L^alphaL
alpha_l = 0.65; # Labor share
AA=1; #Scale productivity
aalpha= alpha_k/(1-alpha_l)


ddelta = 0.08;  # Capital depreciation
ttheta = 0.42   # Collateral requirement
rho_z = 0.767;
var_z = 0.211^2;


# 0.3 shock process
Nz= 10;        # number of z points
(logShocks,zTrans) = tauchen(Nz,rho_z,var_z^0.5); # Process of firm productivity z
shocks= exp(logShocks);
zGrid = AA*exp(shocks);

#Fixed costs are roughly such that about that firms on the lowest third of productivies are not solvent with an interest rate of 0.02
#zzff=zGrid[round(Nz/3)];
#kuncff= ( aalpha*zGrid[end]/(r+ddelta) )^(1/(1-aalpha));
ff= 3*10^6;

# 0.4 Precition parameters
ttolhh= 10^-3.0; #Tolerance for household policy iteration
maxithh=100;   # Max iterations for household policy iteration
ttolfirm= 10^-3.0; #Tolerance for firm value iteration
maxitfirm=100;   # Max iterations for firm value iteration

NbPrime = 2;  #Number of gridpoints for b'
NkPrime = 100;  #Number of gridpoints for k'
# b' is in [theta(1-delta)k' ; (1+r)(k'-w)]. When  the bounds of this interval are very close, I may be that NkPrime gridpoints are too much.
# I save some gridpoint by setting a floor of bPrimeStep for the seteps in the grid.
bPrimeStep = 10^-2.0; #Number of gridpoints for k'
kPrimeStep = 10^-2.0; #Number of gridpoints for k'

