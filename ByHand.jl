
pr.exitprobability #

capital, debt, networth, dividends, investment, z_history_ind = simulation(50, 50,pr,pa; seed=1234);
figure()
plot(capital)

figure()
plot(dividends)

figure()
plot(debt)


# Calibration Board presentation
pa  = init_parameters( H=1.3, ff= 0.015, llambda0=0.02, llambda1= 0.04, ddelta = 0.12, ttheta = 0.42,rhoz= 0.76, ssigmaz= 0.0352);
tau = init_taxes(ttaud =0.12, ttauc= 0.35, ttaui= 0.29, ttaug= 0.12, ttaul=0.28);

# Very nice
pa  = init_parameters( H=1.3, ff= 0.15, llambda0=0.02, llambda1= 0.04, ddelta = 0.14, ttheta = 0.3,rhoz= 0.76, ssigmaz= 0.0385);
comment:"Simulations are decent for large shocks. The selection mechanism is not doing much: productive firms stay and unproductive firms exit,
  independent of size. Investment is a little bit low."

#Decrease support of shocks
comment:"Turnover decreases, decreasing mean and freq equity shares. And increasing investment. We get some equity issuances by incumbents"
# Recalibrate
pa  = init_parameters( H=1.3, ff= 0.15, llambda0=0.02, llambda1= 0.04, ddelta = 0.14, ttheta = 0.3,rhoz= 0.76, ssigmaz= 0.045);
# The wage decreases substancially and profits increase
comment:"Once recalibrated, the moments look good (except the autocorr of profits, which increases but is still too low)"

# Decrease fix cost but increase e
pa  = init_parameters( H=1.3, ff= 0.05, llambda0=0.02, llambda1= 0.04, ddelta = 0.14, ttheta = 0.3,rhoz= 0.76, ssigmaz= 0.045, e=0.08);
comment:"Turnover decreases, decreasing mean and freq equity shares. And increasing investment. We get some equity issuances by incumbents
wage is 0.55"
pa  = init_parameters( H=1.3, ff= 0.05, llambda0=0.02, llambda1= 0.04, ddelta = 0.095, ttheta = 0.3,rhoz= 0.76, ssigmaz= 0.09, e=0.08);
comment:"Moments are very nice: but growing out of constraints mechanism is not present. Wage =0.67, Which helps with ave profits "

pa  = init_parameters( H=1.5, ff= 0.01, llambda0=0.06, llambda1= 0.04, ddelta = 0.095, ttheta = 0.25,rhoz= 0.76, ssigmaz= 0.09, e=0.08);



pa  = init_parameters( H=1.32, bbeta=0.972, ff= 0.15, aalphak=0.23, aalphal=0.64, llambda0=0.05, llambda1= 0.06, ddelta = 0.12, allowance=0.2, ttheta = 0.3,rhoz= 0.76, ssigmaz= 0.085, e=0.1, A=1.0);
tau = init_taxes(ttaud =0.12, ttauc= 0.35, ttaui= 0.29, ttaug= 0.12, ttaul=0.28);
comment:"Turnover is very low "

#First fix leverage and H
pa  = init_parameters( H=1.4, bbeta=0.972, ff= 0.15, aalphak=0.23, aalphal=0.64, llambda0=0.05, llambda1= 0.06, ddelta = 0.12, allowance=0.2, ttheta = 0.25,rhoz= 0.76, ssigmaz= 0.085, e=0.1, A=1.0);
comment:"Turnover is very low. I need more small firm.  "

#Increase lambda0, lambda1
pa  = init_parameters( H=1.4, bbeta=0.972, ff= 0.15, aalphak=0.23, aalphal=0.64, llambda0=0.08, llambda1= 0.08, ddelta = 0.12, allowance=0.2, ttheta = 0.25,rhoz= 0.76, ssigmaz= 0.085, e=0.1, A=1.0);
comment:"lambda0 desn't do much. lambda1 actually decreases turnover "

#Forget about turnover
comment:"I need a new moment to match f"
#and try to match equity isuuances

#THIS IS THE BEST CALIBRATION.
pa  = init_parameters( H=1.4, bbeta=0.972, ff= 0.15, aalphak=0.23, aalphal=0.64, llambda0=0.02, llambda1= 0.04, ddelta = 0.12, allowance=0.2, ttheta = 0.25,rhoz= 0.75, ssigmaz= 0.085, e=0.0, A=1.0);
comment:"Frecuency of equity issuances is slightly low. Everything else looks very nice (except for tax base ratio)"

#decrease lambda0
pa  = init_parameters( H=1.4, bbeta=0.972, ff= 0.15, aalphak=0.23, aalphal=0.64, llambda0=0.01, llambda1= 0.04, ddelta = 0.12, allowance=0.2, ttheta = 0.25,rhoz= 0.75, ssigmaz= 0.085, e=0.0, A=1.0);
comment:"I am happy enough with this calibration"

#Try to match tax base ratio
pa  = init_parameters( H=1.4, bbeta=0.972, ff= 0.15, aalphak=0.23, aalphal=0.64, llambda0=0.01, llambda1= 0.04, ddelta = 0.12, allowance=0.0, ttheta = 0.25,rhoz= 0.75, ssigmaz= 0.085, e=0.0, A=1.0);
comment:"can't match the tax base, even with 0.0 allowance"
