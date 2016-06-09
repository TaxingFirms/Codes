#Presentation

normalized_distr = p.distr/sum(p.distr);
prod = repmat(fp.zgrid,1,fp.Nomega)';
avgprod = sum(normalized_distr.*prod)

avgK =  sum(normalized_distr.*res.kprime)



capital = sum(p.distr.*res.kprime)

p.output

tfp=p.output/(capital^fp.alphak*p.laborsupply^fp.alphal)

p.w
p.consumption
p.laborsupply

delta_welfare = (p.welfare/(1-hp.beta) + 58.465859102488174)/(58.465859102488174)*100

tau
sum([p.collections.i p.collections.c p.collections.d])



##########
# Model taxing liquidations
hp = init_hhparameters();

pr,tau,fp,res,p=load("/Users/danielwillsr/Dropbox/1FirmTaxation/SimpleDiscreteTime/ModelResultsTxExit.jld", "pr","tau","fp","res","p");

pr,tau,fp,res,p=load("/Users/danielwillsr/Dropbox/1FirmTaxation/SimpleDiscreteTime/Reform03TxExit.jld", "pr","tau","fp","res","p");

pr,tau,fp,res,p=load("/Users/danielwillsr/Dropbox/1FirmTaxation/SimpleDiscreteTime/Reform04TxExit.jld", "pr","tau","fp","res","p");



