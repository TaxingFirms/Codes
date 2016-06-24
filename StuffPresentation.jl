#Presentation



pr,tau,fp,res,p=load("/Users/danielwillsr/Dropbox/1FirmTaxation/Results/ModelResults.jld", "pr","tau","fp","res","p");
pr1,tau1,fp1,res1,p1=load("/Users/danielwillsr/Dropbox/1FirmTaxation/Results/counterfactual1.jld", "pr","tau","fp","res","p");
pr1rn,tau1rn,fp1rn,res1rn,p1rn=load("/Users/danielwillsr/Dropbox/1FirmTaxation/Results/counterfactual1RN.jld", "pr","tau","fp","res","p");
pr2,tau2,fp2,res2,p2=load("/Users/danielwillsr/Dropbox/1FirmTaxation/Results/counterfactual2.jld", "pr","tau","fp","res","p");
pr2rn,tau2rn,fp2rn,res2rn,p2rn=load("/Users/danielwillsr/Dropbox/1FirmTaxation/Results/counterfactual2RN.jld", "pr","tau","fp","res","p");


#tfp=p.output/(capital^fp.alphak*p.laborsupply^fp.alphal)
