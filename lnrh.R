
###This script calculates lnRH values and quantifies how the statistic changes from the least to the most divergent genes.


setwd("/home/seb/Documents/repeatability") #set up working directory #set up working directory 

comp5_1 =  read.table("results_123/deb_pet_fst_map",  header = T, stringsAsFactors = F); comp5_1 = cbind(comp5_1,0,0,0);colnames(comp5_1)[(ncol(comp5_1)-2):ncol(comp5_1)] = c("ho1","ho2","lnRH")
comp5_2 =  read.table("results_123/bol_exi_fst_map",  header = T, stringsAsFactors = F); comp5_2 = cbind(comp5_2,0,0,0);colnames(comp5_2)[(ncol(comp5_2)-2):ncol(comp5_2)] = c("ho1","ho2","lnRH")
comp5_3 =  read.table("results_123/ann_arg_fst_map",  header = T, stringsAsFactors = F); comp5_3 = cbind(comp5_3,0,0,0);colnames(comp5_3)[(ncol(comp5_3)-2):ncol(comp5_3)] = c("ho1","ho2","lnRH")

name = rbind(c("_pet","_deb"),c("_EAST","_WEST"),c("_ann","_arg"))
for(i in 1:3)
{
		if(i == 1) comp = comp5_1 #pet-deb
		if(i == 2) comp = comp5_2 #bol-exi
		if(i == 3) comp = comp5_3 #ann-arg
		
		comp1 = comp[,colnames(comp)[regexpr(name[i,1],colnames(comp)) >0]]
		comp2 = comp[,colnames(comp)[regexpr(name[i,2],colnames(comp)) >0]]
		ho  = cbind(rep(0,nrow(comp)),rep(0,nrow(comp)), rep (0,nrow(comp)))
				for(j in 1:nrow(comp))
					{
				
					a1 = substring(comp1[j,],1,1); a2 = substring(comp1[j,],2,2); a1 = c(a1,a2)
					b1 = substring(comp2[j,],1,1); b2 = substring(comp2[j,],2,2); b1 = c(b1,b2)
				
					a1_val = rle(sort(a1))$values[rle(sort(a1))$values != "X"]
					a1_len = rle(sort(a1))$length[rle(sort(a1))$values != "X"]
					b1_val = rle(sort(b1))$values[rle(sort(b1))$values != "X"]
					b1_len = rle(sort(b1))$length[rle(sort(b1))$values != "X"]
					pop1_2pq = (a1_len[1] / sum(a1_len)) * (a1_len[2] / sum(a1_len)) * 2 #expected hetero pop 1
					pop2_2pq = (b1_len[1] / sum(b1_len)) * (b1_len[2] / sum(b1_len)) * 2

					#expected hetero
					ho[j,1] = pop1_2pq #x
					ho[j,2] = pop2_2pq #y
					ho[j,3] = log(((1/(1-pop1_2pq)^2)-1)/ ((1/(1-pop2_2pq)^2)-1)) #   log(((1/(1-x)^2)-1)/ ((1/(1-y)^2)-1)) #lnRH
					if(is.na(ho[j,3])) ho[j,3] = Inf
					if(j %% 100000 == 0) print(paste(j,"of",nrow(comp),Sys.time()))
				}
		if(i == 1) comp5_1[,(ncol(comp5_1)-2):ncol(comp5_1)] = ho
		if(i == 2) comp5_2[,(ncol(comp5_2)-2):ncol(comp5_2)] = ho
		if(i == 3) comp5_3[,(ncol(comp5_3)-2):ncol(comp5_3)] = ho
}

###
###lnrh statistic decreases/increases of the highly divergent genes?
###
results = matrix(0,nrow = 20,ncol = 2)
ccc = comp5_1[comp5_1[,51] != Inf,]; cc_pd = ccc[ccc[,51] != -Inf,];cc_pd = cbind(cc_pd,0) #; t.test(cc[,51])
ccc = comp5_3[comp5_3[,64] != Inf,]; cc_aa = ccc[ccc[,64] != -Inf,];cc_aa = cbind(cc_aa,0) #; t.test(cc[,64])

unique_genes = read.delim("results_123/unique_genes", header = T, sep = " ", stringsAsFactors = F)
unique_genes[is.na(unique_genes[,2]),2] = 0 
unique_genes[is.na(unique_genes[,3]),3] = 0 
unique_genes[is.na(unique_genes[,4]),4] = 0 
outliers_genes = matrix(0,nrow = 20,ncol = 6)
top = seq(0,1,by = 0.05)
un = unique_genes[unique_genes[,2] != 0,]; un = un[un[,4] != 0,];outliers_genes[t,5] = top[2] * top[2] * nrow(un) #e

for(t in 1:20)
{
q1_a = quantile(un[,2],top[t],na.rm = T);q1_b = quantile(un[,2],top[t+1],na.rm = T)
q2_a = quantile(un[,4],top[t],na.rm = T);q2_b = quantile(un[,4],top[t+1],na.rm = T)

cc_pd[(cc_pd[,40] > q1_a) & (cc_pd[,40] < q1_b),52] = t
cc_aa[(cc_aa[,53] > q2_a) & (cc_aa[,53] < q2_b),65] =  t

results[t,1] = mean(cc_pd[(cc_pd[,40] > q1_a) & (cc_pd[,40] < q1_b),51])
results[t,2] = mean(cc_aa[(cc_aa[,53] > q2_a) & (cc_aa[,53] < q2_b),64])
}


###plot lnrh for the quantiles
par(mar = c(6,6,4,2))
boxplot(cc_pd[,51]~ cc_pd[,52], at = seq(0,61, by = 3), xaxt = "n", outline = F,col = "lightblue", ylim = c(-10,10))
boxplot(cc_aa[,64]~ cc_aa[,65], add = T, at = seq(1,62, by = 3), xaxt = "n", outline = F, col = "red")
axis(1,at = seq(0.5,61, by = 3), labels = paste( seq(0,100,by = 5), "%", sep = ""), las = 2)
axis(1, at = 30.5,labels = "Quantiles",font.axis = 2, line = 3,cex.axis = 2, tick = F)
axis(2, at = 0,labels =  "lnRH statistic",font.axis = 2, line = 2.5,cex.axis = 2, tick = F)
legend(30,10, legend = c("H. petiolaris - H. debilis", "H. annuus - H. argophyllus"), bg = "white", text.font = 4, fill = c("lightblue","red"), border = c("black","black"))
abline(h = 0)

dev.print(device=pdf, "results_6species/lnrh.pdf", onefile=FALSE)
dev.off()



