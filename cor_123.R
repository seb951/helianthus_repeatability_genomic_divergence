#!/usr/bin/Rscript --verbose

###This script looks at overall corellation coefficients at 3 level of genome organisation (SNP, gene, genome level). 
###It also looks at the effect of recombination rate on patterns of Fst



####################
######
### MAPPING FST VALUES ###
##########################
setwd("/home/seb/Documents/repeatability") #set up working directory #set up working directory 

################
### PLOTTING ###
################
c1 =  read.table("results_123/deb_pet_fst_map_cluster_nonoverlapping_windows",  header = T, stringsAsFactors = F);c1[c1[,3] == 0,3] = NA #non overlapping windows (to avoid autocorrelation problems...)
c2 =  read.table("results_123/bol_exi_fst_map_cluster_nonoverlapping_windows",  header = T, stringsAsFactors = F);c2[c2[,3] == 0,3] = NA
c3 =  read.table("results_123/ann_arg_fst_map_cluster_nonoverlapping_windows",  header = T, stringsAsFactors = F);c3[c3[,3] == 0,3] = NA

#############
###per region ###
############
par(mfrow = c(2,2))
c12 = cbind(c1[,c(1,3)],c2[,3]) # ; c12 = c12[c12[,2] != 0,]; c12 = c12[c12[,3] != 0,];
plot(c12[,3],c12[,2], xlab = "bolexi", ylab = "debpet", xlim = c(-0.1,0.8), ylim = c(-0.1,1))      #smoothScatter(c1[c1[,3] != 0,3],c2[c1[,3] != 0,3])
text(0.6,0.4, labels = signif(cor.test(c12[,2],c12[,3], alternative = "g")$estimate,2), col = "darkred",cex = 3)

c23 = cbind(c2[,c(1,3)],c3[,3]) #; c23 = c23[c23[,2] != 0,]; c23 = c23[c23[,3] != 0,];
plot(c23[,2],c23[,3], xlab = "bolexi", ylab = "annarg", xlim = c(-0.1,0.8), ylim = c(-0.1,1)) 
text(0.4,0.6, labels = signif(cor.test(c23[,2],c23[,3], alternative = "g")$estimate,2), col = "darkred",cex = 3)

c13 = cbind(c1[,c(1,3)],c3[,3]) #; c13 = c13[c13[,2] != 0,]; c13 = c13[c13[,3] != 0,];
plot(c13[,2],c13[,3], xlab = "debpet", ylab = "annarg", xlim = c(-0.1,0.8), ylim = c(-0.1,1)) 
text(0.6,0.6, labels = signif(cor.test(c13[,2],c13[,3], alternative = "g")$estimate,2), col = "darkred",cex = 3)

### ### ### ### ### ### ###
### per chromosome BARPLOT ###
### ### ### ### ### ### ###
cor_per_chr = NULL

for(comp in 1:3)
	{
	cor_per_chr = NULL
	for(chr in 1:17)
		{
		if(comp == 3) cor_per_chr = rbind(cor_per_chr,c(chr,signif(cor.test(c12[c1[,2] == chr,2],c12[c1[,2] == chr,3], alternative = "g")$estimate,2)))
		if(comp == 2) cor_per_chr = rbind(cor_per_chr,c(chr,signif(cor.test(c23[c1[,2] == chr,2],c23[c1[,2] == chr,3], alternative = "g")$estimate,2)))
		if(comp == 1) cor_per_chr = rbind(cor_per_chr,c(chr,signif(cor.test(c13[c1[,2] == chr,2],c13[c1[,2] == chr,3], alternative = "g")$estimate,2)))
		}
	main = c("deb_pet vs ann_arg","bol_exi vs ann_arg","bol_exi vs deb_pet")
	if(comp == 1) plot((cor_per_chr[,1]-0.2),cor_per_chr[,2], type = "h", ylim = c(-0.1,0.9), ,col = "darkblue", lwd = 4,xaxt = "n",xlab = "linkage group", ylab =  "Pearson's product-moment correlation (r)")
	if(comp == 2) points((cor_per_chr[,1]),cor_per_chr[,2], type = "h", lwd = 4,xaxt = "n",col = "darkgreen")
	if(comp == 3) points((cor_per_chr[,1]+0.2),cor_per_chr[,2], type = "h", lwd = 4,xaxt = "n",col = "darkred")
}
axis(1,c(1:17),c(1:17))
legend(y = 0.9,x = 9,legend = main, fill = c("darkblue","darkgreen","darkred"), cex = 1.0)

dev.print(device=svg, "results_123/cor_per_region2.svg",onefile=FALSE)
dev.off()
############
###per gene ###
############
comp5_1 =  read.table("results_123/deb_pet_fst_map",  header = T, stringsAsFactors = F)
comp5_2 =  read.table("results_123/bol_exi_fst_map",  header = T, stringsAsFactors = F)
comp5_3 =  read.table("results_123/ann_arg_fst_map",  header = T, stringsAsFactors = F)

#################################
#unique_genes = cbind(as.data.frame(unique(sort(c(comp5_1[,1],comp5_2[,1],comp5_3[,1]))), stringsAsFactors = F),NA,NA,NA)
colnames(unique_genes) = c("reference","mean_fst_deb_pet","mean_fst_bol_exi","mean_fst_ann_arg","numberSNP_deb_pet","numberSNP_fst_bol_exi","numberSNP_fst_ann_arg")
unique_genes = cbind(as.data.frame(unique(sort(c(comp5_1[,1],comp5_2[,1],comp5_3[,1]))), stringsAsFactors = F),NA,NA,NA,NA,NA,NA)

#for(i in 6000:6100)

	for(i in 1:nrow(unique_genes)) 
	{
	x1 = comp5_1[comp5_1[,1] == unique_genes[i,1],ncol(comp5_1)-8] #fst
	x1_cleaned = x1[x1 > 0]
	if(length(x1_cleaned) > 0) {unique_genes[i,2] = signif(mean(x1_cleaned),3); unique_genes[i,5] = length(x1_cleaned)}# else  {unique_genes[i,2] = 0} 
	#
	x2 = comp5_2[comp5_2[,1] == unique_genes[i,1],ncol(comp5_2)-8] #fst
	x2_cleaned = x2[x2 > 0]
	if(length(x2_cleaned) > 0) {unique_genes[i,3] = signif(mean(x2_cleaned),3);unique_genes[i,6] = length(x2_cleaned)}  #else  {unique_genes[i,3] = 0} 
	#
	x3 = comp5_3[comp5_3[,1] == unique_genes[i,1],ncol(comp5_3)-8] #fst
	x3_cleaned = x3[x3 > 0]
	if(length(x3_cleaned) > 0) {unique_genes[i,4] = signif(mean(x3_cleaned),3);unique_genes[i,7] = length(x3_cleaned)} #else  {unique_genes[i,4] = 0} 
	#
	if(i %% 1000 == 0) print(paste(i,"of",nrow(unique_genes), Sys.time()))
	}
write.table(unique_genes,"results_123/unique_genes", col.names = T, row.names = F, quote = F)
	
###per gene
unique_genes_ori = read.table("results_123/unique_genes", header = T, stringsAsFactors = F)

#unique_genes[is.na(unique_genes[,2]),2] = 0
#unique_genes[is.na(unique_genes[,3]),3] = 0
#unique_genes[is.na(unique_genes[,4]),4] = 0
#unique_genes[(unique_genes[,2]) ==  0,2] = NA
#unique_genes[(unique_genes[,3]) == 0,3] = NA
#unique_genes[(unique_genes[,4]) == 0,4] = NA
attach(unique_genes)

par(mfrow = c(2,2))
unique_genes1 = unique_genes[!is.na(unique_genes[,2]), ];  unique_genes1 = unique_genes1[!is.na(unique_genes1[,3]), ]; smoothScatter(mean_fst_bol_exi,mean_fst_deb_pet);cor(unique_genes1[,2],unique_genes1[,3])
text(0.8,0.8, labels = signif(cor.test(unique_genes[,2],unique_genes[,3], alternative = "g")$estimate,2), col = "darkred",cex = 3)
unique_genes2 = unique_genes[!is.na(unique_genes[,3]), ];  unique_genes2 = unique_genes2[!is.na(unique_genes2[,4]), ]; smoothScatter(mean_fst_bol_exi,mean_fst_ann_arg);cor(unique_genes2[,3],unique_genes2[,4])
text(0.8,0.8, labels = signif(cor.test(unique_genes[,3],unique_genes[,4], alternative = "g")$estimate,2), col = "darkred",cex = 3)
unique_genes3 = unique_genes[!is.na(unique_genes[,2]), ];  unique_genes3 = unique_genes3[!is.na(unique_genes3[,4]), ]; smoothScatter(mean_fst_deb_pet,mean_fst_ann_arg);cor(unique_genes3[,2],unique_genes3[,4])
text(0.8,0.8, labels = signif(cor.test(unique_genes[,2],unique_genes[,4], alternative = "g")$estimate,2), col = "darkred",cex = 3)

dev.print(device=svg, "results_123/cor_pergene.svg",onefile=FALSE)
dev.off()


### ### ### ###
### per SNP ###
### ### ### ###
if(file.exists("results_123/unique_SNP") == F) {

comp5_1_uniq = paste(comp5_1[,1],comp5_1[,2], sep = "_");
comp5_2_uniq = paste(comp5_2[,1],comp5_2[,2], sep = "_");
comp5_3_uniq = paste(comp5_3[,1],comp5_3[,2], sep = "_");
unique_SNP = cbind(as.data.frame(unique(sort(c(comp5_1_uniq,comp5_2_uniq,comp5_3_uniq))), stringsAsFactors = F),NA,NA,NA);

colnames(unique_SNP) = c("reference","fst_deb_pet","fst_bol_exi","fst_ann_arg");

	for(i in 1:nrow(unique_SNP)) 
#	for(i in c(1:1000))
	{
	x1 = comp5_1[comp5_1_uniq == unique_SNP[i,1],ncol(comp5_1)-8] #fst
	if(length(x1) >1) x1 = x1[1]
	if(length(x1) > 0) unique_SNP[i,2] = x1 # else  {unique_SNP[i,2] = 0} 
	#
	x2 = comp5_2[comp5_2_uniq == unique_SNP[i,1],ncol(comp5_2)-8] #fst
	if(length(x2) >1) x2 = x2[1]
	if(length(x2) > 0) unique_SNP[i,3] = x2  #else  {unique_SNP[i,3] = 0} 
	#
	x3 = comp5_3[comp5_3_uniq == unique_SNP[i,1],ncol(comp5_3)-8] #fst
	if(length(x3) >1) x3 = x3[1]
	if(length(x3) > 0) unique_SNP[i,4] = x3  #else  {unique_SNP[i,4] = 0} 
	#
	if(i %% 100 == 0) print(paste(i,"of",nrow(unique_SNP), " >The time is...",Sys.time()))
	}
	
	write.table(unique_SNP,"results_123/unique_SNP", col.names = T, row.names = F, quote = F)
	}
###plot per SNP
unique_SNP = read.table("results_123/unique_SNP", header = T, stringsAsFactors = F)
attach(unique_SNP)
cc = rbind(c(3,2),c(3,4),c(2,4))

#unique_SNP[is.na(unique_SNP[,2]),2] = 0
#unique_SNP[is.na(unique_SNP[,3]),3] = 0
#unique_SNP[is.na(unique_SNP[,4]),4] = 0

#png("results_123/cor_perSNP.png",height = 2000, width = 2000)

par(mfrow = c(2,2))
#par(mar = c(8, 7, 6, 4))
for(i in 1:3)
{
unique_SNP_small = unique_SNP[!is.na(unique_SNP[,cc[i,1]]),]; unique_SNP_small = unique_SNP_small[!is.na(unique_SNP_small[,cc[i,2]]),];  cor(unique_SNP_small[,cc[i,1]],unique_SNP_small[,cc[i,2]]);smoothScatter(unique_SNP_small[,cc[i,1]],unique_SNP_small[,cc[i,2]],xlab = colnames(unique_SNP_small)[(cc[i,1])],ylab = colnames(unique_SNP_small)[(cc[i,2])], cex.lab = 1, ylim = c(-0.2,1.1), xlim = c(-0.2,1.1))
#abline(lm(unique_SNP_small[,cc[i,2]] ~ unique_SNP_small[,cc[i,1]], data = unique_SNP_small), lwd = 4,col = "darkred", lty = 2)
text(0.5,0.5, labels = signif(cor.test(unique_SNP_small[,cc[i,1]],unique_SNP_small[,cc[i,2]], alternative = "g")$estimate,2), cex = 1)
}

dev.print(device=svg, "results_123/cor_perSNP2.svg",onefile=FALSE)
#dev.print(device=png, "results_123/cor_perSNP.png",height = 2000, width = 2000)
dev.off()

### ### ### ### ### ### ### ### ### 
### plot recombination rates against FST ###
### ### ### ### ### ### ### ### ###
c1 =  read.table("results_123/deb_pet_fst_map_cluster_nonoverlapping_windows",  header = T, stringsAsFactors = F);c1[c1[,3] == 0,3] = NA #non overlapping windows (to avoid autocorrelation problems...)
c2 =  read.table("results_123/bol_exi_fst_map_cluster_nonoverlapping_windows",  header = T, stringsAsFactors = F);c2[c2[,3] == 0,3] = NA
c3 =  read.table("results_123/ann_arg_fst_map_cluster_nonoverlapping_windows",  header = T, stringsAsFactors = F);c3[c3[,3] == 0,3] = NA


unique_map_transcript = read.delim("~/Documents/speciation_islands_dec2011/reference/recombination_rates", header = T, stringsAsFactors = F,sep = " ") ###recombination rates###
main2 = c("deb_pet","bol_exi","ann_arg")	
#color = c("#f8756b25","#7aae0025","#00bdc225") ; color2 = c("#f8756b","#7aae00","#00bdc2")
color = c("black","black","black");color2 = c("black","black","black")
par(mfrow = c(2,2))

plot(c1[,3]~log(unique_map_transcript[,6]), col = color[1], xlim = c(-8,5.5), pch = 21, xlab = "recom_rate (cM/MB)",ylab = "Fst",xaxt = "n");
axis(1, at = log(c(0.001,0.01,0.1,1,10,100)),c(0.001,0.01,0.1,1,10,100))
text(cex = 1.2,x = log(0.01), y =0.6,font = 2, labels = paste("cor. coef. = ",signif(cor.test(c1[,3],log(unique_map_transcript[,6]))$estimate,4)," (",main2[1],")",sep = ""),col = color2[1])
lm_regression_1 = lm(c1[,3]~log(unique_map_transcript[,6]))
#abline(lm_regression_1, lwd = 6,col = color2[1])

plot(c2[,3]~log(unique_map_transcript[,6]), col = color[2], xlim = c(-8,5.5), pch = 21, xlab = "recom_rate (cM/MB)",ylab = "Fst",xaxt = "n");
axis(1, at = log(c(0.001,0.01,0.1,1,10,100)),c(0.001,0.01,0.1,1,10,100))
text(cex = 1.2,x = log(0.01), y =0.6,font = 2, labels = paste("cor. coef. = ",signif(cor.test(c2[,3],log(unique_map_transcript[,6]))$estimate,4)," (",main2[2],")",sep = ""),col = color2[2])
lm_regression_1 = lm(c2[,3]~log(unique_map_transcript[,6]))
#abline(lm_regression_1, lwd = 6,col = color2[2])

plot(c3[,3]~log(unique_map_transcript[,6]), col = color[3], xlim = c(-8,5.5), pch = 21, xlab = "recom_rate (cM/MB)",ylab = "Fst",xaxt = "n");
axis(1, at = log(c(0.001,0.01,0.1,1,10,100)),c(0.001,0.01,0.1,1,10,100))
text(cex = 1.2,x = log(0.01), y =0.6,font = 2, labels = paste("cor. coef. = ",signif(cor.test(c3[,3],log(unique_map_transcript[,6]))$estimate,4)," (",main2[3],")",sep = ""),col = color2[3])
lm_regression_1 = lm(c3[,3]~log(unique_map_transcript[,6]))
#abline(lm_regression_1, lwd = 6,col = color2[3])

dev.print(device=svg, "results_123/cor_recombrate2.svg", onefile=FALSE)
dev.off()



########################################################################################################
###SANDBOX#######SANDBOX#######SANDBOX#######SANDBOX#######SANDBOX#######SANDBOX#######SANDBOX#######SANDBOX##
######################################################################################################

###	###
###	per gene dnds 
###	###
setwd("/home/seb/Documents/repeatability") #set up working directory #set up working directory 

dnds_deb_pet = (read.table("results_123/deb_pet_fst_map_dnds.out",stringsAsFactors = F,header = T))
dnds_bol_exi = (read.table("results_123/bol_exi_fst_map_dnds.out",stringsAsFactors = F, header = T))
dnds_ann_arg = (read.table("results_123/ann_arg_fst_map_dnds.out",stringsAsFactors = F, header = T))

dnds = unique(sort(c(dnds_deb_pet[,1],dnds_bol_exi[,1],dnds_ann_arg[,1])))
dnds = cbind(dnds,0,0,0,0,0,0,0,0,0); dnds = as.data.frame(dnds, stringsAsFactors = F);dnds[,2:10] = 0
colnames(dnds) = c("name","dn1","ds1","dnds1","dn2","ds2","dnds2","dn3","ds3","dnds3")
for(i in 1:nrow(dnds))
{
	temp1= dnds_deb_pet[dnds_deb_pet[,1] == dnds[i,1],2:4]
	temp2= dnds_bol_exi[dnds_bol_exi[,1] == dnds[i,1],2:4]
	temp3= dnds_ann_arg[dnds_ann_arg[,1] == dnds[i,1],2:4]
	
	if(length(temp1) > 0) dnds[i,2:4] = as.numeric(temp1)
	if(length(temp2) > 0) dnds[i,5:7] = as.numeric(temp2)
	if(length(temp3) > 0) dnds[i,8:10] = as.numeric(temp3)
	
	if(i %% 1000 == 0) print(paste(i,"of", nrow(dnds),Sys.time()))
}

dnds[is.na(dnds[,2]),2] = 0
dnds[is.na(dnds[,3]),3] = 0
dnds[is.na(dnds[,4]),4] = 0
dnds[is.na(dnds[,5]),5] = 0
dnds[is.na(dnds[,6]),6] = 0
dnds[is.na(dnds[,7]),7] = 0
dnds[is.na(dnds[,8]),8] = 0
dnds[is.na(dnds[,9]),9] = 0
dnds[is.na(dnds[,10]),10] = 0

dnds_log = dnds
for(i in 2:ncol(dnds))
{dnds_log[,i] = log(dnds[,i])}
MAYBE THE CORRELATIONS WITH BOL EXI ARE DRIVEN BY DS?

###dn
par(mfrow = c(2,2))
plot(dnds_log[(dnds_log[,2] != -Inf) & (dnds_log[,5] != -Inf),2],  dnds_log[(dnds_log[,5] != -Inf) & (dnds_log[,2] != -Inf),5], main = "dnds correlation per gene for bol_exi and pet_deb") 
text(-4,-4,signif(cor(dnds_log[(dnds_log[,5] != -Inf) & (dnds_log[,2] != -Inf),5],  dnds_log[(dnds_log[,5] != -Inf) & (dnds_log[,2] != -Inf),2])), cex = 2)
plot(dnds_log[(dnds_log[,5] != -Inf) & (dnds_log[,8] != -Inf),8],  dnds_log[(dnds_log[,8] != -Inf) & (dnds_log[,5] != -Inf),5], main = "dnds correlation per gene for ann_arg and bol_exi") 
text(-4,-4,signif(cor(dnds_log[(dnds_log[,8] != -Inf) & (dnds_log[,5] != -Inf),5],  dnds_log[(dnds_log[,8] != -Inf) & (dnds_log[,5] != -Inf),8])), cex = 2)
plot(dnds_log[(dnds_log[,8] != -Inf) & (dnds_log[,2] != -Inf),2],  dnds_log[(dnds_log[,8] != -Inf) & (dnds_log[,2] != -Inf),8], main = "dnds correlation per gene for ann_arg and pet_deb") 
text(-4,-4,signif(cor(dnds_log[(dnds_log[,8] != -Inf) & (dnds_log[,2] != -Inf),2],  dnds_log[(dnds_log[,8] != -Inf) & (dnds_log[,2] != -Inf),8])), cex = 2)

}

### ### ### ### ### ### ### ### ###
### per genomic region, on a sliding window scale ###
### ### ### ### ### ### ### ### ###

map = read.delim("reference/ordered_transcriptome_trinity.txt", header = F, stringsAsFactors = F) # new map
map = map[map[,3] != "no",]

#unique position for the markers#
max_position = cbind(c(1:17), 0,0) # the size of each chromosome and the cumulative number of centimorgans.

for(i in 1:17)
	{
	max_position[i,2] = max(as.numeric(map[as.numeric(map[,2]) == i,3]))
	max_position[i,3] = sum(max_position[1:i,2])-max_position[i,2]
	}

if(pair == 1) fst_map_match = read.delim("results_123/deb_pet_fst_map", header = T, sep = " ", stringsAsFactors = F) # ann_pet_fst
if(pair == 2) fst_map_match = read.delim("results_123/bol_exi_fst_map", header = T, sep = " ", stringsAsFactors = F) # pet_deb_fst
if(pair == 3) fst_map_match = read.delim("results_123/ann_arg_fst_map", header = T, sep = " ", stringsAsFactors = F) # ann_deb_fst
c1 =  read.table("results_123/deb_pet_fst_map_cluster",  header = T, stringsAsFactors = F);c1[c1[,3] == 0,3] = NA
rearrang = read.delim("~/Documents/chromosomal_rearrangments_johnjohn/reference/rearrang.txt", header = T) #rearrang
unique_map_transcript = read.delim("~/Documents/speciation_islands_dec2011/reference/recombination_rates", header = T, stringsAsFactors = F,sep = " ") ###recombination rates###
	col = ncol(fst_map_match) 
	fst_map_match[is.na(fst_map_match[,(col-8)]),(col-8):(col-6)] = 0
	top_q = quantile(fst_map_match[,(col-8)],0.97) #top3%
	
	####running window - distance based###
 	d = 5 #distance cutoff in cM
	
	c_cor = c1
	c_cor[,3:9] = 0 #position you are interrogating along with the result of the distance based running mean. 
	colnames(c_cor)[1:8] = c("unique_position","LG","cor12","cor23","cor13","ndata","genedensity_swscale","recombrate")
		for(i in 1:nrow(c_cor)) 
			{
			temp1 = c1[c1[,1] > (c1[i,1]-d) & c1[,1] < (c1[i,1]+d),]; temp1 = temp1[temp1[,2] == c1[i,2],3]
			temp2 = c2[c2[,1] > (c1[i,1]-d) & c2[,1] < (c1[i,1]+d),]; temp2 = temp2[temp2[,2] == c1[i,2],3]
			temp3 = c3[c3[,1] > (c1[i,1]-d) & c3[,1] < (c1[i,1]+d),]; temp3 = temp3[temp3[,2] == c1[i,2],3]
		    temp4 = c1[c1[,1] > (c1[i,1]-d) & c1[,1] < (c1[i,1]+d),]; temp4 = temp4[temp4[,2] == c1[i,2],9]
		    temp5 = unique_map_transcript[c1[,1] > (c1[i,1]-d) & c1[,1] < (c1[i,1]+d),]; temp5 = temp5[temp5[,2] == c1[i,2],6]
			if(length(temp1[!is.na(temp1) & !is.na(temp2)]) > 2) c_cor[i,3] = cor.test(temp1,temp2)$estimate    #cor
			if(length(temp1[!is.na(temp2) & !is.na(temp3)]) > 2) c_cor[i,4] = cor.test(temp2,temp3)$estimate	#cor
			if(length(temp1[!is.na(temp1) & !is.na(temp3)]) > 2) c_cor[i,5] = cor.test(temp1,temp3)$estimate	#cor
			c_cor[i,6] = length(temp1[!is.na(temp1)]) #nb unique regions in correlation
			c_cor[i,7] = sum(temp4) #nb nucleotide
			c_cor[i,8] = mean(temp5,na.rm = T) #mean recomb. rate
			if(length(temp1[!is.na(temp1)]) > 2) if(i %% 1000 == 0) print(paste(i, Sys.time()))
			}
			
	main = c("deb_pet vs bol_exi","bol_exi vs ann_arg","deb_pet vs ann_arg")		
###plotting
###
par(mfrow = c(3,1))
for(pair in c(1,2,3))
{
color = c("#f8756b","#7aae00","#00bdc2")
###fst outlier regions

for(i in 1:17) #17 linkage groups
{
if(pair == 1) cc1= c_cor[c_cor[,2] == i, c(1,3,7,8)]
if(pair == 2) cc1= c_cor[c_cor[,2] == i, c(1,4,7,8)]
if(pair == 3) cc1= c_cor[c_cor[,2] == i, c(1,5,7,8)]
r_plot = rearrang[rearrang[,1] == i,]###rearrangments###

	if(i != 17) cc2 = c_cor[c_cor[,2] == (i+1), ] else {cc2 = cc1[cc1[,1] == max(cc1[,1]),]; cc2 = rbind(cc2,cc2)}
	if(i == 1) plot(0,1,  type = "l", xlim = c(0,max(c1[,1])+1), ylim = c(-1.09,1.09), xaxt = "n", lwd = 4,col = "black", main = main[pair], xlab = "", ylab = "")
	
	if(i == 1) rect(-1,-1,min(cc2[cc2[,3]>0,1]),1, border = F, col = ifelse(i %% 2 == 0,colors()[606],colors()[334]))
	if(i != 1) rect(min(cc1[,1]),-1,min(cc2[,1]),1, border = F, col = ifelse(i %% 2 == 0,colors()[606],colors()[334]))
	
	#points(cc1[,1],1-(cc1[,4])/max((c_cor[,8]),na.rm = T),type = "h", col = "#00009925") #number of nucleotide
	
	if(nrow(r_plot) ==1) { abline(v= min(cc1[,1]) + r_plot[1,2],lty = 2);abline(v= min(cc1[,1]) + r_plot[1,3],lty = 2)}
if(nrow(r_plot) ==2) { abline(v= min(cc1[,1]) + r_plot[1,2],lty = 2);abline(v= min(cc1[,1]) + r_plot[1,3],lty = 2);
abline(v= min(cc1[,1]) + r_plot[2,2],lty = 2);abline(v= min(cc1[,1]) + r_plot[2,3],lty = 2)} #plot rearrangments

	
	points(cc1[,1],cc1[,2], type = "l", xlim = c(0,max(c1[,1])+1), ylim = c(-1.09,1.09), lwd  = 4, col =  color[pair])
	
	points(cc1[,1],(cc1[,3])/max((c_cor[,7]),na.rm = T),type = "h", col = "#99000025") #recombrate
	#points(c[,1],(c[,10])/max((c[,4]),na.rm = T),type = "h", col = "#99000025") #plot ds
	#points(c[,1],-(c[,10])/max((c[,10]),na.rm = T),type = "h", col = "#00990025") #plot ds

}
axis(1, at = max_position[,3], label = c(1:17))
}

dev.print(device=svg, "results_123/cor_sliding_rearrangments_recombrate2.svg", onefile=FALSE)
dev.off()




