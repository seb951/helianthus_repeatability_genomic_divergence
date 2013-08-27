-#!/usr/bin/Rscript --verbose

###This script does the quantile approach to look at parallel patterns of divergence.
###In addition, it also looks at over-representation of GO terms in the list of the least or most divergence genes.
###





####################
######
### MAPPING FST VALUES ###
##########################
setwd("/home/seb/Documents/repeatability") #set up working directory #set up working directory 

###per gene
comp5_1 =  read.table("results_123/deb_pet_fst_map",  header = T, stringsAsFactors = F)
comp5_2 =  read.table("results_123/bol_exi_fst_map",  header = T, stringsAsFactors = F)
comp5_3 =  read.table("results_123/ann_arg_fst_map",  header = T, stringsAsFactors = F)

###more top3% outliers than expected###

comp5_1_out = comp5_1[comp5_1[,ncol(comp5_1)] == 1, ]
comp5_2_out = comp5_2[comp5_2[,ncol(comp5_2)] == 1, ]
comp5_3_out = comp5_3[comp5_3[,ncol(comp5_3)] == 1, ]

comp5_1_uniq = paste(comp5_1_out[,1],comp5_1_out[,2], sep = "_")
comp5_2_uniq = paste(comp5_2_out[,1],comp5_2_out[,2], sep = "_")
comp5_3_uniq = paste(comp5_3_out[,1],comp5_3_out[,2], sep = "_")

temp = c(comp5_1_uniq,comp5_2_uniq,comp5_3_uniq)
unique_out = cbind(unique(sort(temp)),0,0,0)

for(i in 1:nrow(unique_out))
	{
	x1 = comp5_1_out[comp5_1_uniq == unique_out[i,1],1] 
	x2 = comp5_2_out[comp5_2_uniq == unique_out[i,1],1] 
	x3 = comp5_3_out[comp5_3_uniq == unique_out[i,1],1] 

	if(length(x1) == 1) unique_out[i,2] = 1
	if(length(x2) == 1) unique_out[i,3] = 1
	if(length(x3) == 1) unique_out[i,4] = 1
	}
###
comp5_1_uniq = paste(comp5_1[,1],comp5_1[,2], sep = "_")
comp5_2_uniq = paste(comp5_2[,1],comp5_2[,2], sep = "_")
comp5_3_uniq = paste(comp5_3[,1],comp5_3[,2], sep = "_")

unique_SNP = cbind(as.data.frame(unique(sort(c(comp5_1_uniq,comp5_2_uniq,comp5_3_uniq))), stringsAsFactors = F),NA,NA,NA)

colnames(unique_SNP) = c("reference","fst_deb_pet","fst_bol_exi","fst_ann_arg")



### ### ###
### per SNP ###
### ### ###
unique_SNP = read.table("results_123/unique_SNP", header = T, stringsAsFactors = F)

unique_SNP[is.na(unique_SNP[,2]),2] = 0 #nas back to zeros#
unique_SNP[is.na(unique_SNP[,3]),3] = 0
unique_SNP[is.na(unique_SNP[,4]),4] = 0

top = 0.03; outliers = rep(0,6)

un = unique_SNP[unique_SNP[,2] != 0,]; un = un[un[,3] != 0,];outliers[1] = top *top * nrow(un) #expected polym in both comp...
q1 = quantile(un[,2],1- top,na.rm = T); q2 = quantile(un[,3],1-top,na.rm = T); outliers[2] = length(un[(un[,2] > q1) & (un[,3] > q2),1])#top

un = unique_SNP[unique_SNP[,3] != 0,]; un = un[un[,4] != 0,];outliers[3] = 0.03 *0.03 * nrow(un) #expected polym in both comp...
q1 = quantile(un[,3],1- top,na.rm = T); q3 = quantile(un[,4],1-top,na.rm = T);outliers[4] =  length(un[(un[,3] > q1) & (un[,4] > q3),1])#top

un = unique_SNP[unique_SNP[,2] != 0,]; un = un[un[,4] != 0,];outliers[5] = 0.03 *0.03 * nrow(un) #expected polym in both comp...
q2 = quantile(un[,2],1- top,na.rm = T); q3 = quantile(un[,4],1-top,na.rm = T); outliers[6] = length(un[(un[,2] > q2) & (un[,4] > q3),1])#top

### ### ### ###
### per islands ###
### ### ### ###
c1 =  read.table("results_123/deb_pet_fst_map_cluster",  header = T, stringsAsFactors = F)
c2 =  read.table("results_123/bol_exi_fst_map_cluster",  header = T, stringsAsFactors = F)
c3 =  read.table("results_123/ann_arg_fst_map_cluster",  header = T, stringsAsFactors = F)

top = 0.03; outliers_window = rep(0,6)

q1 = quantile(c1[,3], 1-top)
q2 = quantile(c1[,3],  1-top)
q3 = quantile(c1[,3],  1-top)

outliers_window[2] = length(c1[(c1[,3] > q1) & (c2[,3] > q2),1]) #observed
outliers_window[4] = length(c1[(c2[,3] > q2) & (c3[,3] > q3),1]) #observed
outliers_window[6] = length(c1[(c1[,3] > q1) & (c3[,3] > q3),1]) #observed

outliers_window[1] = top * top * nrow(c1)#expected
outliers_window[3] = top * top  * nrow(c1)#expected
outliers_window[5] = top * top * nrow(c1)#expected

### ### ### ###
### per gene ### 
### ### ### ###
unique_genes = read.delim("results_123/unique_genes", header = T, sep = " ", stringsAsFactors = F)
	
unique_genes[is.na(unique_genes[,2]),2] = 0 
unique_genes[is.na(unique_genes[,3]),3] = 0 
unique_genes[is.na(unique_genes[,4]),4] = 0 

top = 0.03; outliers_genes = rep(0,6)

un = unique_genes[unique_genes[,2] != 0,]; un = un[un[,3] != 0,];outliers_genes[1] = top *top * nrow(un) #expected polym in both comp...
q1 = quantile(un[,2],1- top,na.rm = T); q2 = quantile(un[,3],1-top,na.rm = T); outliers_genes[2] = length(un[(un[,2] > q1) & (un[,3] > q2),1])#top

un = unique_genes[unique_genes[,3] != 0,]; un = un[un[,4] != 0,];outliers_genes[3] = top * top * nrow(un) #expected polym in both comp...
q1 = quantile(un[,3],1- top,na.rm = T); q3 = quantile(un[,4],1-top,na.rm = T);outliers_genes[4] =  length(un[(un[,3] > q1) & (un[,4] > q3),1])#top

un = unique_genes[unique_genes[,2] != 0,]; un = un[un[,4] != 0,];outliers_genes[5] = top * top * nrow(un) #expected polym in both comp...
q2 = quantile(un[,2],1- top,na.rm = T); q3 = quantile(un[,4],1-top,na.rm = T); outliers_genes[6] = length(un[(un[,2] > q2) & (un[,4] > q3),1])#top

### ### ### ### ### ###
### plot outliers OBS EXP ###
### ### ### ### ### ###
plot(c(1,1.3,2,2.3,3,3.3,8,8.3,9,9.3,10,10.3,15,15.3,16,16.3,17,17.3), c(outliers,outliers_genes,outliers_window), type = "h",lwd = 10,xlab = "comparisons",ylab = "number of outliers (top3%)", ylim = c(0,140),col = c(rep("black",2),rep("darkblue",2),rep("darkred",2)), xaxt = "n")
axis(1,c(1,1.3,2,2.3,3,3.3,8,8.3,9,9.3,10,10.3,15,15.3,16,16.3,17,17.3),rep(c("Exp","Obs"),9),las = 2)

text(c(1.15,2.15,3.15),c(128,128,128),srt = 90, label = paste(signif(outliers[c(2,4,6)] / outliers[c(1,3,5)],2),"X",sep = ""),font =2,cex = 1.5 )
text(c(8.15,9.15,10.15),c(75,75,75),srt = 90, label = paste(signif(outliers_genes[c(2,4,6)] / outliers_genes[c(1,3,5)],2),"X",sep = ""),font = 2,cex = 1.5)
text(c(15.15,16.15,17.15),c(75,75,75),srt = 90, label = paste(signif(outliers_window[c(2,4,6)] / outliers_window[c(1,3,5)],2),"X",sep = ""),font = 2,cex = 1.5)

dev.print(device=svg, "results_123/outliers_obsexp.svg",onefile=FALSE)
dev.off()





### ### ### ### ### ###
### BLAST parallel genes ### 
### ### ### ### ### ###
pet_deb_ann_arg = (unique_genes[(unique_genes[,5] == 1) & (unique_genes[,7] == 1),1])
pet_deb_ann_arg = cbind(pet_deb_ann_arg,0,0)

for(seq in 1:nrow(pet_deb_ann_arg))

#for(seq in 1:3)
	{
	grep_cmd = paste("grep '",pet_deb_ann_arg[seq,1],"' -A1 reference/HA412_trinity_noAltSplice_400bpmin.fa >seq_out" ,sep = "")
	system(grep_cmd)
	blast_cmd = system("blastx -evalue 1e-10 -max_target_seqs 5 -query seq_out -db ~/Documents/blast_database/ncbi_protein_nr/nr -num_threads 8 -outfmt 5 -out seq_blast.out")
	system("egrep -i '</BlastOutput>|<Hit_def>|<Hsp_evalue>' seq_blast.out > hit_eval")
	hit_def = as.matrix(read.delim("hit_eval", header = F))
	
		if(length(hit_def) > 1) 
			{
				hit_vector = c(1:nrow(hit_def))[regexpr("<Hit_def>",hit_def, fixed = T) > 0]
	
				hit_only = hit_def[hit_vector,1]

				hit_matrix = cbind(hit_only, hit_def[(hit_vector+1),1],c(1:length(hit_only)))

				hit_matrix[regexpr("unnamed|PREDICTED|hypothetical|predicted|unknown",hit_matrix[,1]) >0,3] = 1000 #discard if possible the hits with unnamed functions.

				best_hit = hit_matrix[as.numeric(hit_matrix[,3]) == min(as.numeric(hit_matrix[,3])),1:2] 
	
				if(length(best_hit) > 2) best_hit = best_hit[1,]
	
				hit = substring(best_hit[1],regexpr(">", best_hit[1],fixed = T)+1,gregexpr("<", best_hit[1],fixed = T)[[1]][2]-1) #best hit
				hit = gsub(",","",hit, fixed = T) # get rid of commas in the annotation.
				eval = as.numeric(substring(best_hit[2],regexpr(">", best_hit[2],fixed = T)+1,gregexpr("<", best_hit[2],fixed = T)[[1]][2]-1)) #best hit evaue
			} 
		else {hit = "no hits"; eval = "1"}

	pet_deb_ann_arg_blast[seq,2:3] = c(hit, eval)
	print(paste(seq,"time is:",Sys.time()))
}
write.table(pet_deb_ann_arg_blast,"pet_deb_ann_arg_parallel.blast", row.names = F, quote = T, col.names = F)


##############SANDBOX

### ### ###
### per SNP ###
### ### ###
unique_SNP = read.table("results_123/unique_SNP", header = T, stringsAsFactors = F)

unique_SNP[is.na(unique_SNP[,2]),2] = 0 #nas back to zeros#
unique_SNP[is.na(unique_SNP[,3]),3] = 0
unique_SNP[is.na(unique_SNP[,4]),4] = 0

outliers = matrix(0,nrow = 20,ncol = 6)
top = seq(0,1,by = 0.05)

for(t in 1:20)
{
un = unique_SNP[unique_SNP[,2] != 0,]; un = un[un[,3] != 0,];outliers[t,1] = top[2] * top[2] * nrow(un) #expected polym in both comp...
q1_a = quantile(un[,2],top[t],na.rm = T);q1_b = quantile(un[,2],top[t+1],na.rm = T)
q2_a = quantile(un[,3],top[t],na.rm = T);q2_b = quantile(un[,3],top[t+1],na.rm = T)
outliers[t,2] = length(un[(un[,2] >= q1_a) & (un[,2] < q1_b) & (un[,3] >= q2_a) & (un[,3] < q2_b)         ,2]) #how many in same quantile in comparison1 and 2? 

un = unique_SNP[unique_SNP[,3] != 0,]; un = un[un[,4] != 0,];outliers[t,3] = top[2] * top[2] * nrow(un) #expected polym in both comp...
q1_a = quantile(un[,3],top[t],na.rm = T);q1_b = quantile(un[,3],top[t+1],na.rm = T)
q2_a = quantile(un[,4],top[t],na.rm = T);q2_b = quantile(un[,4],top[t+1],na.rm = T)
outliers[t,4] = length(un[(un[,3] >= q1_a) & (un[,3] < q1_b) & (un[,4] >= q2_a) & (un[,4] < q2_b)         ,2]) #how many in same quantile in comparison1 and 2? 

un = unique_SNP[unique_SNP[,2] != 0,]; un = un[un[,4] != 0,];outliers[t,5] = top[2] * top[2] * nrow(un) #expected polym in both comp...
q1_a = quantile(un[,2],top[t],na.rm = T);q1_b = quantile(un[,2],top[t+1],na.rm = T)
q2_a = quantile(un[,4],top[t],na.rm = T);q2_b = quantile(un[,4],top[t+1],na.rm = T)
outliers[t,6] = length(un[(un[,2] >= q1_a) & (un[,2] < q1_b) & (un[,4] >= q2_a) & (un[,4] < q2_b)         ,2]) #how many in same quantile in comparison1 and 2? 
}


### ### ###
### per gene ###
### ### ###unique_genes = read.delim("results_123/unique_genes", header = T, sep = " ", stringsAsFactors = F)
unique_genes = read.delim("results_123/unique_genes", header = T, sep = " ", stringsAsFactors = F)

unique_genes[is.na(unique_genes[,2]),2] = 0 
unique_genes[is.na(unique_genes[,3]),3] = 0 
unique_genes[is.na(unique_genes[,4]),4] = 0 

outliers_genes = matrix(0,nrow = 20,ncol = 6)
top = seq(0,1,by = 0.05)

for(t in 1:20)
{
#q1_b = 1/20 * t; q1_a = 1/20 * (t-1) #for fixed quantile of a value of 0.05 fst, from 0 to 1.
#q2_b = 1/20 * t; q2_a = 1/20 * (t-1) #for fixed quantile of a value of 0.05 fst, from 0 to 1.
un = unique_genes[unique_genes[,2] != 0,]; un = un[un[,3] != 0,];outliers_genes[t,1] = top[2] * top[2] * nrow(un) #expected polym in both comp...
#un = unique_genes[unique_genes[,2] != 0,]; un = un[un[,3] != 0,];outliers_genes[t,1] = (length(un[(un[,2] >= q1_a) & (un[,2] < q1_b),2]) *  length(un[(un[,3] >= q2_a) & (un[,3] < q2_b),3]) )  / nrow(un) #expected polym in both comp for fixed quantile of a value of 0.05 fst, from 0 to 1.

q1_a = quantile(un[,2],top[t],na.rm = T);q1_b = quantile(un[,2],top[t+1],na.rm = T)
q2_a = quantile(un[,3],top[t],na.rm = T);q2_b = quantile(un[,3],top[t+1],na.rm = T)
outliers_genes[t,2] = length(un[(un[,2] >= q1_a) & (un[,2] < q1_b) & (un[,3] >= q2_a) & (un[,3] < q2_b),2]) #how many in same quantile in comparison1 and 2? 
un = unique_genes[unique_genes[,3] != 0,]; un = un[un[,4] != 0,];outliers_genes[t,3] = top[2] * top[2] * nrow(un) #expected polym in both comp...
#un = unique_genes[unique_genes[,3] != 0,]; un = un[un[,4] != 0,];outliers_genes[t,3] = (length(un[(un[,4] >= q1_a) & (un[,4] < q1_b),4]) *  length(un[(un[,3] >= q2_a) & (un[,3] < q2_b),3]) )  /nrow(un) #expected polym in both comp...

q1_a = quantile(un[,3],top[t],na.rm = T);q1_b = quantile(un[,3],top[t+1],na.rm = T)
q2_a = quantile(un[,4],top[t],na.rm = T);q2_b = quantile(un[,4],top[t+1],na.rm = T)
outliers_genes[t,4] = length(un[(un[,3] >= q1_a) & (un[,3] < q1_b) & (un[,4] >= q2_a) & (un[,4] < q2_b),2]) #how many in same quantile in comparison1 and 2? 
un = unique_genes[unique_genes[,2] != 0,]; un = un[un[,4] != 0,];outliers_genes[t,5] = top[2] * top[2] * nrow(un) #expected polym in both comp...
#un = unique_genes[unique_genes[,2] != 0,]; un = un[un[,4] != 0,];outliers_genes[t,5] = (length(un[(un[,2] >= q1_a) & (un[,2] < q1_b),2]) *  length(un[(un[,4] >= q2_a) & (un[,4] < q2_b),4]) )  /nrow(un) 
#expected polym in both comp...
q1_a = quantile(un[,2],top[t],na.rm = T);q1_b = quantile(un[,2],top[t+1],na.rm = T)
q2_a = quantile(un[,4],top[t],na.rm = T);q2_b = quantile(un[,4],top[t+1],na.rm = T)
outliers_genes[t,6] = length(un[(un[,2] >= q1_a) & (un[,2] < q1_b) & (un[,4] >= q2_a) & (un[,4] < q2_b),2]) #how many in same quantile in comparison1 and 2? 
}

### ### ### ###
###per islands###
### ### ### ###
c1 =  read.table("results_123/deb_pet_fst_map_cluster",  header = T, stringsAsFactors = F)
c2 =  read.table("results_123/bol_exi_fst_map_cluster",  header = T, stringsAsFactors = F)
c3 =  read.table("results_123/ann_arg_fst_map_cluster",  header = T, stringsAsFactors = F)

c1_2_3 = cbind(c1[,c(1,3)],c2[,3],c3[,3])

outliers_window = matrix(0,nrow = 20,ncol = 6)
top = seq(0,1,by = 0.05)

for(t in 1:20)
{
un = c1_2_3[c1_2_3[,2] != 0,]; un = un[un[,3] != 0,];outliers_window[t,1] = top[2] * top[2] * nrow(un) #expected polym in both comp...
q1_a = quantile(un[,2],top[t],na.rm = T);q1_b = quantile(un[,2],top[t+1],na.rm = T)
q2_a = quantile(un[,3],top[t],na.rm = T);q2_b = quantile(un[,3],top[t+1],na.rm = T)
outliers_window[t,2] = length(un[(un[,2] >= q1_a) & (un[,2] < q1_b) & (un[,3] >= q2_a) & (un[,3] < q2_b),2]) #how many in same quantile in comparison1 and 2? 

un = c1_2_3[c1_2_3[,3] != 0,]; un = un[un[,4] != 0,];outliers_window[t,3] = top[2] * top[2] * nrow(un) #expected polym in both comp...
q1_a = quantile(un[,3],top[t],na.rm = T);q1_b = quantile(un[,3],top[t+1],na.rm = T)
q2_a = quantile(un[,4],top[t],na.rm = T);q2_b = quantile(un[,4],top[t+1],na.rm = T)
outliers_window[t,4] = length(un[(un[,3] >= q1_a) & (un[,3] < q1_b) & (un[,4] >= q2_a) & (un[,4] < q2_b),2]) #how many in same quantile in comparison1 and 2? 

un = c1_2_3[c1_2_3[,2] != 0,]; un = un[un[,4] != 0,];outliers_window[t,5] = top[2] * top[2] * nrow(un) #expected polym in both comp...
q1_a = quantile(un[,2],top[t],na.rm = T);q1_b = quantile(un[,2],top[t+1],na.rm = T)
q2_a = quantile(un[,4],top[t],na.rm = T);q2_b = quantile(un[,4],top[t+1],na.rm = T)
outliers_window[t,6] = length(un[(un[,2] >= q1_a) & (un[,2] < q1_b) & (un[,4] >= q2_a) & (un[,4] < q2_b),2]) #how many in same quantile in comparison1 and 2? 
}


### per gene
### PLOTTING per gene better version ###
### per gene
par(mar = c(6,6,2,2))  
plot(c(2:21),outliers_genes[,6]/outliers_genes[,5],col = "darkblue", ylim= c(0,4.5),ylab = "", xlab = "", xaxt = "n",yaxt = "n",lwd = 3,font = 2)
#plot(c(2:21)[outliers_genes[,6] != 0],outliers_genes[outliers_genes[,6] != 0,6]/outliers_genes[outliers_genes[,6] != 0,5],xlim = c(2,21),col = "darkblue", ylim= c(0,11),ylab = "", xlab = "", xaxt = "n",yaxt = "n",lwd = 3,font = 2)
points(c(2:21)[outliers_genes[,4] != 0],outliers_genes[outliers_genes[,4] != 0,4]/outliers_genes[outliers_genes[,4] != 0,3],lwd = 3, pch = 8,col = "darkgreen")
points(c(2:21)[outliers_genes[,2] != 0],outliers_genes[outliers_genes[,2] != 0,2]/outliers_genes[outliers_genes[,2] != 0,1],lwd = 3, pch = 2,col = "#00000099")
lm_xxy = lm((outliers_genes[,6]/outliers_genes[,5])~c(2:21)+I(c(2:21)^2));qftn <- function(x) lm_xxy$coefficients[1] + lm_xxy$coefficients[2]*x + lm_xxy$coefficients[3]*x^2;curve(qftn, 2, 21,col = "darkblue",,add = T, lwd = 3)

axis(1,at = c(1:21)+0.5, labels = paste( seq(0,100,by = 5), "%", sep = ""), las = 2, pos = -0.3)
#axis(1,at = c(1:21)+0.5, labels = seq(0,1,by = 0.05), las = 2, pos = -0.3)
#axis(1,at = c(1:21)+0.5, labels = seq(0,1,by = 0.05), las = 2, pos = -0.6)
axis(1, at = 12,labels = "Quantiles",font.axis = 2, line = 3.5,cex.axis = 2, tick = F)
axis(2,at = c(0,1,2,3,4), labels = c(0,1,2,3,4), las = 2, pos = 0.9,las = 2)
#axis(2,at = c(0,1,2,3,4,5,6,7,8,9,10,11), labels = c(0,1,2,3,4,5,6,7,8,9,10,11), las = 2, pos = 0.9,las = 2)
axis(2,at  = 2, labels =  "Observed / Expected",font.axis = 2, line = 3,cex.axis = 2, tick = F)
axis(2, at = 2,labels =  "shared genes",font.axis = 2, line = 1.4,cex.axis = 2, tick = F)
#axis(2,at  = 5, labels =  "Observed / Expected",font.axis = 2, line = 3,cex.axis = 2, tick = F)
#axis(2, at = 5,labels =  "shared genes",font.axis = 2, line = 1.4,cex.axis = 2, tick = F)

main = c("H. petiolaris - H. debilis versus H. bolanderi - H. exilis","H. bolanderi - H. exilis versus H. annuus - H. argophyllus","H. annuus - H. argophyllus versus H. petiolaris - H. debilis")
legend(y = 4.5,x = 4.5,legend = main, fill = c("#00000099","darkgreen","darkblue"), cex = 0.8,text.font = 4, )
#legend(y = 11,x = 3.5,legend = main, fill = c("#00000099","darkgreen","darkblue"), cex = 0.8,text.font = 4, )
dev.print(device=pdf, "figures_2013-03-07/4.quadratic_distrib.pdf",onefile=FALSE)
#dev.print(device=pdf, "figures_2013-03-07/S1.quadratic_distrib_fixed_quantiles.pdf",onefile=FALSE)
dev.off()



###
###do some GO analysis to find out what is going on with the most and least divergent genes
###
setwd("/home/seb/Documents/repeatability") #set up working directory #set up working directory 

all_genes = read.table("~/Documents/speciation_islands_dec2011/reference/go_annotations/out3",stringsAsFactors = F, sep = "\t",skip = 1)
outliers_genes = matrix(0,nrow = 20,ncol = 6)
expression = cbind(all_genes[,1],0)
top = seq(0,1,by = 0.05)

for(i in 1:nrow(unique_genes))
	{	
	x =  expression[,1] %in% unique_genes[i,1]
	if(length(x[x==T]) ==1) expression[x ==T,2] = 0.9 #expressed genes, the non expressed stay at zero.
	if(i %% 1000 == 0) print(paste(Sys.time(),i))
	}

a = 2; b = 4#which comparison are you looking at

for(t in 1:20)
	{
	un = unique_genes[unique_genes[,a] != 0,]; un = un[un[,b] != 0,];outliers_genes[t,5] = top[2] * top[2] * nrow(un) #expected polym in both comp...
	q1_a = quantile(un[,a],top[t],na.rm = T);q1_b = quantile(un[,a],top[t+1],na.rm = T)
	q2_a = quantile(un[,b],top[t],na.rm = T);q2_b = quantile(un[,b],top[t+1],na.rm = T)
	temp = un[(un[,a] >= q1_a) & (un[,a] < q1_b) & (un[,b] >= q2_a) & (un[,b] < q2_b),1] #how many in same quantile in comparison1 and 2? 

	for(ti in 1:length(temp)) {
	x =  expression[,1] %in% temp[ti]
	if(length(x[x==T]) ==1) expression[x ==T,2] = t / 200} #top 5% most divergent get a t/200 "pvalue"
	}
	
expression2 = expression[all_genes[,2] != "no_hits",] #remove genes with no GO hits
expression3 = expression2[expression2[,2] != "0",] #remove unexpressed genes.

write.table(expression3,"~/Documents/speciation_islands_dec2011/reference/go_annotations/expressed_genes.txt",row.names = F, col.names = F, quote = F, sep = "\t")
#write.table(all_genes[all_genes[,2] != "no_hits",],"~/Documents/speciation_islands_dec2011/reference/go_annotations/out4",row.names = F, col.names = F, quote = F, sep = "\t") #must add a white space in first column.

system("echo >~/Documents/speciation_islands_dec2011/reference/go_annotations/out5")
system("cat ~/Documents/speciation_islands_dec2011/reference/go_annotations/out4 >>~/Documents/speciation_islands_dec2011/reference/go_annotations/out5")
#system("awk '{gsub(/NA/,"")}; 1' out5 >out4")
#system("rm out5")


write.table(unique_genes[!is.na(unique_genes[,4]),c(1,4)],"~/Documents/speciation_islands_dec2011/reference/go_annotations/expressed_genes.txt",row.names = F, col.names = F, quote = F, sep = "\t")
write.table(all_genes[all_genes[,2] != "no_hits",],"~/Documents/speciation_islands_dec2011/reference/go_annotations/out4",row.names = F, col.names = F, quote = F, sep = "\t") #must add a white space in first column.

