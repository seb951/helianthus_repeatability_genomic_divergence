setwd("/home/seb/Documents/repeatability") #set up working directory #set up working directory 

###This script looks at the correlation coefficient for genomic windows using a permuted random dataset. This is to make sure that there is no effect of average many SNP over few genomic regions. 


####prep the map
map = read.delim("reference/ordered_transcriptome_trinity.txt", header = F, stringsAsFactors = F) # new map
map = map[map[,3] != "no",]
map = cbind(map,0,0,0,0)
colnames(map) = c("name","LG","map_centiMorgan","unique_position","nbnuc","dn","ds")

ref = read.delim("reference/HA412_trinity_noAltSplice_400bpmin.fa", header = F, stringsAsFactors = F) # new reference assembly
len = cbind(gsub(">","",ref[seq(from = 1, to = nrow(ref), by = 2),]),apply(ref, 1,nchar)[seq(from = 2, to = nrow(ref), by = 2)])

#unique position for the markers#
max_position = cbind(c(1:17), 0,0) # the size of each chromosome and the cumulative number of centimorgans.

for(i in 1:17)
	{
	max_position[i,2] = max(as.numeric(map[as.numeric(map[,2]) == i,3]))
	max_position[i,3] = sum(max_position[1:i,2])-max_position[i,2]
	}

for(i in 1:nrow(map))
#for(i in 1:200)
	{
	for(j in 1:17)
		{
		if(as.numeric(map[i,2]) == max_position[j,1]) map[i,4] = (as.numeric(map[i,3]) + sum(max_position[1:j,2]) - sum(max_position[j,2]) + 1)
		}
		map[i,5] = len[map[i,1] == len[,1],2]
	if(i %% 1000 == 0) print(paste(i,Sys.time()))
}
max_position[,3] = (max_position[,3] +1)

##########################
### MAPPING FST VALUES at random! ###
##########################

perm_res = matrix(0,nrow = 10, ncol = 3)

for(perm in 1:10)
{

for(pair in 1:3)
{

if(pair == 1) fst_map_match = read.delim("results_123/deb_pet_fst_map", header = T, sep = " ", stringsAsFactors = F) # ann_pet_fst
if(pair == 2) fst_map_match = read.delim("results_123/bol_exi_fst_map", header = T, sep = " ", stringsAsFactors = F) # pet_deb_fst
if(pair == 3) fst_map_match = read.delim("results_123/ann_arg_fst_map", header = T, sep = " ", stringsAsFactors = F) # ann_deb_fst

gmean=function(y=y,x=x, q=q,d=d){  # y=values, x=snp pos vector, q=query pos, d=dist:    Rose's gaussian average
	weights=rep(NA,length(x))
	weights[abs(x-q)<3*d]=exp(-(x[abs(x-q)<3*d]-q)^2/(2*d^2))
	weighted=weights*y/sum(weights,na.rm=T)
	return(sum(weighted,na.rm=T))
}

	col = ncol(fst_map_match) 
	fst_map_match[is.na(fst_map_match[,(col-8)]),(col-8):(col-6)] = 0
	#fst_map_match[,(col-8)] = jitter(fst_map_match[,(col-8)],0.1)  #add/remove between -/+ 2e-05. to all fst values 	
	top_q = quantile(fst_map_match[,(col-8)],0.97) #top3%
	
	####running window - distance based###
 	d = 0.001 #distance cutoff in cM
	
	c = cbind(unique(as.numeric(map[,4])),0,0,0,0,0,0,0,0) #position you are interrogating along with the result of the distance based running mean. 
	
	for(i in 1:nrow(c)){temp = map[as.numeric(map[,4]) == c[i,1],];c[i,2] = temp[,2][1]}

	c = matrix(as.numeric(c), ncol = 9)
	colnames(c) = c("unique_position","LG","mean_sliding_fst", "top3_window","bottom97_window", "number_genes_in_SW","resampling_mean","resampling_top3","nbnuc")

	for(chr in 1:17)
		{
		chromo_temp = fst_map_match[(as.numeric(fst_map_match[,col-4]) == chr),] # subset of chromosome a
	
		for(i in c(1:nrow(c))[c[,2] == chr]) 
			{
			temp = chromo_temp[chromo_temp[,(col-2)] > (c[i,1]-d) & chromo_temp[,(col-2)] < (c[i,1]+d),] # mean fst at unique position in a sliding window of X cM
				#as.numeric(sample(temp[,5],6))
			#if(nrow(temp) >0) c[i,3] = gmean(temp[,(col-8)], temp[,(col-2)],c[i,1],d)	 #gaussian mean
			if(nrow(temp) > 0)	c[i,3] = mean(as.numeric(sample(fst_map_match[,col-8],nrow(temp))))
			}

		#print(paste(chr, Sys.time()))
}

if(pair == 1) write.table(c,"results_permutted/deb_pet_fst_map_cluster_nonoverlapping_windows", row.names = F)
if(pair == 2) write.table(c,"results_permutted/bol_exi_fst_map_cluster_nonoverlapping_windows", row.names = F)
if(pair == 3) write.table(c,"results_permutted/ann_arg_fst_map_cluster_nonoverlapping_windows", row.names = F)

}
################
### PLOTTING ###
################
c1 =  read.table("results_permutted/deb_pet_fst_map_cluster_nonoverlapping_windows",  header = T, stringsAsFactors = F);c1[c1[,3] == 0,3] = NA #non overlapping windows (to avoid autocorrelation problems...)
c2 =  read.table("results_permutted/bol_exi_fst_map_cluster_nonoverlapping_windows",  header = T, stringsAsFactors = F);c2[c2[,3] == 0,3] = NA
c3 =  read.table("results_permutted/ann_arg_fst_map_cluster_nonoverlapping_windows",  header = T, stringsAsFactors = F);c3[c3[,3] == 0,3] = NA

#############
###per region ###
############
c12 = cbind(c1[,c(1,3)],c2[,3])
c23 = cbind(c2[,c(1,3)],c3[,3])
c13 = cbind(c1[,c(1,3)],c3[,3])

perm_res[perm,1] = signif(cor.test(c12[,2],c12[,3], alternative = "g")$estimate,2)
perm_res[perm,2] = signif(cor.test(c23[,2],c23[,3], alternative = "g")$estimate,2)
perm_res[perm,3] = signif(cor.test(c13[,2],c13[,3], alternative = "g")$estimate,2)

print(paste(perm,Sys.time()))
}


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




