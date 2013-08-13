#!/usr/bin/Rscript --verbose
args = commandArgs(TRUE)
pair = as.numeric(args[1])


###This script parses the raw SNP table produced by mpileup and bcftools to only keep a small fraction of high quality SNPs. 
#It also calculates Fst per SNP and per genomic position.



setwd("/home/seb/Documents/repeatability") #set up working directory 

###kick out SNP which have more than 10% XX data missing. Allow more for the annuus comparisons given that there is more 454 libraries in annuus which contain much more missing data.strrrr
individuals = read.delim("reference/all_species_mar2013.txt", header = T, stringsAsFactors = F) #individuals of interest
individuals = individuals[-c(11,12,14,16),] #Ames7109_bol, BOL2436_bol,BOL1024_bol,BOL1023_bol are removed according to the splistree result


if(pair == 1) {comp = read.table("mpileup_123/snp_table_1", header = T, stringsAsFactors = F);individuals = individuals[regexpr("deb|pet",individuals[,2])>0,]}
if(pair == 2) {comp = read.table("mpileup_123/snp_table_2", header = T, stringsAsFactors = F);individuals =individuals[regexpr("EAST|WEST",individuals[,2])>0,]} #according to the splistree, ind were labelled as east/west. Some were removed. 
if(pair == 3) {comp = read.table("mpileup_123/snp_table_3", header = T, stringsAsFactors = F);individuals =individuals[regexpr("ann|arg",individuals[,2])>0,]}

#comp = comp[c(1:10000),]
colnames(comp)[4:ncol(comp)] = paste(colnames(comp)[4:ncol(comp)],individuals[,2],sep = "_")

too_much_missing = c(1:nrow(comp))


	for(i in 1:nrow(comp))
		{
		n_ind = (ncol(comp)-3)
		too_much_missing[i] = length(c(1:n_ind)[grepl("XX",comp[i,4:ncol(comp)]) == T]) / n_ind
		if(i %% 10000 == 0) print(paste(i,Sys.time()))
		}

if(pair == 1) comp2 = comp[too_much_missing < 0.2,]
if(pair == 2) comp2 = comp[too_much_missing < 0.3,]
if(pair == 3) comp2 = comp[too_much_missing < 0.2,]



##############################
### trim based on 2pq (Ht) ###
##############################
comp3 =  c(1)

#counting the alleles.
	counter_ACGTX = matrix(0,nrow = nrow(comp2), ncol = 5)
	colnames(counter_ACGTX) = c("A","C","G","T","X")
	ht = c(1:nrow(comp2))
	n_ind = ncol(comp2)

	for(i in 1:nrow(comp2))
		{
		x = paste(comp2[i,4:n_ind], collapse = "")
		xx = strsplit(x, split = "")
		
		counter_ACGTX[i,1] = length(c(1:((n_ind-3)*2))[grepl("A", xx[[1]])])
		counter_ACGTX[i,2] = length(c(1:((n_ind-3)*2))[grepl("C", xx[[1]])])
		counter_ACGTX[i,3] = length(c(1:((n_ind-3)*2))[grepl("G", xx[[1]])])
		counter_ACGTX[i,4] = length(c(1:((n_ind-3)*2))[grepl("T", xx[[1]])])
		counter_ACGTX[i,5] = length(c(1:((n_ind-3)*2))[grepl("X", xx[[1]])])
		p = sort(counter_ACGTX[i,1:4])[4] 
		q = sort(counter_ACGTX[i,1:4])[3] 
		al = sum(counter_ACGTX[i,1:4])
		ht[i] = 2 * (p / al) * (q / al)
		}
	comp3 = comp2[ht > 0.2,]

#comp2 = comp = NULL

###################################
### trim based on Ho (paralogs) ###
###################################
comp4 = c(1)
ho  = rep(0,nrow(comp3))

for(i in 1:nrow(comp3))
	{
		a1 = substring(comp3[i,4:ncol(comp3)],1,1)
		a2 = substring(comp3[i,4:ncol(comp3)],2,2)		
		for(h in 1:length(a1))
		{
		if((a1[h] != a2[h]) & (a1[h] != "X") & (a2[h] != "X") & !is.na(a1[h]) & !is.na(a2[h])) ho[i] = (ho[i] + 1) #count the heterozygotes
		}
	}
	ho = ho / (ncol(comp3)-3) # observed heterozygosity
	comp4 = comp3[ho < 0.6,]


########################
### hierfstat format ###
########################
#takes about 1min per 1000 snps. 
library(hierfstat)

fstat_results = cbind(c(1:nrow(comp4)),0,0,0)

colnames(fstat_results) = c("snp","fst","fit","fis")

	for(i in 1:nrow(fstat_results))
	{
	hf = t( rbind( rep("empty",length( comp4[i,4:ncol(comp4)])),comp4[i,4:ncol(comp4)]))

	if(pair == 1) { hf[regexpr("_deb",rownames(hf)) >0,1] = "1";  hf[regexpr("_pet",rownames(hf)) >0,1] = "2"}
	if(pair == 2) { hf[regexpr("_EAST",rownames(hf)) >0,1] = "1";  hf[regexpr("_WEST",rownames(hf)) >0,1] = "2"}
	if(pair == 3) { hf[regexpr("_ann",rownames(hf)) >0,1] = "1"; hf[regexpr("_arg",rownames(hf)) >0,1] = "2"}	
	hf  = gsub("U|R|Y|M|K|S|W|B|D|H|V|N","X",hf)
	hf = gsub("X","",hf)
	hf = gsub("0","",hf)
	hf = gsub("A",1, hf)
	hf = gsub("C",2, hf)
	hf = gsub("G",3, hf)
	hf = gsub("T",4, hf)

	colnames(hf) = c("pop","snp1")

	if(length(hf[hf[,2] == "",2])/length(hf[,2]) < 0.5) x = varcomp(data.frame(as.integer(hf[,1]), as.integer(hf[,2])))
	fstat_results[i,2] = x$F[1,1] #Fst
	fstat_results[i,3] = x$F[1,2] #Fit
	fstat_results[i,4] = x$F[2,2] #Fis
	if(i %% 1000 == 0) print(paste(i,"of",nrow(fstat_results),Sys.time()))
	}

comp4 = cbind(comp4, fstat_results)
####################
### MAPPING FST VALUES ###
##########################

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


###keep info seperate per snp.
	map = map[,1:5]

	comp5 =  cbind(comp4,0,0,0,0,0)
	colnames(comp5) = c(colnames(comp4),colnames(map))
	map_temp = matrix(0,nrow = 100, ncol = 5)
	ceiling = ceiling(nrow(comp5)/100)
	map_all = NULL
	for(j in 1:ceiling)
		{
			if(j != ceiling) map_temp = matrix(0,nrow = 100, ncol = 5) else map_temp = matrix(0, nrow = nrow(comp5) %% 100, ncol = 5)
			if(length(map_temp) == 0) map_temp  = matrix(0,nrow = 100, ncol = 5)
			for(z in 1:nrow(map_temp))
			{
			x = map[,1] %in% comp5[(((j -1)*100)+ z),1]
			if(length(x[x == T]) == 1)  map_temp[z,] = as.matrix(map[(x == T),])
			if(length(x[x == T]) > 1) print(map[(x == T),])
			}			
			map_all = rbind(map_all,map_temp)
			if(j %% 100 == 0) print(paste(j,"of ceiling",ceiling,Sys.time()))		
			}		
			comp5[,(ncol(comp4)+1):(ncol(comp4)+5)] = map_all
	

###add a column to flag whether the SNP is in the top 3% quantile.###

	col = ncol(comp5) 
	top_quantile = quantile(as.numeric(comp5[,col-7]), 0.97,na.rm = T)
	comp5 = cbind(comp5,0)
	colnames(comp5)[col+1] = "top-quant"
	for(i in 1:nrow(comp5))
		{if((as.numeric(comp5[i,col-7]) > top_quantile)   & !is.na(comp5[i,col-7])  ) comp5[i,col +1] = 1}

comp5[,(col-7)] = jitter(comp5[,(col-7)],0.1)  #add/remove between -/+ 2e-05. to all fst values 


cc = Sys.time()

if(pair == 1) write.table(comp5,"results_123/deb_pet_fst_map", row.names = F)
if(pair == 2) write.table(comp5,"results_123/bol_exi_fst_map", row.names = F)
if(pair == 3) write.table(comp5,"results_123/ann_arg_fst_map", row.names = F)
	
###################
#### CLUSTERING ###
###################
gmean=function(y=y,x=x, q=q,d=d){  # y=values, x=snp pos vector, q=query pos, d=dist:    Rose's gaussian average
	weights=rep(NA,length(x))
	weights[abs(x-q)<3*d]=exp(-(x[abs(x-q)<3*d]-q)^2/(2*d^2))
	weighted=weights*y/sum(weights,na.rm=T)
	return(sum(weighted,na.rm=T))
}


if(pair == 1) fst_map_match = read.delim("results_123/deb_pet_fst_map", header = T, sep = " ", stringsAsFactors = F) # ann_pet_fst
if(pair == 2) fst_map_match = read.delim("results_123/bol_exi_fst_map", header = T, sep = " ", stringsAsFactors = F) # pet_deb_fst
if(pair == 3) fst_map_match = read.delim("results_123/ann_arg_fst_map", header = T, sep = " ", stringsAsFactors = F) # ann_deb_fst

	col = ncol(fst_map_match) 
	fst_map_match[is.na(fst_map_match[,(col-8)]),(col-8):(col-6)] = 0
	#fst_map_match[,(col-8)] = jitter(fst_map_match[,(col-8)],0.1)  #add/remove between -/+ 2e-05. to all fst values 	
	top_q = quantile(fst_map_match[,(col-8)],0.97) #top3%
	
	####running window - distance based###
 	d = 0.001 #distance cutoff in cM
	
	c = cbind(unique(as.numeric(map[,4])),0,0,0,0,0,0,0,0) #position you are interrogating along with the result of the distance based running mean. 
	
	for(i in 1:nrow(c))
		{		temp = map[as.numeric(map[,4]) == c[i,1],]
				c[i,2] = temp[,2][1]
				c[i,9] =  sum(as.numeric(temp[,5]))

	}

	c = matrix(as.numeric(c), ncol = 9)
	colnames(c) = c("unique_position","LG","mean_sliding_fst", "top3_window","bottom97_window", "number_genes_in_SW","resampling_mean","resampling_top3","nbnuc")
#	colnames(c) = c("unique_position","LG","mean_sliding_fst", "top3_window","bottom97_window", "number_SNPs_in_SW","resampling_mean","resampling_top3","nbnuc")

	for(chr in 1:17)
		{
		chromo_temp = fst_map_match[(as.numeric(fst_map_match[,col-4]) == chr),] # subset of chromosome a
	
		for(i in c(1:nrow(c))[c[,2] == chr]) 
			{
			temp = chromo_temp[chromo_temp[,(col-2)] > (c[i,1]-d) & chromo_temp[,(col-2)] < (c[i,1]+d),] # mean fst at unique position in a sliding window of X cM
			
			if(nrow(temp) >0) c[i,3] = gmean(temp[,(col-8)], temp[,(col-2)],c[i,1],d)	 #gaussian mean
			if(nrow(temp) >0) c[i,4] = length(temp[temp[,(col-8)] > top_q,(col-8)]) # observed fixed SNP per sliding window...
			if(nrow(temp) >0) c[i,5] = length(temp[temp[,(col-8)] < top_q,(col-8)])  # observed non fixed SNP per sliding window...
			if(nrow(temp) >0) c[i,6] = length(unique(temp[,1])) # to get number of genes per region
#			if(nrow(temp) >0) c[i,6] = nrow(temp) # to get number of snp per region

			
			perm_100 = c(0,0) # mean and fixed resamplers.
			for(s in 1:10000) # 1000 resampling
				{			
				temp_perm = sample(fst_map_match[,(col-8)],nrow(temp))
				if((nrow(temp) >0) &  (c[i,3] > mean(temp_perm) ))	 perm_100[1] = perm_100[1] + 1
				if((nrow(temp) >0) & (c[i,4] > length(temp_perm[temp_perm >top_q ])))	 perm_100[2] = perm_100[2] + 1
				}
			
			if(perm_100[1] > 900 | perm_100[2] > 900) 	for(s in 1:99000) # 10000 resampling
				{			
				temp_perm = sample(fst_map_match[,(col-8)],length(temp))
				if(c[i,3] >	mean(temp_perm))	 perm_100[1] = perm_100[1] + 1
				if(c[i,4] > length(temp_perm[temp_perm > top_q ]))	 perm_100[2] = perm_100[2] + 1
				} 

			if(perm_100[1] < 1000)   c[i,7] = perm_100[1] / 1000 else c[i,7] = perm_100[1] / 100000
			if(perm_100[2] < 1000)   c[i,8] = perm_100[2] / 1000 else c[i,8] = perm_100[2] / 100000
			perm_100 = NULL
			}

		print(paste(chr, Sys.time()))
}

bb = Sys.time()
if(pair == 1) write.table(c,"results_123/deb_pet_fst_map_cluster_nonoverlapping_windows", row.names = F)
if(pair == 2) write.table(c,"results_123/bol_exi_fst_map_cluster_nonoverlapping_windows", row.names = F)
if(pair == 3) write.table(c,"results_123/ann_arg_fst_map_cluster_nonoverlapping_windows", row.names = F)




