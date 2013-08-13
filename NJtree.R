#!/usr/bin/Rscript --slave

###This script is used to take a subset of SNP to generate a file ready to be loaded into SPLISTREE.

setwd("/home/seb/Documents/repeatability") #set up working directory 

###############################
#### NEIGHBOUR JOINING TREE ###
###############################
snp_table = as.matrix(read.delim("results_6species/snp_table_all", header = T, sep = " ", stringsAsFactors = F)) # this is a subset of about 120 00 snps for the first XYZ genes###
###snp_table = as.matrix(read.delim("results_6species/snp_table_bol_exi", header = T, sep = " ", stringsAsFactors = F)) 
snp_table[,2] = gsub(" ","",snp_table[,2])
individuals = read.delim("reference/all_species_nov2012_cleaned.txt", header = T, stringsAsFactors = F) #individuals of interest (454 individuals removed, weird individuals based on NatCom paper also removed#
###individuals = individuals[regexpr("exi|bol",individuals[,2])>0,] #GregO only wants exi and bolanderi

colnames(snp_table) = gsub("X14","14",colnames(snp_table)); colnames(snp_table) = gsub("X2O","2O",colnames(snp_table))
#colnames(snp_table)[4:ncol(snp_table)] = paste(colnames(snp_table)[4:ncol(snp_table)],individuals[,2],sep = "_") #dont need this if you are using a single species pair. 

################
### TRIMMING ###
################
all_comparisons = snp_table
###kick out SNP which have more than 20% XX data missing
too_much_missing = c(1:nrow(all_comparisons))

	for(i in 1:nrow(all_comparisons))
		{
		n_ind = ncol(all_comparisons)-3
		too_much_missing[i] = length(c(1:n_ind)[grepl("XX",all_comparisons[i,4:ncol(all_comparisons)]) == T]) / n_ind
		}

all_comparisons_2 = all_comparisons[too_much_missing < 0.2,]

##############################
### trim based on 2pq (Ht) ###
##############################
all_comparisons_3 = 1

#counting the alleles.
	counter_ACGTX = matrix(0,nrow = nrow(all_comparisons_2), ncol = 5)
	colnames(counter_ACGTX) = c("A","C","G","T","X")
	ht = c(1:nrow(all_comparisons_2))
	n_ind = ncol(all_comparisons_2)
	
	for(i in 1:nrow(all_comparisons_2))
		{
		x = paste(all_comparisons_2[i,4:n_ind], collapse = "")
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
	all_comparisons_3 = all_comparisons_2[ht > 0.2,]


###################################
### trim based on Ho (paralogs) ###
###################################

all_comparisons_4 = list(1)

ho  = rep(0,nrow(all_comparisons_3))

for(i in 1:nrow(all_comparisons_3))
	{
		a1 = substring(all_comparisons_3[i,4:ncol(all_comparisons_3)],1,1)
		a2 = substring(all_comparisons_3[i,4:ncol(all_comparisons_3)],2,2)		
		for(h in 1:length(a1))
		{
		if((a1[h] != a2[h]) & (a1[h] != "X") & (a2[h] != "X") & !is.na(a1[h]) & !is.na(a2[h])) ho[i] = (ho[i] + 1) #count the heterozygotes
		}
	}
	ho = ho / (ncol(all_comparisons_3)-3) # observed heterozygosity
	all_comparisons_4 = all_comparisons_3[ho < 0.6,]

###################################
### create a single consensus of both alleles. ###
###################################	
	
col = ncol(all_comparisons_4)
all_comparisons_4_nj_preformat = all_comparisons_4

for(i in 1:nrow(all_comparisons_4))
	{
	temp1 = all_comparisons_4[i,4:col]
	temp2 = temp1
	for(j in 1:length(temp2))
		{
		if(substring(temp1[j],1,1) == substring(temp1[j],2,2)) temp2[j] = substring(temp1[j],1,1)
		if(substring(temp1[j],1,1) == "X" & substring(temp1[j],2,2) != "X" ) temp2[j] = substring(temp1[j],2,2)
		if(substring(temp1[j],2,2) == "X" & substring(temp1[j],1,1) != "X" ) temp2[j] = substring(temp1[j],1,1)
	
		if((length(grep("X",temp1[j])) < 1) & ((substring(temp1[j],1,1) == "A" & substring(temp1[j],2,2) == "C") | (substring(temp1[j],1,1) == "C" & substring(temp1[j],2,2) == "A"))) temp2[j] = "M"
		if((length(grep("X",temp1[j])) < 1) & ((substring(temp1[j],1,1) == "A" & substring(temp1[j],2,2) == "G") | (substring(temp1[j],1,1) == "G" & substring(temp1[j],2,2) == "A"))) temp2[j] = "T"
		if((length(grep("X",temp1[j])) < 1) & ((substring(temp1[j],1,1) == "A" & substring(temp1[j],2,2) == "T") | (substring(temp1[j],1,1) == "T" & substring(temp1[j],2,2) == "A"))) temp2[j] = "W"
		if((length(grep("X",temp1[j])) < 1) & ((substring(temp1[j],1,1) == "C" & substring(temp1[j],2,2) == "G") | (substring(temp1[j],1,1) == "G" & substring(temp1[j],2,2) == "C"))) temp2[j] = "S"
		if((length(grep("X",temp1[j])) < 1) & ((substring(temp1[j],1,1) == "C" & substring(temp1[j],2,2) == "T") | (substring(temp1[j],1,1) == "T" & substring(temp1[j],2,2) == "C"))) temp2[j] = "Y"
		if((length(grep("X",temp1[j])) < 1) & ((substring(temp1[j],1,1) == "G" & substring(temp1[j],2,2) == "T") | (substring(temp1[j],1,1) == "T" & substring(temp1[j],2,2) == "G"))) temp2[j] = "K"
		
		if(temp2[j] == "X") temp2[j] = "N"
		}
	if((i %% 500) == 0) print(paste(i, "of",nrow(all_comparisons_4),"the time is now:",Sys.time()))
	all_comparisons_4_nj_preformat[i,4:col] = temp2
	}
write.table(all_comparisons_4_nj_preformat,"results_6species/all_comparisons_4_nj_preformat", row.names = F, col.names = T, quote = T)	

##############
####NJ TREE###
##############	
library(ape)
all_comparisons_4_nj_preformat = as.matrix(read.delim("results_6species/all_comparisons_4_nj_preformat", header  = T, sep = " "))

all_comparisons_4_nj_format = t(tolower(all_comparisons_4_nj_preformat))[4:ncol(all_comparisons_4_nj_preformat),]
	
colnames(all_comparisons_4_nj_format) = apply(cbind(all_comparisons_4_nj_preformat[,1],as.numeric(all_comparisons_4_nj_preformat[,2])),1,paste,collapse = "_")
	
all_comparisons_4_nj_format_1 = as.DNAbin(all_comparisons_4_nj_format) # transform DNA matrix to DNAbin object
all_comparisons_4_nj_format_2 = dist.dna(all_comparisons_4_nj_format_1, pairwise.deletion = T) # distance matrix

#all_comparisons_4_nj_format_2[(all_comparisons_4_nj_format_2) == Inf] = 10
#all_comparisons_4_nj_format_2[is.na(all_comparisons_4_nj_format_2)] = 10

all_comparisons_4_nj_format_3 = nj(all_comparisons_4_nj_format_2) # neighbour joining tree

###################
###PLOTTING TREE###
###################

##############
####NEXUS FILE FOR SPLITSTREE###
##############
all_comparisons_4_nj_preformat = as.matrix(read.delim("results_6species/all_comparisons_4_nj_preformat", header  = T, sep = " "))
individuals = read.delim("reference/all_species_nov2012_cleaned.txt", header = T, stringsAsFactors = F) 
sequences = c(1:(ncol(all_comparisons_4_nj_preformat) - 3))
colnames(all_comparisons_4_nj_preformat) = gsub(".454Reads","",colnames(all_comparisons_4_nj_preformat),fixed = T)
colnames(all_comparisons_4_nj_preformat) = gsub(".white","",colnames(all_comparisons_4_nj_preformat),fixed = T)
colnames(all_comparisons_4_nj_preformat)[4:108] = paste(colnames(all_comparisons_4_nj_preformat)[4:108],individuals[,2], sep = "_") #species colnames


for(j in 4:ncol(all_comparisons_4_nj_preformat))
{
temp = colnames(all_comparisons_4_nj_preformat)[j]
temp = paste(temp,paste(all_comparisons_4_nj_preformat[,j], collapse = ""), sep = "     ") #create one consensus concatenated sequence per individual. 
sequences[j-3] = temp
}

nrow(all_comparisons_4_nj_preformat)
ncol(all_comparisons_4_nj_preformat)

all_helianthus_speciation_island = c("#NEXUS","BEGIN taxa;",paste("DIMENSIONS ntax=",ncol(all_comparisons_4_nj_preformat)-3,";", sep = ""),"TAXLABELS",colnames(all_comparisons_4_nj_preformat)[4:ncol(all_comparisons_4_nj_preformat)],";","END;","BEGIN characters;",paste("DIMENSIONS nchar=",nrow(all_comparisons_4_nj_preformat),";",sep = ""),"FORMAT","datatype=DNA","missing=N","gap=-","symbols=\"A C G T\"","labels","interleave",";","MATRIX",sequences, ";","END;") #CRAZY splitstree formating

write.table(all_helianthus_speciation_island,"results_6species/all_helianthus_speciation_island_12k.nex", sep = "", quote = F, row.names = F,col.names = F)

############
###SANDBOX### 
############

