#!/usr/bin/Rscript --verbose
args = commandArgs(TRUE)
from = as.numeric(args[1]) # Specify which sequences in "list_ind" file you want to count stats on. You should do all that you aligned otherwise, bug!
to = as.numeric(args[2])


###This script is to get a few summary statistics###


setwd("/home/seb/Documents/repeatability") #set up working directory 


##########################################
### count the number of aligned sequences (samtools idxstats) ###
##########################################

idxstats = function(list_ind= "reference/all_species_mar2013.txt",total_idxstats = NULL, from = 1, to = 1) {

if(from == 1) individuals = read.delim(list_ind, header = T, stringsAsFactors = F)# reference file
if(from == 1) for(i in 1:nrow(individuals)) {individuals[i,5] = strsplit(individuals[i,1], split = "/")[[1]][length(strsplit(individuals[i,1], split = "/")[[1]])]}
if(from == 1) individuals = cbind(individuals,0); individuals[,5] = gsub(".fastq.gz","",individuals[,5])

	for(i in from:to) {

	#step 1: generate the alignment statistics#
	command_idxstats = paste("samtools idxstats alignments/",individuals[i,5],".sorted.bam"," >idxstats_results",  sep = "")
	system(command_idxstats)
	
	#step 2: read the covage and the alignment stats and then update the result_matrix file
	idxstats_results = read.delim("idxstats_results", stringsAsFactors = F, header = F) #don't read last line.
	individuals[i,6] = sum(idxstats_results[,3]) / 1000000
	print(i)
	}
	
	system("rm idxstats_results")
	write.table(individuals,"results_6species/number_total_aligned",row.names = F)
}

############
### running ###
############
idxstats(from = 1, to = 105) ### count the number of aligned sequences


