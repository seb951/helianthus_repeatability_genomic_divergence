#!/usr/bin/Rscript --verbose
args = commandArgs(TRUE)
pair = as.numeric(args[1])
cpu =  as.numeric(args[2])# how many cpus can you afford (note that it may sometimes run (cpu + 2) CPUs, so allow for some buffering...

###This script does the the SNP calling using samtools mpileup and bcftools. It produces a large table of putative SNP which need to be filtered. 


#
###
###calling SNP###
###
#what to use: ANN DEB PET ARG BOL EXI
aa = Sys.time()
setwd("/home/seb/Documents/repeatability") #set up working directory 


#individuals of interest (454 individuals removed, weird individuals based on NatCom paper also removed#
#also remove the three ''weird" annuus bolanderii and the duplicated exilis (Ames) sample. this info was identified based on a splitstree with 20k SNP.
#individuals = individuals[-c(11,12,14,16),] #Ames7109_bol, BOL2436_bol,BOL1024_bol,BOL1023_bol
#if(pair == 2) individuals = individuals[regexpr("EAST|WEST",individuals[,2])>0,]
individuals = read.delim("reference/all_species_mar2013.txt", header = T, stringsAsFactors = F) 
individuals = cbind(individuals,0)
for(i in 1:nrow(individuals)) {individuals[i,5] = strsplit(individuals[i,1], split = "/")[[1]][length(strsplit(individuals[i,1], split = "/")[[1]])]}

#individuals of interest (454 individuals removed, weird individuals based on NatCom paper also removed#
#also remove the three ''weird" annuus bolanderii and the duplicated exilis (Ames) sample. this info was identified based on a splitstree with 20k SNP.

#individuals = read.delim("reference/all_species_nov2012_cleaned.txt", header = T, stringsAsFactors = F) #individuals of interest (454 individuals removed, weird individuals based $
#individuals[,1] = gsub("/home/transfer/Documents/new/bam_files_jan2012/","/SciBorg/array0/renaut/speciation_islands_individuals/alignments/",individuals[,1])
#individuals[1:16,1] = gsub("/SciBorg/array0/renaut/speciation_islands_individuals/alignments/","/SciBorg/array0/renaut/repeatability/alignments/",individuals[1:16,1])
#individuals[100:105,1] = gsub("/home/seb/Documents/repeatability/alignments/","/SciBorg/array0/renaut/repeatability/alignments/",individuals[100:105,1])

#individuals = individuals[-c(11,12,14,16),] #Ames7109_bol, BOL2436_bol,BOL1024_bol,BOL1023_bol
if(pair == 1) individuals = individuals[regexpr("deb|pet",individuals[,2])>0,]
if(pair == 2) individuals = individuals[regexpr("EAST|WEST",individuals[,2])>0,]
if(pair == 3) individuals = individuals[regexpr("ann|arg",individuals[,2])>0,]

#for(i in 1:nrow(individuals)) {individuals[i,1] = gsub("^","alignments/",strsplit(individuals[i,1], split = "/")[[1]][7])}

all_final = paste("alignements/",individuals[,5],".sorted.bam",sep = "",collapse = " ")

###mpileup###
split = read.table("reference/split", stringsAsFactors = F) #run mpileup gene per gene, this way you can parallelize.#

mpileup_dir = paste("mpileup_",pair,sep = "")
mpileup_cmd_dir =  paste("mpileup_",pair,"/cmd",sep = "")
system(paste("mkdir",mpileup_dir))
system(paste("mkdir",mpileup_cmd_dir))

setwd(paste("~/Documents/repeatability/",mpileup_dir,sep = "")) #set up new working directory
user = "seb"

for(j in 1:nrow(split))
#for(j in 1:200)#2000 genes only
{
system("ps -u renaut | grep 'samtools' | wc -l >prog")
mpileup = paste("samtools mpileup -C50 -I -ugf ../reference/HA412_trinity_noAltSplice_400bpmin.fa ", all_final, " -r ",split[j,1]," | bcftools view -cvg - > variants",j,".raw.vcf", sep = "")
Sys.sleep(ifelse(read.table("prog") <= cpu,0, (read.table("prog") - cpu)^4    ))
mpileup_exe = paste("cmd/mpileup_",j,sep = "")
write.table(mpileup,mpileup_exe, row.names = F, col.names = F, quote = F)
system(paste("chmod +x",mpileup_exe))
system(paste("nohup ./", mpileup_exe, ">log&",sep = ""))
}
print(paste("It took",Sys.time() - aa, "to run 5000 SNP with mpileup", "FUCK YEAH"))
#cat the variants_j.raw.vcf files. 
Sys.sleep(300) #tidy up everything
system("echo -n '' > 0_51k_variants.raw.vcf") #final file
system("cat variants*  >>0_51k_variants.raw.vcf") # cat verything
system("grep '#' -v  0_51k_variants.raw.vcf >0_51k_variants.clean.vcf") #tidy up
system("rm variants*")
system("rm cmd/*")


###Pre-cleaning step###
	system("wc -l 0_51k_variants.clean.vcf >wordcount1")
	wordcount = read.table("wordcount1")
	system("awk 'NR == 27' 0_51k_variants.raw.vcf  >header")
	header = read.delim("header", header = F, stringsAsFactors = F)
	header = gsub("alignments/","",header[10:length(header)], fixed = T)
	header = gsub(".fq.sorted.bam","",header, fixed = T)
	header = gsub(".sorted.bam","",header, fixed = T)
	header = gsub("_R1.fq.bz2","",header, fixed = T)
	header = gsub("/SciBorg/array0/renaut/repeatability/","",header, fixed = T)
	header = gsub("/SciBorg/array0/renaut/speciation_islands_individuals/","", header, fixed = T)
	write.table(t(as.matrix(c("reference", "position","ref_allele",header))),"snp_table",row.names = F, col.names = F, quote = F, sep = " ")

temp = NULL
	mis = NULL
	con = file("0_51k_variants.clean.vcf")
	open(con)
#	for(i in 1:10100) 
	for(i in 1:as.numeric(wordcount[1]))
	{
		x = readLines(con,1) #you will now read the mega big 0-51k_variants.clean.vcf file line by line in R. 
		xx = strsplit(x, split = "\t")[[1]]
		xxx = xx[(length(xx)-(length(header) - 1)):length(xx)]
		
		for(q in 1:length(xxx))	# this is to replace low quality calls (maximum Phred-scaled genotype likelihoods below 20). ie. essentially, you need at least 2 reads to call a SNP!
		{
		temp = strsplit(xxx[q], split = ":|,")[[1]]
		if(max(as.numeric(temp[2:4])) < 20) temp[1] = "XX" else if(min(as.numeric(temp[2:4])) == temp[2]) temp[1] = "RR" else if(min(as.numeric(temp[2:4])) == temp[3]) temp[1] = "AR" else if(min(as.numeric(temp[2:4])) == temp[4]) temp[1] = "AA"
		xxx[q] = temp[1]
		}
		xxx = gsub("R",xx[4],xxx)
		xxx = gsub("A",xx[5],xxx)
		xxx = c(xx[c(1,2,4)],xxx)
		if(nchar(xx[5]) > 1) xxx[4:length(xxx)] = rep("XX",length(xxx)-3)
		
		if((length(xxx[xxx == "XX"]) / length(xxx[4:length(xxx)])) < 0.5)	cat(t(as.matrix(xxx)),file = "snp_table",append = T, fill = F, "\n") else mis = c(mis, i) #only cat loci with less than 40% missing data. 
	#	if(regexpr(xxx[3],paste(xxx[4:21],collapse = "")) < 0) temp = rbind(temp,xxx)
		if(i %% 1000 == 0) print(paste(i,"of", wordcount[1], Sys.time()))
		}
close(con)
system(paste("cat snp_table | sed 's/[ \t]*$//' >snp_table_",pair,sep = ""))
system("rm snp_table") #### rm snp_table








