#!/usr/bin/Rscript --verbose
args = commandArgs(TRUE)
from = as.numeric(args[1])
to = as.numeric(args[2])

###This scrip does the actual aligments using bwa and indexing the files with samtools 


###
###Step 1 Prepare things up 
###
setwd("/home/seb/Documents/repeatability") #set up working directory 

individuals = read.delim("reference/all_species_mar2013.txt", header = T, stringsAsFactors = F) #individuals of interest
individuals = cbind(individuals,0)
for(i in 1:nrow(individuals)) {individuals[i,5] = strsplit(individuals[i,1], split = "/")[[1]][length(strsplit(individuals[i,1], split = "/")[[1]])]}

#system("bwa index reference/HA412_trinity_noAltSplice_400bpmin.fa") #index bwa reference
#system("samtools faidx reference/HA412_trinity_noAltSplice_400bpmin.fa") #index samtools reference

###
###Step 2. BWA SAMTOOLS.
###

#for(i in 1: nrow(individuals))
for(i in from:to)
{
	if(length(grep("bz2",individuals[i,5])) == 1) #1. unzip bzip2 files
		{
			#1. unzip 
			unzip1= paste("bzip2 -dk", individuals[i,1])
			unzip2= paste("bzip2 -dk",sub("_R1.fq.bz2","_R2.fq.bz2",individuals[i,1], fixed = T))
			fastq_1 = sub("_R1.fq.bz2","_R1.fq",individuals[i,1], fixed = T);fastq_1_res = sub("_R1.fq.bz2","_R1.fq",individuals[i,5], fixed = T)
			fastq_2 = sub("_R1.fq.bz2","_R1.fq",individuals[i,1], fixed = T);fastq_2_res = sub("_R1.fq.bz2","_R2.fq",individuals[i,5], fixed = T)
			system(unzip1)
			system(unzip2)
		}
	
	if(length(grep("gz",individuals[i,5])) == 1) #1. unzip gz files
		{
			#1. unzip 
			unzip1 = paste("gzip -dc ",sub(".fastq.gz","_1.fastq.gz",individuals[i,1], fixed = T)," >", sub(".fastq.gz","_1.fastq",individuals[i,1]), sep = "")
			unzip2 = paste("gzip -dc ",sub(".fastq.gz","_2.fastq.gz",individuals[i,1], fixed = T)," >", sub(".fastq.gz","_2.fastq",individuals[i,1]), sep = "")
			fastq_1 = sub(".fastq.gz","_1.fastq",individuals[i,1], fixed = T);fastq_1_res = sub(".fastq.gz","_1.fastq",individuals[i,5], fixed = T)
			fastq_2 = sub(".fastq.gz","_2.fastq",individuals[i,1], fixed = T);fastq_2_res = sub(".fastq.gz","_2.fastq",individuals[i,5], fixed = T)
			system(unzip1)
			system(unzip2)
		}			
		
	if(length(grep("Linux",individuals[i,5])) == 1) #1. unzip gz files
		{
			#1. unzip 
			fastq_1 = sub(".fq","_1.fq",individuals[i,1], fixed = T);fastq_1_res = sub(".fq","_1.fq",individuals[i,5], fixed = T)
			fastq_2 = sub(".fq","_2.fq",individuals[i,1], fixed = T);fastq_2_res = sub(".fq","_2.fq",individuals[i,5], fixed = T)
			}			
		
#2.align
bwa_aln1 = bwa_aln2 = "ls"
	if((individuals[i,3] == "ill") & (individuals[i,4] == "sanger")) # new illumina quality format
{bwa_aln1 = paste("bwa aln -t 2 -q 20 reference/HA412_trinity_noAltSplice_400bpmin.fa ",fastq_1,"  >alignments/",fastq_1_res, ".sai", sep = "");
bwa_aln2 = paste("bwa aln -t 2 -q 20 reference/HA412_trinity_noAltSplice_400bpmin.fa ",fastq_2,"  >alignments/",fastq_2_res, ".sai", sep = "")}

	if((individuals[i,3] == "ill") & (individuals[i,4] == "Ill1.3")) # old illumina quality format
{bwa_aln1 = paste("bwa aln -t 2 -I -q 20 reference/HA412_trinity_noAltSplice_400bpmin.fa ",fastq_1,"  >alignments/",fastq_1_res, ".sai", sep = "");
bwa_aln2 = paste("bwa aln -t 2 -I -q 20 reference/HA412_trinity_noAltSplice_400bpmin.fa ",fastq_2,"  >alignments/",fastq_2_res, ".sai", sep = "")}

#3. Sampe
bwa_sampe = "ls"
bwa_sampe = paste("bwa sampe reference/HA412_trinity_noAltSplice_400bpmin.fa alignments/",fastq_1_res, ".sai"," alignments/",fastq_2_res, ".sai ",fastq_1," ",fastq_2," >alignments/",individuals[i,5], ".sam", sep = "")

#4. samview
sam_view = paste("samtools view -bS -o alignments/",individuals[i,5], ".bam", " alignments/",individuals[i,5], ".sam", sep = "") #sam to bam

#5. samsort
sam_sort = paste("samtools sort alignments/",individuals[i,5], ".bam ","alignments/",individuals[i,5], ".sorted",sep = "") #sort the bam file

#6. samindex
sam_index = paste("samtools index alignments/",individuals[i,5],".sorted.bam",sep = "") #index bam file

#run commands
system(bwa_aln1)
system(bwa_aln2) 
system(bwa_sampe)
system(sam_view)
system(sam_sort)
system(sam_index)
}





