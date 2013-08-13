#!/usr/bin/Rscript --verbose

args = commandArgs(TRUE)
aa = as.numeric(args[1])

#This script is to annotate SNP into the reference sequences and then prepare the sequences to get fed into PAML. PAML then calculates dn and ds. And back in R to get a final matrix of dn, ds and dn/ds ratio. 


### loading required files ###
###Step  set up working directory and indiviudals of interest###
setwd("/home/seb/Documents/repeatability") #set up working directory 
snp_table = c("results_123/deb_pet_fst_map","results_123/bol_exi_fst_map","results_123/ann_arg_fst_map")
#snp_table = c("mpileup_123/snp_table_1","mpileup_123/snp_table_2","mpileup_123/snp_table_3")

individuals = read.delim("reference/all_species_nov2012_cleaned.txt", header = T, stringsAsFactors = F) 
individuals = individuals[-c(11,12,14,16),] #Ames7109_bol, BOL2436_bol,BOL1024_bol,BOL1023_bol

if(aa == 1) individuals = individuals[regexpr("deb|pet",individuals[,2])>0,]
if(aa == 2) individuals =individuals[regexpr("EAST|WEST",individuals[,2])>0,]
if(aa == 3) individuals =individuals[regexpr("ann|arg",individuals[,2])>0,]



comp5 =  read.table(snp_table[aa],  header = T, stringsAsFactors = F)

#if(aa != 2) colnames(comp5)[4:(ncol(comp5)-10)] = paste(colnames(comp5)[4:(ncol(comp5)-10)],individuals[,2],sep = "_")
#colnames(comp5)[4:(ncol(comp5))] = paste(colnames(comp5)[4:(ncol(comp5))],individuals[,2],sep = "_")
#comp5 = comp5[1:10000,]
###
#create a table of the consensus call for annuus and debilis.
###
counter_ACGTX = matrix(0,nrow = nrow(comp5), ncol = 5)
counter_ACGTX1 = matrix(0,nrow = nrow(comp5), ncol = 5)
colnames(counter_ACGTX) = c("A","C","G","T","X")
colnames(counter_ACGTX1) = c("A","C","G","T","X")
consensus = cbind(comp5[,c(1,2)],0,0,0,0,0,0,0,0)
colnames(consensus) = c("contig","site","fst","one","two", "syn","nonsyn","shared","private","noncoding")

for(i in 1:nrow(comp5))
#for(i in 1:10000)
		{
		if(aa == 1) {one = c(1:ncol(comp5))[regexpr("_deb",colnames(comp5)) >0]; two = c(1:ncol(comp5))[regexpr("_pet",colnames(comp5)) >0]}
		if(aa == 2) {one = c(1:ncol(comp5))[regexpr("_EAST",colnames(comp5)) >0]; two = c(1:ncol(comp5))[regexpr("_WEST",colnames(comp5)) >0]}
		if(aa == 3) {one = c(1:ncol(comp5))[regexpr("_ann",colnames(comp5)) >0]; two = c(1:ncol(comp5))[regexpr("_arg",colnames(comp5)) >0]}
	
		#majority base call for ONE. 	
		x = paste(comp5[i,one], collapse = "")
		xx = strsplit(x, split = "")
		counter_ACGTX[i,1] = length(c(1:nchar(x))[grepl("A", xx[[1]])])
		counter_ACGTX[i,2] = length(c(1:nchar(x))[grepl("C", xx[[1]])])
		counter_ACGTX[i,3] = length(c(1:nchar(x))[grepl("G", xx[[1]])])
		counter_ACGTX[i,4] = length(c(1:nchar(x))[grepl("T", xx[[1]])])
		counter_ACGTX[i,5] = length(c(1:nchar(x))[grepl("X", xx[[1]])])
		if(consensus[i,3] == 1) consensus[i,4] = names(sort(counter_ACGTX[i,1:4])[4])
		if(sort(counter_ACGTX[i,1:4])[3] == 0) consensus[i,4] =names(sort(counter_ACGTX[i,1:4])[4])
		if(sort(counter_ACGTX[i,1:4])[3] != 0) consensus[i,4] =names(sort(counter_ACGTX[i,1:4])[4])
		
		#majority base call for TWO. 
		x1 = paste(comp5[i,two], collapse = "")
		xx1 = strsplit(x1, split = "")
		counter_ACGTX1[i,1] = length(c(1:nchar(x1))[grepl("A", xx1[[1]])])
		counter_ACGTX1[i,2] = length(c(1:nchar(x1))[grepl("C", xx1[[1]])])
		counter_ACGTX1[i,3] = length(c(1:nchar(x1))[grepl("G", xx1[[1]])])
		counter_ACGTX1[i,4] = length(c(1:nchar(x1))[grepl("T", xx1[[1]])])
		counter_ACGTX1[i,5] = length(c(1:nchar(x1))[grepl("X", xx1[[1]])])
		if(consensus[i,3] == 1) consensus[i,5] =names(sort(counter_ACGTX1[i,1:4])[4])
		if(sort(counter_ACGTX1[i,1:4])[3] == 0) consensus[i,5] = names(sort(counter_ACGTX1[i,1:4])[4])
		if(sort(counter_ACGTX1[i,1:4])[3] != 0) consensus[i,5] = names(sort(counter_ACGTX1[i,1:4])[4])
		if((i %% 1000) == 0) print(paste(i,"of",nrow(comp5)," >The time is...",  Sys.time()))
		#shared or private snps?
		}

if(aa == 1) {consensus[,3] = comp5[,ncol(comp5)-8]; write.table(consensus, "results_123/consensus_final_deb_pet", row.names = F, col.names = T, quote = F, sep = "\t")}
if(aa == 2) {consensus[,3] = comp5[,ncol(comp5)-8]; write.table(consensus, "results_123/consensus_final_bol_exi", row.names = F, col.names = T, quote = F, sep = "\t")}
if(aa == 3) {consensus[,3] = comp5[,ncol(comp5)-8]; write.table(consensus, "results_123/consensus_final_ann_arg", row.names = F, col.names = T, quote = F, sep = "\t")}

consensus_table = c("results_123/consensus_final_deb_pet","results_123/consensus_final_bol_exi","results_123/consensus_final_ann_arg")

### Sébastien Renaut
### University of British Columbia, Vancouver 
### sebastien.renaut@gmail.com

#This scripts does five main things. Use is at your own risks!
#It should work fine on Linux or Mac. Maybe not on a PC.
#######1. It annotates all the reference sequences with all your SNP
#######2. It finds all the best OpenReadingFrames
#######3. It places the annotated reference sequences in the right orientation and reading frame. It removes alternate stop codons. 
#######4. It prepares the PAML files for codeml.
#######5. It parses the output files of codeml. 


###############################
###define some required parameters and file###
###############################
#wd = "/home/seb/Documents/helianthus_analysis/92samples_new_alignments_aug12/results" # working directory where files are. 
ref = "/home/seb/Documents/repeatability/reference/HA412_trinity_noAltSplice_400bpmin.fa" #raw reference transcriptome
getorf = "getorf" #where is the program "getorf" used to call OpenReadingFrames? If you already have ORF, you may simply supply an object (matrix) which contains name of contig, start position and end position as in the object "orf_positions"
#snp_table = "snp_table" #the SNP table. 
aa_code = "/home/seb/Documents/helianthus_analysis/92samples_new_alignments_aug12/reference/aa_code" #table with amino acid codes and codons correspondences.
		
#########################################
### UPLOAD AND ANNOTATE REFERENCE SEQUENCES#######
#########################################

###Replace names with unique numerical ID### 
reference_transcriptome = as.matrix(read.delim(ref, header = F,sep = "\t"))
reference_transcriptome_unique_ID = reference_transcriptome
reference_transcriptome_unique_ID[(regexpr(">", reference_transcriptome_unique_ID[,1],fixed = T) > 0),1] = paste(">",c(1:length(reference_transcriptome_unique_ID[(regexpr(">", reference_transcriptome_unique_ID[,1],fixed = T) > 0),1])), sep = "") 
write.table(reference_transcriptome_unique_ID,"unique_ID.fa", row.names = F, col.names = F, quote = F)

###create a objet reference_matrix which contains the annotated consensus sequences. 
reference_vector = c(1:nrow(reference_transcriptome))[(regexpr(">", reference_transcriptome[,1],fixed = T) > 0)] #where do sequences start#
reference_vector = cbind(reference_transcriptome[(regexpr(">", reference_transcriptome[,1],fixed = T) > 0),1], reference_vector,0)
reference_vector = rbind(reference_vector,c("null",nrow(reference_transcriptome),0)) #extra line

for(i in 1:(nrow(reference_vector)-1)) #add the sequence on one row
{reference_vector[i,3] = paste(reference_transcriptome[(as.numeric(reference_vector[i,2])+1): (as.numeric(reference_vector[(i+1),2])-1),1],collapse = "")} 

reference_matrix = cbind(gsub(">","",reference_vector[1:nrow(reference_vector),1]),reference_vector[,3],reference_vector[,3],reference_vector[,3]) #have two versions of the sequence for 2 pseudo haplotypes.
colnames(reference_matrix) = c("name","raw_sequence","ONE_haplo1","TWO_haplo1")

### annotate reference matrix with the consensus object (ie. a SNP table).
snp = read.delim(consensus_table[aa], header = T, stringsAsFactors = F)

for(i in 1:nrow(reference_matrix))
#for(i in 1:2000)
	{
	x = snp[,1] %in% reference_matrix[i,1]
	y = snp[x == T, ]

	if(nrow(y) == 1)  #only one SNP in contig
	{
	substring(reference_matrix[i,3], as.numeric(y[2]), as.numeric(y[2])) = substring(y[4],1,1);
	substring(reference_matrix[i,4], as.numeric(y[2]), as.numeric(y[2])) = substring(y[5],1,1)
	}

	if(nrow(y) > 1) for(j in 1:nrow(y)) #more than one SNP in contig
	{
	substring(reference_matrix[i,3], as.numeric(y[j,2]), as.numeric(y[j,2])) = substring(y[j,4],1,1)
	substring(reference_matrix[i,4], as.numeric(y[j,2]), as.numeric(y[j,2])) = substring(y[j,5],1,1)
	}

	if(i %% 500 == 0) print(paste(i, "of",nrow(reference_matrix),"sequences annotated, time is:",  Sys.time())) #progress report
}

#############################
### FIND BEST ORF using getorf#######
#############################
###Alternatively, if you already know were your ORF are, you can simply skip this section and supply a file which contains name of contig, start position and end position as in the object "orf_positions" created at the end of this section.

###RUN GETORF PACKAGE FROM THE EMBOSS PIPELINE (FOUND HERE: http://emboss.sourceforge.net/) ###
system(paste(getorf,"-sequence unique_ID.fa -minsize 300 -find 2 -outseq unique_ID_orf.fa"))

### PARSE THE OUTPUT OF GETORF###
out = as.matrix(read.delim("unique_ID_orf.fa", header = F, sep = " "))
out = out[,1:5]
out[,2] = gsub("^.","", out[,2]);out[,4] = gsub(".$","", out[,4]);out = gsub("SENSE)","", out, fixed = T);out = cbind(gsub("(REVERSE","", out, fixed = T), 0);out = rbind(out,">")

x = c(1:nrow(out))[(regexpr(">",out[,1],fixed = T) > 0)] #sequence starting lines
for(i in 1: (length(x)-1)) {out[x[i],5]  = paste(out[(x[i]+1): (x[(i+1)]-1),1], collapse = "")} #add full sequence to column 5

out_seq = out[(regexpr(">",out[,1],fixed = T) > 0),] #remove uneeded rows
out_seq[,3] =  abs(as.numeric(out_seq[,2]) - as.numeric(out_seq[,4]))  #sequence length
out_seq[,6] = apply(out_seq,2,substring,(regexpr(">", out_seq[,1],fixed = T)+1),(regexpr("_", out_seq[,1],fixed = T)-1))[,1] #unique ID 
unique_gene = unique(out_seq[,6])
unique_orf = NULL #These will be the longest ORF.

for(i in 1:(length(unique_gene)-1)) # this loop is to keep only the longest ORF in case there are several for one contig.
	{
	x = out_seq[,6] %in% unique_gene[i]
	temp = out_seq[x == T, ]
	if(length(temp) == 6) longest_temp = temp[as.numeric(temp[3]) == max(as.numeric(temp[3]))] else longest_temp = temp[as.numeric(temp[,3]) == max(as.numeric(temp[,3])),]
	if(length(longest_temp) == 6) unique_orf = rbind(unique_orf, longest_temp) else unique_orf = rbind(unique_orf,longest_temp[1,])
	}

### get the proper names back from the original file
	names = reference_transcriptome[(regexpr(">", reference_transcriptome[,1],fixed = T) > 0),1] 
	for(i in 1:nrow(unique_orf)) {unique_orf[i,1] = names[as.numeric(unique_orf[i,6])]}
	unique_orf[,6] = gsub(">","",unique_orf[,1], fixed = T)
	colnames(unique_orf) = c("contig","start","length","stop","sequence","contigs")

#These are the position of the ORF in the reference contigs.
orf_positions = unique_orf[,c(6,2,4)]


##########################################
### create a new matrix which contains the ANNOTATED ORF ###
##########################################
symbols =      c("A","C","G","T","M","R","W","S","Y","K","V","H","D","B")
replacements = c("t","g","c","a","k","y","w","s","r","m","b","d","h","v")

orf_consensus = cbind(reference_matrix[,1],0,0)#This object will contain the SNP annotated ORF
colnames(orf_consensus) = c("reference","ann","deb")

for(i in 1:nrow(reference_matrix))
{

x = unique_orf[,6] %in% reference_matrix[i,1]
x_positions = unique_orf[x == T, c(1,2,4)]

if((length(x_positions) == 3) & (as.numeric(x_positions[2]) < as.numeric(x_positions[3]))) orf_consensus[i,2:3] =  substring(reference_matrix[i,c(3,4)],as.numeric(x_positions[2]),as.numeric(x_positions[3]))
if((length(x_positions) == 3) & (as.numeric(x_positions[2]) > as.numeric(x_positions[3]))) # the orf for the reverse complement case.
	{
		orf_consensus[i,2:3] = substring(reference_matrix[i,c(3,4)],as.numeric(x_positions[3]),as.numeric(x_positions[2])) #substring the ORF
		orf_consensus[i,2] = paste(rev(strsplit(orf_consensus[i,2],"")[[1]]),collapse = "")  # ORF in the reverse order
		orf_consensus[i,3] = paste(rev(strsplit(orf_consensus[i,3],"")[[1]]),collapse = "")  # ORF in the reverse order
		for(s in 1:length(symbols)) orf_consensus[i,2:3] = gsub(symbols[s], replacements[s], orf_consensus[i,2:3]) #Complement sequence. 
	}
}
orf_consensus[,2:3] = toupper(orf_consensus[,2:3]) #gsub for small caps to capital letters. 

###This is to check for stop codons and cut the sequence until the STOP codon appears### ####STOP CODONS: TAA TRA TAG TGA#### 
premature_stop = NULL #list of premature STOP. If there are STOP, sequence is cut until the stop codon.

for(i in 1: nrow(orf_consensus))
#for(i in 1:2000)
{
	for(t in 1: (nchar(orf_consensus[i,3])/3))
		{
			for(k in 2:3)
			{
		if((substring(orf_consensus[i,k],((t*3)-2),(t*3)) == "TAA") | (substring(orf_consensus[i,k],((t*3)-2),(t*3)) == "TRA")  | (substring(orf_consensus[i,k],((t*3)-2),(t*3)) == "TAG") | (substring(orf_consensus[i,k],((t*3)-2),(t*3)) == "TGA")) {  premature_stop = c(premature_stop, orf_consensus[i,1]); orf_consensus[i,2:3] = substring(orf_consensus[i,2:3],1,((t-1)*3)) }
  			}
  		}	
  if(i %% 1000 == 0) print(paste(i, "of",nrow(orf_consensus),"ORF checked, time is:",Sys.time())) #progress report
}

########################################
###FORMAT PAML INPUT FILE, RUN PAML, PARSE OUTPUT###
########################################
temp = c(rbind(colnames(orf_consensus)[2:3], orf_consensus[i,2:3]))

orf_length = nchar(orf_consensus[,3]) #number of character in each sequence
orf_paml = NULL # this will be the object which contains all the sequences
name_paml = NULL # this will contain all the names 

for(i in 1: nrow(orf_consensus))
	{
	if(orf_length[i] > 300) temp = c(paste(2,orf_length[i]),c(rbind(colnames(orf_consensus)[2:3], orf_consensus[i,2:3])))
	if(orf_length[i] > 300) orf_paml = c(orf_paml, temp)
	if(orf_length[i] > 300) name_paml = c(name_paml, orf_consensus[i,1])
	}

write.table(c(orf_paml,"//end", name_paml),paste(snp_table[aa], "_file_for_CODEML",sep = ""), row.names = F, col.names = F, quote = F)

###modify the control file for codeml
codeml.ctl = read.table("codeml.ctl", stringsAsFactors = F, header = F, sep = "\t")
codeml.ctl[1,1] = paste("ndata =",length(orf_paml[orf_paml == "ann"])) #how many sequences
codeml.ctl[2,1] = paste("seqfile = ",snp_table[aa], "_file_for_CODEML",sep = "")  #sequence file
codeml.ctl[3,1] = "treefile = ann_deb.tree" #tree file
codeml.ctl[4,1] = paste("outfile = ",snp_table[aa], "_out",sep = "")

write.table(codeml.ctl,paste(snp_table[aa], "_codeml.ctl",sep = ""), row.names = F, col.names = F, quote = F)

###RUN PAML########
system(paste("codeml ",snp_table[aa], "_codeml.ctl",sep = ""))

###Parse the output of PAML###
ds = as.matrix(read.delim("2NG.dS", header = F))
dn = as.matrix(read.delim("2NG.dN", header = F))
dnds = matrix(0, nrow = (length(ds)/3), ncol = 4)
colnames(dnds) = c("name","dn","ds","dnds")

for(i in 1: (length(ds)/3))
	{
	dnds[i,2] = tail(strsplit(dn[(i*3),1],split = " ")[[1]],1)
	dnds[i,3] =  tail(strsplit(ds[(i*3),1],split = " ")[[1]],1)
	}

	dnds[,4] = signif(as.numeric(dnds[,2]) / as.numeric(dnds[,3]),4)
	dnds[,1] = name_paml
	dnds[,4] = gsub("Inf","NaN", dnds[,4])

write.table(dnds,paste(snp_table[aa], "_dnds.out", sep = ""), row.names = F, col.names = T, quote = F, sep = "\t") ###dnds results###


###tidy up
system(paste("rm 2NG* 4fold.nuc rub rst1 rst lnf unique_ID_orf.fa unique_ID.fa ",snp_table[aa], "_file_for_CODEML ",snp_table[aa], "_codeml.ctl ",snp_table[aa], "_out",sep = ""))

print(paste("We're done boyz. Let's go get a beer!", Sys.time()))



