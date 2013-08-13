
###This script generates various plots (genome scans along the chromosomes)

setwd("/home/seb/Documents/repeatability") #set up working directory #set up working directory 

################
### PLOTTING ###
################

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





par(mfrow = c(3,1))
c1 = read.table("results_123/deb_pet_fst_map_cluster",  header = T, stringsAsFactors = F);main = "DEB_PET"
c2 = read.table("results_123/bol_exi_fst_map_cluster",  header = T, stringsAsFactors = F);main = "BOL_EXI"
c3 = read.table("results_123/ann_arg_fst_map_cluster",  header = T, stringsAsFactors = F);main = "ANN_ARG"
 color = c("#f8756b","#7aae00","#00bdc2")
 
plot(0,1,  type = "l", xlim = c(0,max(c1[,1])+1), ylim = c(0,1.09), xaxt = "n",yaxt = "n", lwd = 1,col = "black", main = main, xlab = "", ylab = "")
 points(c1[c1[,7] < 0.001,1],rep(1,length(c1[c1[,7] < 0.001,1])),col = color[1], pch = 19, cex = 0.5)
 points(c2[c2[,7] < 0.001,1],rep(0.9,length(c2[c2[,7] < 0.001,1])), col = color[2], pch = 19, cex = 0.5)
 points(c3[c3[,7] < 0.001,1],rep(0.8,length(c3[c3[,7] < 0.001,1])), col = color[3], pch = 19, cex = 0.5)


for(pair in c(1,2,3))
{
################
### PLOTTING ###
################
#if(a == 7) c1 =  read.table("results_6species/ann_arg_fst_map_cluster",  header = T, stringsAsFactors = F)
if(pair == 1) {c1 = read.table("results_123/deb_pet_fst_map_cluster_nonoverlapping_windows",  header = T, stringsAsFactors = F);main = "DEB_PET";c1[c1[,3] == 0,3] = NA}
if(pair == 2) {c2 = read.table("results_123/bol_exi_fst_map_cluster_nonoverlapping_windows",  header = T, stringsAsFactors = F);main = "BOL_EXI";c2[c2[,3] == 0,3] = NA}
if(pair == 3) {c3 = read.table("results_123/ann_arg_fst_map_cluster_nonoverlapping_windows",  header = T, stringsAsFactors = F);main = "ANN_ARG";c3[c3[,3] == 0,3] = NA}

###fst outlier regions
for(i in 1:17) #17 linkage groups
{
	cc1= c1[c1[,2] == i, ]

	if(i != 17) cc2 = c1[c1[,2] == (i+1), ] else {cc2 = cc1[cc1[,1] == max(cc1[,1]),]; cc2 = rbind(cc2,cc2)}
	if(i == 1) plot(0,1,  type = "l", xlim = c(0,max(c1[,1])+1), ylim = c(0,1.09), xaxt = "n", lwd = 1,col = "black", main = main, xlab = "", ylab = "")
	if(i == 1) rect(0,0,min(cc2[cc2[,3]>0,1],na.rm = T),1, border = F, col = ifelse(i %% 2 == 0,colors()[606],colors()[334]))
	if(i != 1) rect(min(cc1[cc1[,3]>0,1],na.rm = T),0,min(cc2[cc2[,3]>0,1],na.rm = T),1, border = F, col = ifelse(i %% 2 == 0,colors()[606],colors()[334]))

	points(cc1[cc1[,3]>0,1],cc1[cc1[,3] > 0 ,3], type = "l", xlim = c(0,max(c1[,1])+1), ylim = c(0,1.09), lwd  = 0.3, col =  color[pair])
	
	#gray coloring ### colors()[153],colors()[230] 
	#colours:  ifelse(i %%2 == 0,"darkred","darkblue")
	if(i == 17) abline(v = cc1[132,1], lwd = 1, col = "darkgreen")
}
#points(c1[c1[,7] < 0.001,1],rep(1.09,length(c1[c1[,7] < 0.001,1])), col = color[pair], pch = 19)
axis(1, at = max_position[,3], label = c(1:17))
}

dev.print(device=svg, "~/Desktop/islands_mean_fst.svg", onefile=FALSE)
dev.off()
################
### PLOT only certain chromosomes ###
################
#if(a == 7) c1 =  read.table("results_6species/ann_arg_fst_map_cluster",  header = T, stringsAsFactors = F)
c1 = read.table("results_123/deb_pet_fst_map_cluster",  header = T, stringsAsFactors = F)#;c1[c1[,3] == 0,3] = NA
c2 = read.table("results_123/bol_exi_fst_map_cluster",  header = T, stringsAsFactors = F)#;c2[c2[,3] == 0,3] = NA
c3 = read.table("results_123/ann_arg_fst_map_cluster",  header = T, stringsAsFactors = F)#;c3[c3[,3] == 0,3] = NA

pos_per_chr = c(1:3047)

for(i in 1:3047)
{
temp = map[signif(map[,4],6) == signif(c1[i,1],6),3]
pos_per_chr[i] = as.numeric(temp[1])
}

main = c("DEB_PET","BOL_EXI", "ANN_ARG")

par(mfrow = c(3,1))
for(chr in c(1,3,7))
{


plot(pos_per_chr[c1[,2] == chr],c1[c1[,2] == chr,3],xaxt = "n",col = color[1], type = "l", ylim = c(-0.1,1), lwd = 3,xlab = "cM",ylab = "Fst")
points(pos_per_chr[c2[,2] == chr],c2[c2[,2] == chr,3],xaxt = "n",col = color[2], type = "l", lwd = 3)
points(pos_per_chr[c3[,2] == chr],c3[c3[,2] == chr,3],xaxt = "n",col = color[3], type = "l", lwd = 3)

axis(1, at = pos_per_chr[c1[,2] == chr],label = rep("",length(c1[c1[,2] == chr,1])))
axis(1, at = c(0,20,40,60,80),label = c(0,20,40,60,80))
legend(x = 40,y = 0.8,fill = color, legend = main)

}

dev.print(device=svg, "results_123/chr1_3_7.svg", onefile=FALSE)
dev.off()



