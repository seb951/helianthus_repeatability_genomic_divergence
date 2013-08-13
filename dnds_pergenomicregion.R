#!/usr/bin/Rscript --verbose

###This script looks at the correlation between dn, ds and fst per genes and genomic regions



### loading required files ###
###Step  set up working directory and indiviudals of interest###
setwd("/home/seb/Documents/repeatability") #set up working directory 

c1 =  read.table("results_123/deb_pet_fst_map_cluster_nonoverlapping_windows",  header = T, stringsAsFactors = F)
c2 =  read.table("results_123/bol_exi_fst_map_cluster_nonoverlapping_windows",  header = T, stringsAsFactors = F)
c3 =  read.table("results_123/ann_arg_fst_map_cluster_nonoverlapping_windows",  header = T, stringsAsFactors = F)
dnds_deb_pet = as.matrix(read.table("results_123/deb_pet_fst_map_dnds.out",stringsAsFactors = F,header = T))
dnds_bol_exi = as.matrix(read.table("results_123/bol_exi_fst_map_dnds.out",stringsAsFactors = F, header = T))
dnds_ann_arg = as.matrix(read.table("results_123/ann_arg_fst_map_dnds.out",stringsAsFactors = F, header = T))

dnds_deb_pet[as.numeric(dnds_deb_pet[,2]) == 0,2] = NA
dnds_deb_pet[as.numeric(dnds_deb_pet[,3]) == 0,3] = NA
dnds_bol_exi [as.numeric(dnds_bol_exi [,4]) == 0,4] = NA
dnds_bol_exi [as.numeric(dnds_bol_exi [,2]) == 0,2] = NA
dnds_bol_exi [as.numeric(dnds_bol_exi [,3]) == 0,3] = NA
dnds_ann_arg[as.numeric(dnds_ann_arg[,4]) == 0,4] = NA
dnds_ann_arg[as.numeric(dnds_ann_arg[,2]) == 0,2] = NA
dnds_ann_arg[as.numeric(dnds_ann_arg[,3]) == 0,3] = NA
dnds_ann_arg[as.numeric(dnds_ann_arg[,4]) == 0,4] = NA

###	###
###	per gene dnds 
###	###
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

###	###	###	###
###	add map info	###
###	###	###	###
map = read.delim("reference/ordered_transcriptome_trinity.txt", header = F, stringsAsFactors = F) # new map
map = map[map[,3] != "no",];map[,2] = as.numeric(map[,2]);map[,3] = as.numeric(map[,3])
map = cbind(map,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA)
colnames(map) = c("name","LG","map_centiMorgan","unique_position","nbnuc","dn_deb_pet","ds_deb_pet","dnds_deb_pet","dn_bol_exi","ds_bol_exi","dnds_bol_exi","dn_ann_arg","ds_ann_arg","dnds_ann_arg")

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
	if(i %% 1000 == 0) print(paste(i,Sys.time()))
}
max_position[,3] = (max_position[,3] +1)
#add dnds
for(i in 1:nrow(map))
#for(i in 1:200)
{
x = dnds[,1] %in% map[i,1]

if(nrow(dnds[x == T,]) == 1)	map[i,6:14] = dnds[x == T,2:10]

if(i %% 1000 == 0) print(i)
}

write.table(dnds,"results_123/dnds_all3",col.names = T, row.names = F, quote = F)

###dnds per genomic position

#### CLUSTERING ###
###################
gmean=function(y=y,x=x, q=q,d=d){  # y=values, x=snp pos vector, q=query pos, d=dist:    Rose's gaussian average
	weights=rep(NA,length(x))
	weights[abs(x-q)<3*d]=exp(-(x[abs(x-q)<3*d]-q)^2/(2*d^2))
	weighted=weights*y/sum(weights,na.rm=T)
	return(sum(weighted,na.rm=T))
}

####running window - distance based###
 	d = 0.001 #distance cutoff in cM
	
	c = cbind(unique(as.numeric(map[,4])),0,0,0,0,0,0,0,0,0,0) #position you are interrogating along with the result of the distance based running mean. 
	
	c = matrix(as.numeric(c), ncol = 11)
	colnames(c) =  c("unique_position","LG","dn_deb_pet","ds_deb_pet","dnds_deb_pet","dn_bol_exi","ds_bol_exi","dnds_bol_exi","dn_ann_arg","ds_ann_arg","dnds_ann_arg")

	for(i in 1:nrow(c))
	{
	c_temp = map[as.numeric(map[,4]) == c[i,1],]
	c[i,2] = c_temp[,2][1]
	}

for(pair in 1:3)
{
	for(chr in 1:17)
		{
		chromo_temp = map[(as.numeric(map[,2]) == chr),] # subset of chromosome a
	
		for(i in c(1:nrow(c))[c[,2] == chr]) 
			{
			temp = chromo_temp[chromo_temp[,4] > (c[i,1]-d) & chromo_temp[,4] < (c[i,1]+d),] # mean fst at unique position in a sliding window of X cM
			
			#if(nrow(temp) >0) c[i,2] = gmean(temp[,((pair*3)+3)], temp[,4],c[i,1],d)	 #gaussian mean dn
			#if(nrow(temp) >0) c[i,3] = gmean(temp[,((pair*3)+4)], temp[,4],c[i,1],d)	 #gaussian mean ds 
			#if(nrow(temp) >0) c[i,4] = gmean(temp[,((pair*3)+5)], temp[,4],c[i,1],d)	 #gaussian mean dnds
			
			if(nrow(temp) >0) c[i,(pair*3)] = mean(temp[,((pair*3)+3)],na.rm = T)	 #regular mean dn
			if(nrow(temp) >0) c[i,((pair*3)+1)] = mean(temp[,((pair*3)+4)],na.rm = T)	 #regular mean ds 
			if(nrow(temp) >0) c[i,((pair*3)+2)] = mean(temp[,((pair*3)+5)],na.rm = T)	 #regulat mean dnds		
			}

		print(paste(chr, Sys.time()))
}
}
write.table(c,"results_123/dnds_all3_pergomicregion",col.names = T, row.names = F, quote = F)

bb = Sys.time()

### ### ### ### ### ### ### ###
###dn ds and dnds correlated with fst###
### ### ### ### ### ### ### ###
if(file.exists("results_123/dnds_all3_pergomicregion") ==T) c = read.delim("results_123/dnds_all3_pergomicregion",header = T, stringsAsFactors = T,sep = " ")

par(mfrow = c(2,2))
plot(log(c[,4]),log(c[,7]),main = "ds-fst",ylab = "ds",xlab = "ds");text(0.8,0.02, "cor = 0.24");cor.test(log(c[(log(c[,4]) != -Inf) & (log(c[,7]) != -Inf)  ,4]),log(c[(log(c[,4]) != -Inf) & (log(c[,7]) != -Inf) ,7])) #ds1 versus ds2
plot(log(c[,7]),log(c[,10]),main = "ds-fst",ylab = "ds",xlab = "ds");text(0.8,0.06, "cor = 0.36");cor.test(log(c[(log(c[,7]) != -Inf) & (log(c[,10]) != -Inf)  ,7]),log(c[(log(c[,7]) != -Inf) & (log(c[,10]) != -Inf) ,10])) #ds3 versus ds4
plot(log(c[,4]),log(c[,10]),main = "ds-fst",ylab = "ds",xlab = "ds");text(0.8,2, "cor = ns");cor.test(log(c[(log(c[,4]) != -Inf) & (log(c[,10]) != -Inf)  ,4]),log(c[(log(c[,4]) != -Inf) & (log(c[,10]) != -Inf) ,10])) #ds2 versus ds4
dev.print(device=svg, "results_123/ds_per_region.svg",onefile=FALSE)
dev.off()

###ds and fst for each species pair.
par(mfrow = c(2,2))
plot(log(c1[,3]),log(c[,4]),main = "ds-fst",ylab = "ds",xlab = "fst",xlim = c(-4,0));text(-3,-6, "cor = 0.27");cor.test(log(c1[(log(c1[,3]) != -Inf) & (log(c[,4]) != -Inf)  ,3]),log(c[(log(c1[,3]) != -Inf) & (log(c[,4]) != -Inf) ,4])) #ds1 versus ds2
plot(log(c2[,3]),log(c[,7]),main = "ds-fst",ylab = "ds",xlab = "fst", xlim = c(-4,0));text(-3,-6, "cor = 0.19");cor.test(log(c2[(log(c2[,3]) != -Inf) & (log(c[,7]) != -Inf)  ,3]),log(c[(log(c2[,3]) != -Inf) & (log(c[,7]) != -Inf) ,7])) #ds3 versus ds4
plot(log(c3[,3]),log(c[,10]),main = "ds-fst",ylab = "ds",xlab = "fst",xlim = c(-4,0));text(-3,-6, "cor = 0.32");cor.test(log(c3[(log(c3[,3]) != -Inf) & (log(c[,10]) != -Inf)  ,3]),log(c[(log(c3[,3]) != -Inf) & (log(c[,10]) != -Inf) ,10])) #ds2 versus ds4
dev.print(device=svg, "results_123/ds_versus_fst_per_comparison.svg",onefile=FALSE)
dev.off()




