# replicate2 - methylation data, follow ash criteria
#      age.stable/methylated.genes - unt1-unt10 <.05, csc10-unt10 > .2 #these genes had methylation changes because of cig smoke and not in untreated because of aging
#      age.intermediate - unt10-unt1 < .2, csc10-unt10 > .2 #these genes had methylation changes because of cig smoke and not in untreated because of aging
# 
# age-related methylation changes for chip seq
# - took methylated genes
# - plot genes methylated by aging alone
# - pick genes that showed methylated changes only in untreated cells (unt10-unt1) > .3
#      + treatment did not cause methylation (unt10-unt1) > .3, (csc10-csc1) < .1
# - genes methylated only by aging and beta values not changed at all due to treatment
# 
# pick up are probes, so gene X, some probes within this criteria and probes that dont fit criteria 
# genes taht have increased beta values with age, but slight increase in beta values from smoke treatment
# - that criteria, y-axis
# - pick pure age methylation genes to use in age.stable methylation analysis (age.specfic, treatment.specific)
# 
# age.related(ash criteria)- to see how consistent replicate experiments were # not used for chipseq analysis
#      - methylation.list: unt1-unt10 <= .05, csc10-unt10 >= .2
# age only (chipseq) - unt 10-unt1, csc10-csc1


##########  
# Purpose:
##########
# To generate genelists from raw idat files
# This particular script is specific to Michelle's untreated vs CSC-treated lung cells at 4 timepoints (1M, 6M, 10M, 15M)
print("MOST CURRENT")
#Load libraries
# source("http://www.bioconductor.org/biocLite.R")
# biocLite("IlluminaHumanMethylation450kmanifest")
library(limma)
library(gplots)
library(heatmap.plus)
library(IlluminaHumanMethylation450k.db)
library(minfi)
library(IlluminaHumanMethylation450kmanifest)
system.dir = "/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/"
# system.dir = "/Volumes/onc-analysis$/users/shwang26/"
#Load images created by Ash
load(paste0( system.dir, "/Michelle/ash/jhu-methylation-scripts/R_image/all.probes.hg19.RData"))
load(paste0( system.dir, "/Michelle/ash/jhu-methylation-scripts/R_image/cgi.xy.filter.RData"))
load(paste0( system.dir, "/Michelle/ash/jhu-methylation-scripts/R_image/probesWithin5000TSS.RData"))
load(paste0( system.dir, "/Michelle/ash/jhu-methylation-scripts/R_image/probe.to.gene.RData"))
load(paste0( system.dir, "/Michelle/ash/jhu-methylation-scripts/R_image/annotation.RData"))
#Load object, new beta tables generated with new data
load(file=paste0(system.dir, "/Michelle/Robjects/michelle_new-annotation_old.RData")) #REP1
# load(file=paste0(system.dir, "/Michelle/Robjects/michelle_new-annotation.RData")) #REP2
#load functions
source(paste0( system.dir, "/Michelle/Rscripts/methylation450k_functions.R"))
date <- gsub(" \\d+\\:\\d+\\:\\d+.*", "", Sys.time())
# Create new working directory
newdir <- paste0(system.dir, "Michelle/MethylationData/methylatedGeneLists/", date, "/rep1_agerelated/")
newdir
dir.create(newdir, recursive=TRUE)
setwd(newdir)

# #comment
# #Get paths to input files and folders
# idat <- paste0( system.dir, "/Michelle/MethylationData/idat_old")
# sheet_genelist <- paste0( system.dir, "/Michelle/ash/sheet_genelist.csv")
# sheet_deltaplot <- paste0( system.dir, "/Michelle/ash/sheet_deltaplot.csv")
# folder <- paste0( system.dir, "/Michelle/ash/new_iteration_replicate1/")
# dir.create(folder)
# #Fuction for reading and normalizing data
# beta.tab <- data_methylation(idat, folder)
# 
# 
# #Retain the probes which are +/-1500 from TSS (loaded from annotation.Rdata image (for file annotation.R))
# #Subset the beta.tab table to retain CpG-island probes.
# dim(beta.tab)
# head(beta.tab)
# beta.tab_CpGI <- beta.tab[rownames(beta.tab) %in% cgi.xy.filter, ]
# beta.tab_CpGI <- na.omit(beta.tab_CpGI)
# dim(beta.tab_CpGI)
# 
# ########################################################
# #Subset the beta.tab table to retain gene body probes.
# body <- beta.tab[names(body.probes), ]
# body <- na.omit(body)
# dim(body)
# 
# #Subset the beta.tab table to retain shelf probes.
# shelf <- beta.tab[names(shelf.probes), ]
# shelf <- na.omit(shelf)
# dim(shelf)
# 
# #Subset the beta.tab table to retain shore probes.
# shore <- beta.tab[names(shore.probes), ]
# shore <- na.omit(shore)
# dim(shore)
# ##########################################################
# 
# # GET CPGI PROBES
# #get +/- 1500 TSS 
# head(probesWithin5000TSS)
# tss.1500.probes <- as.character(probesWithin5000TSS[abs(probesWithin5000TSS$Dist.to.TSS) <= 1500, "ProbeName"])
# beta.tab_CpGI.PlusMinus1500_TSS <- beta.tab_CpGI[rownames(beta.tab_CpGI) %in% tss.1500.probes,] #these are the CpGI probes within +/-1500 bp from TSS
# beta.tab_CpGI.PlusMinus1500_TSS <- na.omit(beta.tab_CpGI.PlusMinus1500_TSS)
# beta.tab_CpGI<-data.frame(beta.tab_CpGI)
# beta.tab_CpGI.PlusMinus1500_TSS<-data.frame(beta.tab_CpGI.PlusMinus1500_TSS)
# dim(beta.tab_CpGI.PlusMinus1500_TSS)
# 
# #Write TSS +/-1500 probes to file with gene names
# if (length(rownames(beta.tab_CpGI.PlusMinus1500_TSS)) > 0) {
#   beta.tab_CpGI.PlusMinus1500_TSS <- cbind.data.frame(beta.tab_CpGI.PlusMinus1500_TSS, "Gene"=rep(NA, times=nrow(beta.tab_CpGI.PlusMinus1500_TSS)))
#   for(j in 1:nrow(beta.tab_CpGI.PlusMinus1500_TSS)) {
#     probe.id <- rownames(beta.tab_CpGI.PlusMinus1500_TSS)[j]
#     gene <- probe.to.gene$gene[probe.to.gene$probe==probe.id]
#     if(length(gene) == 1) { beta.tab_CpGI.PlusMinus1500_TSS$Gene[j] <- as.character(gene) }
#   }
#   write.csv(beta.tab_CpGI.PlusMinus1500_TSS,paste(folder, "agerelated_tss_CpGI_probes_old.csv",sep="") , quote=F, row.names=T) #row.names set to T to get the probe ids.
# }
# 
# save(beta.tab_CpGI,
#      beta.tab_CpGI.PlusMinus1500_TSS,
#      file=paste0(system.dir, "/Michelle/Robjects/michelle_new-annotation_old_agerelated.RData"))
# #comment
# transition
subdir <- "combined"
timepoint="combined"
type <- "TSS"
head(beta.tab_CpGI.PlusMinus1500_TSS)
tss <- beta.tab_CpGI.PlusMinus1500_TSS
tss_reorder <- cbind.data.frame(tss$unt.1, tss$csc.1, tss$unt.6, tss$csc.6, tss$unt.10, tss$csc.10, tss$unt.15, tss$csc.15)
colnames(tss_reorder) <- c("c1", "csc1", "c6", "csc6", "c10", "csc10", "c15", "csc15")
tss <- tss_reorder
colnames(tss)
rownames(tss) <- rownames(beta.tab_CpGI.PlusMinus1500_TSS)
dim(tss)

tss <- tss[which(tss$c1 <= .2),]
dim(tss)
head(tss)
#############################################################################################################################
# Generate genelists
#############################################################################################################################

######################################################################################################################################################################################
# STABLE HYPERMETHYLATION
######################################################################################################################################################################################
#####
# 6M HYPERMETHYLATION
#####
getwd()
#age.stable, treatment.specific methylation
print("checkpoint1")

# (tss$c6-tss$c1)>=.4 & (tss$csc6-tss$csc1)<=.05
age.hypermethylation.6M <- rownames(tss[which((tss$c6-tss$c1)>=.4 & abs(tss$csc6-tss$csc1)<=.05),])
genelist.name1="age-related.05.hypermethylation_6M_rep1_sb"
age.hypermethylation.6M <- tss[which((tss$c6-tss$c1)>=.4 & abs(tss$csc6-tss$csc1)<=.05),]
if(dim(age.hypermethylation.6M)[1]>0){
age.hypermethylation.6M.genes <- get.genelist(age.hypermethylation.6M, genelist.nam=genelist.name1)
head(age.hypermethylation.6M)
}
genelist.name2="csc-related.05.hypermethylation_6M_rep1_sb"
csc.hypermethylation.6M <- tss[which(abs(tss$c6-tss$c1)<=.05 & (tss$csc6-tss$csc1)>=.4),]
if(dim(csc.hypermethylation.6M)[1]>0){
     csc.hypermethylation.6M.genes <- get.genelist(csc.hypermethylation.6M, genelist.nam=genelist.name2)
     head(csc.hypermethylation.6M)
}
# length(age.hypermethylation.6M.genes)
# length(csc.hypermethylation.6M.genes)

intersect_age.csc <- intersect(age.hypermethylation.6M.genes, csc.hypermethylation.6M.genes)
age.related.stable <- age.hypermethylation.6M.genes
age.related.stable <- setdiff(age.hypermethylation.6M.genes, intersect_age.csc)
write.table(age.related.stable, paste(newdir,"/TSS_", genelist.name1, "_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.
csc.related.stable <- csc.hypermethylation.6M.genes
csc.related.stable <- setdiff(csc.hypermethylation.6M.genes, intersect_age.csc)
write.table(csc.related.stable, paste(newdir,"/TSS_", genelist.name2, "_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.

#####
# 10M HYPERMETHYLATION
#####
getwd()
#age.stable, treatment.specific methylation
print("checkpoint1")

# (tss$c10-tss$c1)>=.4 & (tss$csc10-tss$csc1)<=.05
# age.hypermethylation.10M <- rownames(tss[which((tss$c10-tss$c1)>=.4 & abs(tss$csc10-tss$csc1)<=.05),])
genelist.name1="age-related.05.hypermethylation_10M_rep1_sb"
age.hypermethylation.10M <- tss[which((tss$c10-tss$c1)>=.4 & abs(tss$csc10-tss$csc1)<=.05),]
age.hypermethylation.10M.genes <- get.genelist(age.hypermethylation.10M, genelist.nam=genelist.name1)
head(age.hypermethylation.10M)
genelist.name2="csc-related.05.hypermethylation_10M_rep1_sb"
csc.hypermethylation.10M <- tss[which(abs(tss$c10-tss$c1)<=.05 & (tss$csc10-tss$csc1)>=.4),]
csc.hypermethylation.10M.genes <- get.genelist(csc.hypermethylation.10M, genelist.nam=genelist.name2)
head(csc.hypermethylation.10M)

length(age.hypermethylation.10M.genes)
length(csc.hypermethylation.10M.genes)
intersect_age.csc <- intersect(age.hypermethylation.10M.genes, csc.hypermethylation.10M.genes)
age.related.stable <- setdiff(age.hypermethylation.10M.genes, intersect_age.csc)
write.table(age.related.stable, paste(newdir,"/TSS_", genelist.name1, "_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.
csc.related.stable <- setdiff(csc.hypermethylation.10M.genes, intersect_age.csc)
write.table(csc.related.stable, paste(newdir,"/TSS_", genelist.name2, "_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.

#####
# 15M HYPERMETHYLATION
#####
getwd()
#age.stable, treatment.specific methylation
print("checkpoint2")

# (tss$c15-tss$c1)>=.4 & (tss$csc15-tss$csc1)<=.05
# age.hypermethylation.15M <- rownames(tss[which((tss$c15-tss$c1)>=.4 & abs(tss$csc15-tss$csc1)<=.05),])
genelist.name1="age-related.05.hypermethylation_15M_rep1_sb"
age.hypermethylation.15M <- tss[which((tss$c15-tss$c1)>=.4 & abs(tss$csc15-tss$csc1)<=.05),]
if(dim(age.hypermethylation.15M)[1]>0){
     age.hypermethylation.15M.genes <- get.genelist(age.hypermethylation.15M, genelist.nam=genelist.name1)
}
head(age.hypermethylation.15M)
genelist.name2="csc-related.05.hypermethylation_15M_rep1_sb"
csc.hypermethylation.15M <- tss[which(abs(tss$c15-tss$c1)<=.05 & (tss$csc15-tss$csc1)>=.4),]
if(dim(csc.hypermethylation.15M)[1]>0){
     csc.hypermethylation.15M.genes <- get.genelist(csc.hypermethylation.15M, genelist.nam=genelist.name2)
}
head(csc.hypermethylation.15M)

length(age.hypermethylation.15M.genes)
length(csc.hypermethylation.15M.genes)
intersect_age.csc <- intersect(age.hypermethylation.15M.genes, csc.hypermethylation.15M.genes)
age.related.stable <- setdiff(age.hypermethylation.15M.genes, intersect_age.csc)
write.table(age.related.stable, paste(newdir,"/TSS_", genelist.name1, "_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.
csc.related.stable <- setdiff(csc.hypermethylation.15M.genes, intersect_age.csc)
write.table(csc.related.stable, paste(newdir,"/TSS_", genelist.name2, "_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.

######################################################################################################################################################################################
# TRUE INTERMEDIATE HYPERMETHYLATION
######################################################################################################################################################################################

#####
# 6M HYPERMETHYLATION
#####
print("checkpoint3")
genelist.name1="age-related.2.true.hypermethylation_6M_rep1_sb"
age.hypermethylation.6M <- tss[which((tss$c6-tss$c1)>=.4 & abs(tss$csc6-tss$csc1)>=.05 & abs(tss$csc6-tss$csc1)<=.2),]
if(dim(age.hypermethylation.6M)[1]>0){
     age.hypermethylation.6M.genes <- get.genelist(age.hypermethylation.6M, genelist.nam=genelist.name1)
}
head(age.hypermethylation.6M)
genelist.name2="csc-related.2.true.hypermethylation_6M_rep1_sb"
csc.hypermethylation.6M <- tss[which(abs(tss$c6-tss$c1)>=.05 & abs(tss$c6-tss$c1)<=.2 & (tss$csc6-tss$csc1)>=.4),]
if(dim(csc.hypermethylation.6M)[1]>0){
     csc.hypermethylation.6M.genes <- get.genelist(csc.hypermethylation.6M, genelist.nam=genelist.name2)
}
head(csc.hypermethylation.6M)

length(age.hypermethylation.6M.genes)
length(csc.hypermethylation.6M.genes)
intersect_age.csc <- intersect(age.hypermethylation.6M.genes, csc.hypermethylation.6M.genes)
age.related.stable <- setdiff(age.hypermethylation.6M.genes, intersect_age.csc)
write.table(age.related.stable, paste(newdir,"/TSS_", genelist.name1, "_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.
csc.related.stable <- setdiff(csc.hypermethylation.6M.genes, intersect_age.csc)
write.table(csc.related.stable, paste(newdir,"/TSS_", genelist.name2, "_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.


#####
# 10M HYPERMETHYLATION
#####
print("checkpoint3")
genelist.name1="age-related.2.true.hypermethylation_10M_rep1_sb"
age.hypermethylation.10M <- tss[which((tss$c10-tss$c1)>=.4 & abs(tss$csc10-tss$csc1)>=.05 & abs(tss$csc10-tss$csc1)<=.2),]
if(dim(age.hypermethylation.10M)[1]>0){
     age.hypermethylation.10M.genes <- get.genelist(age.hypermethylation.10M, genelist.nam=genelist.name1)
}
head(age.hypermethylation.10M)
genelist.name2="csc-related.2.true.hypermethylation_10M_rep1_sb"
csc.hypermethylation.10M <- tss[which(abs(tss$c10-tss$c1)>=.05 & abs(tss$c10-tss$c1)<=.2 & (tss$csc10-tss$csc1)>=.4),]
if(dim(csc.hypermethylation.10M)[1]>0){
     csc.hypermethylation.10M.genes <- get.genelist(csc.hypermethylation.10M, genelist.nam=genelist.name2)
}
head(csc.hypermethylation.10M)

length(age.hypermethylation.10M.genes)
length(csc.hypermethylation.10M.genes)
intersect_age.csc <- intersect(age.hypermethylation.10M.genes, csc.hypermethylation.10M.genes)
age.related.stable <- setdiff(age.hypermethylation.10M.genes, intersect_age.csc)
write.table(age.related.stable, paste(newdir,"/TSS_", genelist.name1, "_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.
csc.related.stable <- setdiff(csc.hypermethylation.10M.genes, intersect_age.csc)
write.table(csc.related.stable, paste(newdir,"/TSS_", genelist.name2, "_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.

#####
# 15M HYPERMETHYLATION
#####
print("checkpoint4")
genelist.name1="age-related.2.true.hypermethylation_15M_rep1_sb"
age.hypermethylation.15M <- tss[which((tss$c15-tss$c1)>=.4 & abs(tss$csc15-tss$csc1)>=.05 & abs(tss$csc15-tss$csc1)<=.2),]
age.hypermethylation.15M.genes <- get.genelist(age.hypermethylation.15M, genelist.nam=genelist.name1)
head(age.hypermethylation.15M)
genelist.name2="csc-related.2.true.hypermethylation_15M_rep1_sb"
csc.hypermethylation.15M <- tss[which(abs(tss$c15-tss$c1)<=.05 & abs(tss$c15-tss$c1)<=.2 & (tss$csc15-tss$csc1)>=.4),]
csc.hypermethylation.15M.genes <- get.genelist(csc.hypermethylation.15M, genelist.nam=genelist.name2)
head(csc.hypermethylation.15M)

length(age.hypermethylation.15M.genes)
length(csc.hypermethylation.15M.genes)
intersect_age.csc <- intersect(age.hypermethylation.15M.genes, csc.hypermethylation.15M.genes)
age.related.stable <- setdiff(age.hypermethylation.15M.genes, intersect_age.csc)
write.table(age.related.stable, paste(newdir,"/TSS_", genelist.name1, "_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.
csc.related.stable <- setdiff(csc.hypermethylation.15M.genes, intersect_age.csc)
write.table(csc.related.stable, paste(newdir,"/TSS_", genelist.name2, "_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.


# ######################################################################################################################################################################################
# # INTERMEDIATE HYPERMETHYLATION
# ######################################################################################################################################################################################
# 
# 
# #####
# # 6M HYPERMETHYLATION
# #####
# print("checkpoint3")
# genelist.name1="age-related.2.hypermethylation_6M_rep1_sb"
# age.hypermethylation.6M <- tss[which((tss$c6-tss$c1)>=.4 & abs(tss$csc6-tss$csc1)<=.2),]
# if(dim(age.hypermethylation.6M)[1]>0){
#      age.hypermethylation.6M.genes <- get.genelist(age.hypermethylation.6M, genelist.nam=genelist.name1)
# }
# head(age.hypermethylation.6M)
# genelist.name2="csc-related.2.hypermethylation_6M_rep1_sb"
# csc.hypermethylation.6M <- tss[which(abs(tss$c6-tss$c1)<=.2 & (tss$csc6-tss$csc1)>=.4),]
# if(dim(csc.hypermethylation.6M)[1]>0){
#      csc.hypermethylation.6M.genes <- get.genelist(csc.hypermethylation.6M, genelist.nam=genelist.name2)
# }
# head(csc.hypermethylation.6M)
# 
# length(age.hypermethylation.6M.genes)
# length(csc.hypermethylation.6M.genes)
# intersect_age.csc <- intersect(age.hypermethylation.6M.genes, csc.hypermethylation.6M.genes)
# age.related.stable <- setdiff(age.hypermethylation.6M.genes, intersect_age.csc)
# write.table(age.related.stable, paste(newdir,"/TSS_", genelist.name1, "_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.
# csc.related.stable <- setdiff(csc.hypermethylation.6M.genes, intersect_age.csc)
# write.table(csc.related.stable, paste(newdir,"/TSS_", genelist.name2, "_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.
# 
# 
# #####
# # 10M HYPERMETHYLATION
# #####
# print("checkpoint3")
# genelist.name1="age-related.2.hypermethylation_10M_rep1_sb"
# age.hypermethylation.10M <- tss[which((tss$c10-tss$c1)>=.4 & abs(tss$csc10-tss$csc1)<=.2),]
# if(dim(age.hypermethylation.10M)[1]>0){
#      age.hypermethylation.10M.genes <- get.genelist(age.hypermethylation.10M, genelist.nam=genelist.name1)
# }
# head(age.hypermethylation.10M)
# genelist.name2="csc-related.2.hypermethylation_10M_rep1_sb"
# csc.hypermethylation.10M <- tss[which(abs(tss$c10-tss$c1)<=.2 & (tss$csc10-tss$csc1)>=.4),]
# if(dim(csc.hypermethylation.10M)[1]>0){
#      csc.hypermethylation.10M.genes <- get.genelist(csc.hypermethylation.10M, genelist.nam=genelist.name2)
# }
# head(csc.hypermethylation.10M)
# 
# length(age.hypermethylation.10M.genes)
# length(csc.hypermethylation.10M.genes)
# intersect_age.csc <- intersect(age.hypermethylation.10M.genes, csc.hypermethylation.10M.genes)
# age.related.stable <- setdiff(age.hypermethylation.10M.genes, intersect_age.csc)
# write.table(age.related.stable, paste(newdir,"/TSS_", genelist.name1, "_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.
# csc.related.stable <- setdiff(csc.hypermethylation.10M.genes, intersect_age.csc)
# write.table(csc.related.stable, paste(newdir,"/TSS_", genelist.name2, "_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.
# 
# #####
# # 15M HYPERMETHYLATION
# #####
# print("checkpoint4")
# genelist.name1="age-related.2.hypermethylation_15M_rep1_sb"
# age.hypermethylation.15M <- tss[which((tss$c15-tss$c1)>=.4 & abs(tss$csc15-tss$csc1)<=.2),]
# age.hypermethylation.15M.genes <- get.genelist(age.hypermethylation.15M, genelist.nam=genelist.name1)
# head(age.hypermethylation.15M)
# genelist.name2="csc-related.2.hypermethylation_15M_rep1_sb"
# csc.hypermethylation.15M <- tss[which(abs(tss$c15-tss$c1)<=.2 & (tss$csc15-tss$csc1)>=.4),]
# csc.hypermethylation.15M.genes <- get.genelist(csc.hypermethylation.15M, genelist.nam=genelist.name2)
# head(csc.hypermethylation.15M)
# 
# length(age.hypermethylation.15M.genes)
# length(csc.hypermethylation.15M.genes)
# intersect_age.csc <- intersect(age.hypermethylation.15M.genes, csc.hypermethylation.15M.genes)
# age.related.stable <- setdiff(age.hypermethylation.15M.genes, intersect_age.csc)
# write.table(age.related.stable, paste(newdir,"/TSS_", genelist.name1, "_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.
# csc.related.stable <- setdiff(csc.hypermethylation.15M.genes, intersect_age.csc)
# write.table(csc.related.stable, paste(newdir,"/TSS_", genelist.name2, "_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.



#compare probes between age-related
stable <- read.table("~/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/Michelle/MethylationData/methylatedGeneLists/2015-06-01/rep1_agerelated/TSS_age-related.05.hypermethylation_10M_rep1_sb_beta.values.txt")
intermediate <- read.table("~/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/Michelle/MethylationData/methylatedGeneLists/2015-06-01/rep1_agerelated/TSS_age-related.2.true.hypermethylation_10M_rep1_sb_beta.values.txt")
head(stable)
intersect(rownames(stable), rownames(intermediate))
sort(intersect(unlist(strsplit(as.matrix(stable$Gene),";")), unlist(strsplit(as.matrix(intermediate$Gene),";"))))

#compare probes between csc-related
stable <- read.table("~/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/Michelle/MethylationData/methylatedGeneLists/2015-06-01/rep1_agerelated/TSS_csc-related.05.hypermethylation_10M_rep1_sb_beta.values.txt")
intermediate <- read.table("~/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/Michelle/MethylationData/methylatedGeneLists/2015-06-01/rep1_agerelated/TSS_csc-related.2.true.hypermethylation_10M_rep1_sb_beta.values.txt")
head(stable)
intersect(rownames(stable), rownames(intermediate))
sort(intersect(unlist(strsplit(as.matrix(stable$Gene),";")), unlist(strsplit(as.matrix(intermediate$Gene),";"))))



#############################################################################################################################
# HYPOMETHYLATION
#############################################################################################################################

######################################################################################################################################################################################
# STABLE HYPOMETHYLATION
######################################################################################################################################################################################

#####
# 6M HYPOMETHYLATION
#####
getwd()
#age.stable, treatment.specific methylation
print("checkpoint1")

# (tss$c6-tss$c1)>=.4 & (tss$csc6-tss$csc1)<=.05
# age.hypomethylation.6M <- rownames(tss[which((tss$c6-tss$c1)>=.4 & abs(tss$csc6-tss$csc1)<=.05),])
genelist.name1="age-related.05.hypomethylation_6M_rep1_sb"
age.hypomethylation.6M <- tss[which((tss$c6-tss$c1)<=-.4 & abs(tss$csc6-tss$csc1)<=.05),]
if(dim(age.hypomethylation.6M)[1]>0){
     age.hypomethylation.6M.genes <- get.genelist(age.hypomethylation.6M, genelist.nam=genelist.name1)
     head(age.hypomethylation.6M)
}
genelist.name2="csc-related.05.hypomethylation_6M_rep1_sb"
csc.hypomethylation.6M <- tss[which(abs(tss$c6-tss$c1)<=.05 & (tss$csc6-tss$csc1)<=-.4),]
if(dim(csc.hypomethylation.6M)[1]>0){
     csc.hypomethylation.6M.genes <- get.genelist(csc.hypomethylation.6M, genelist.nam=genelist.name2)
     head(csc.hypomethylation.6M)
}

length(age.hypomethylation.6M.genes)
length(csc.hypomethylation.6M.genes)
intersect_age.csc <- intersect(age.hypomethylation.6M.genes, csc.hypomethylation.6M.genes)
age.related.stable <- setdiff(age.hypomethylation.6M.genes, intersect_age.csc)
write.table(age.related.stable, paste(newdir,"/TSS_", genelist.name1, "_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.
csc.related.stable <- setdiff(csc.hypomethylation.6M.genes, intersect_age.csc)
write.table(csc.related.stable, paste(newdir,"/TSS_", genelist.name2, "_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.

#####
# 10M HYPOMETHYLATION
#####
getwd()
#age.stable, treatment.specific methylation
print("checkpoint1")

# (tss$c10-tss$c1)>=.4 & (tss$csc10-tss$csc1)<=.05
# age.hypomethylation.10M <- rownames(tss[which((tss$c10-tss$c1)>=.4 & abs(tss$csc10-tss$csc1)<=.05),])
genelist.name1="age-related.05.hypomethylation_10M_rep1_sb"
age.hypomethylation.10M <- tss[which((tss$c10-tss$c1)<=-.4 & abs(tss$csc10-tss$csc1)<=.05),]
age.hypomethylation.10M.genes <- get.genelist(age.hypomethylation.10M, genelist.nam=genelist.name1)
head(age.hypomethylation.10M)
genelist.name2="csc-related.05.hypomethylation_10M_rep1_sb"
csc.hypomethylation.10M <- tss[which(abs(tss$c10-tss$c1)<=.05 & (tss$csc10-tss$csc1)<=-.4),]
csc.hypomethylation.10M.genes <- get.genelist(csc.hypomethylation.10M, genelist.nam=genelist.name2)
head(csc.hypomethylation.10M)

length(age.hypomethylation.10M.genes)
length(csc.hypomethylation.10M.genes)
intersect_age.csc <- intersect(age.hypomethylation.10M.genes, csc.hypomethylation.10M.genes)
age.related.stable <- setdiff(age.hypomethylation.10M.genes, intersect_age.csc)
write.table(age.related.stable, paste(newdir,"/TSS_", genelist.name1, "_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.
csc.related.stable <- setdiff(csc.hypomethylation.10M.genes, intersect_age.csc)
write.table(csc.related.stable, paste(newdir,"/TSS_", genelist.name2, "_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.

#####
# 15M HYPOMETHYLATION
#####
getwd()
#age.stable, treatment.specific methylation
print("checkpoint2")

# (tss$c15-tss$c1)>=.4 & (tss$csc15-tss$csc1)<=.05
# age.hypomethylation.15M <- rownames(tss[which((tss$c15-tss$c1)>=.4 & abs(tss$csc15-tss$csc1)<=.05),])
genelist.name1="age-related.05.hypomethylation_15M_rep1_sb"
age.hypomethylation.15M <- tss[which((tss$c15-tss$c1)<=-.4 & abs(tss$csc15-tss$csc1)<=.05),]
if(dim(age.hypomethylation.15M)[1]>0){
     age.hypomethylation.15M.genes <- get.genelist(age.hypomethylation.15M, genelist.nam=genelist.name1)
}
head(age.hypomethylation.15M)
genelist.name2="csc-related.05.hypomethylation_15M_rep1_sb"
csc.hypomethylation.15M <- tss[which(abs(tss$c15-tss$c1)<=.05 & (tss$csc15-tss$csc1)<=-.4),]
if(dim(csc.hypomethylation.15M)[1]>0){
     csc.hypomethylation.15M.genes <- get.genelist(csc.hypomethylation.15M, genelist.nam=genelist.name2)
}
head(csc.hypomethylation.15M)

length(age.hypomethylation.15M.genes)
length(csc.hypomethylation.15M.genes)
intersect_age.csc <- intersect(age.hypomethylation.15M.genes, csc.hypomethylation.15M.genes)
age.related.stable <- setdiff(age.hypomethylation.15M.genes, intersect_age.csc)
write.table(age.related.stable, paste(newdir,"/TSS_", genelist.name1, "_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.
csc.related.stable <- setdiff(csc.hypomethylation.15M.genes, intersect_age.csc)
write.table(csc.related.stable, paste(newdir,"/TSS_", genelist.name2, "_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.

######################################################################################################################################################################################
# TRUE INTERMEDIATE HYPOMETHYLATION
######################################################################################################################################################################################

# #####
# # 6M HYPOMETHYLATION
# #####
# print("checkpoint3")
# genelist.name1="age-related.2.true.hypomethylation_6M_rep1_sb"
# age.hypomethylation.6M <- tss[which((tss$c6-tss$c1)<=-.4 & abs(tss$csc6-tss$csc1)>=.05 & abs(tss$csc6-tss$csc1)<=.2),]
# if(dim(age.hypomethylation.6M)[1]>0){
#      age.hypomethylation.6M.genes <- get.genelist(age.hypomethylation.6M, genelist.nam=genelist.name1)
# }
# head(age.hypomethylation.6M)
# genelist.name2="csc-related.2.true.hypomethylation_6M_rep1_sb"
# csc.hypomethylation.6M <- tss[which(abs(tss$c6-tss$c1)>=.05 & abs(tss$c6-tss$c1)<=.2 & (tss$csc6-tss$csc1)<=-.4),]
# if(dim(csc.hypomethylation.6M)[1]>0){
#      csc.hypomethylation.6M.genes <- get.genelist(csc.hypomethylation.6M, genelist.nam=genelist.name2)
# }
# head(csc.hypomethylation.6M)
# 
# length(age.hypomethylation.6M.genes)
# length(csc.hypomethylation.6M.genes)
# intersect_age.csc <- intersect(age.hypomethylation.6M.genes, csc.hypomethylation.6M.genes)
# age.related.stable <- setdiff(age.hypomethylation.6M.genes, intersect_age.csc)
# write.table(age.related.stable, paste(newdir,"/TSS_", genelist.name1, "_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.
# csc.related.stable <- setdiff(csc.hypomethylation.6M.genes, intersect_age.csc)
# write.table(csc.related.stable, paste(newdir,"/TSS_", genelist.name2, "_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.
# 
# #####
# # 10M HYPOMETHYLATION
# #####
# print("checkpoint3")
# genelist.name1="age-related.2.true.hypomethylation_10M_rep1_sb"
# age.hypomethylation.10M <- tss[which((tss$c10-tss$c1)<=-.4 & abs(tss$csc10-tss$csc1)>=.05 & abs(tss$csc10-tss$csc1)<=.2),]
# if(dim(age.hypomethylation.10M)[1]>0){
#      age.hypomethylation.10M.genes <- get.genelist(age.hypomethylation.10M, genelist.nam=genelist.name1)
# }
# head(age.hypomethylation.10M)
# genelist.name2="csc-related.2.true.hypomethylation_10M_rep1_sb"
# csc.hypomethylation.10M <- tss[which(abs(tss$c10-tss$c1)>=.05 & abs(tss$c10-tss$c1)<=.2 & (tss$csc10-tss$csc1)<=-.4),]
# if(dim(csc.hypomethylation.10M)[1]>0){
#      csc.hypomethylation.10M.genes <- get.genelist(csc.hypomethylation.10M, genelist.nam=genelist.name2)
# }
# head(csc.hypomethylation.10M)
# 
# length(age.hypomethylation.10M.genes)
# length(csc.hypomethylation.10M.genes)
# intersect_age.csc <- intersect(age.hypomethylation.10M.genes, csc.hypomethylation.10M.genes)
# age.related.stable <- setdiff(age.hypomethylation.10M.genes, intersect_age.csc)
# write.table(age.related.stable, paste(newdir,"/TSS_", genelist.name1, "_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.
# csc.related.stable <- setdiff(csc.hypomethylation.10M.genes, intersect_age.csc)
# write.table(csc.related.stable, paste(newdir,"/TSS_", genelist.name2, "_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.

#####
# 15M HYPOMETHYLATION
#####
print("checkpoint4")
genelist.name1="age-related.2.true.hypomethylation_15M_rep1_sb"
age.hypomethylation.15M <- tss[which((tss$c15-tss$c1)<=-.4 & abs(tss$csc15-tss$csc1)>=.05 & abs(tss$csc15-tss$csc1)<=.2),]
age.hypomethylation.15M.genes <- get.genelist(age.hypomethylation.15M, genelist.nam=genelist.name1)
head(age.hypomethylation.15M)
genelist.name2="csc-related.2.true.hypomethylation_15M_rep1_sb"
csc.hypomethylation.15M <- tss[which(abs(tss$c15-tss$c1)<=.05 & abs(tss$c15-tss$c1)<=.2 & (tss$csc15-tss$csc1)<=-.4),]
csc.hypomethylation.15M.genes <- get.genelist(csc.hypomethylation.15M, genelist.nam=genelist.name2)
head(csc.hypomethylation.15M)

length(age.hypomethylation.15M.genes)
length(csc.hypomethylation.15M.genes)
intersect_age.csc <- intersect(age.hypomethylation.15M.genes, csc.hypomethylation.15M.genes)
age.related.stable <- setdiff(age.hypomethylation.15M.genes, intersect_age.csc)
write.table(age.related.stable, paste(newdir,"/TSS_", genelist.name1, "_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.
csc.related.stable <- setdiff(csc.hypomethylation.15M.genes, intersect_age.csc)
write.table(csc.related.stable, paste(newdir,"/TSS_", genelist.name2, "_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.


# ######################################################################################################################################################################################
# # INTERMEDIATE HYPOMETHYLATION
# ######################################################################################################################################################################################
# 
# 
# #####
# # 6M HYPOMETHYLATION
# #####
# print("checkpoint3")
# genelist.name1="age-related.2.hypomethylation_6M_rep1_sb"
# age.hypomethylation.6M <- tss[which((tss$c6-tss$c1)<=-.4 & abs(tss$csc6-tss$csc1)<=.2),]
# if(dim(age.hypomethylation.6M)[1]>0){
#      age.hypomethylation.6M.genes <- get.genelist(age.hypomethylation.6M, genelist.nam=genelist.name1)
# }
# head(age.hypomethylation.6M)
# genelist.name2="csc-related.2.hypomethylation_6M_rep1_sb"
# csc.hypomethylation.6M <- tss[which(abs(tss$c6-tss$c1)<=.2 & (tss$csc6-tss$csc1)<=-.4),]
# if(dim(csc.hypomethylation.6M)[1]>0){
#      csc.hypomethylation.6M.genes <- get.genelist(csc.hypomethylation.6M, genelist.nam=genelist.name2)
# }
# head(csc.hypomethylation.6M)
# 
# length(age.hypomethylation.6M.genes)
# length(csc.hypomethylation.6M.genes)
# intersect_age.csc <- intersect(age.hypomethylation.6M.genes, csc.hypomethylation.6M.genes)
# age.related.stable <- setdiff(age.hypomethylation.6M.genes, intersect_age.csc)
# write.table(age.related.stable, paste(newdir,"/TSS_", genelist.name1, "_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.
# csc.related.stable <- setdiff(csc.hypomethylation.6M.genes, intersect_age.csc)
# write.table(csc.related.stable, paste(newdir,"/TSS_", genelist.name2, "_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.
# 
# #####
# # 10M HYPOMETHYLATION
# #####
# print("checkpoint3")
# genelist.name1="age-related.2.hypomethylation_10M_rep1_sb"
# age.hypomethylation.10M <- tss[which((tss$c10-tss$c1)<=-.4 & abs(tss$csc10-tss$csc1)<=.2),]
# if(dim(age.hypomethylation.10M)[1]>0){
#      age.hypomethylation.10M.genes <- get.genelist(age.hypomethylation.10M, genelist.nam=genelist.name1)
# }
# head(age.hypomethylation.10M)
# genelist.name2="csc-related.2.hypomethylation_10M_rep1_sb"
# csc.hypomethylation.10M <- tss[which(abs(tss$c10-tss$c1)<=.2 & (tss$csc10-tss$csc1)<=-.4),]
# if(dim(csc.hypomethylation.10M)[1]>0){
#      csc.hypomethylation.10M.genes <- get.genelist(csc.hypomethylation.10M, genelist.nam=genelist.name2)
# }
# head(csc.hypomethylation.10M)
# 
# length(age.hypomethylation.10M.genes)
# length(csc.hypomethylation.10M.genes)
# intersect_age.csc <- intersect(age.hypomethylation.10M.genes, csc.hypomethylation.10M.genes)
# age.related.stable <- setdiff(age.hypomethylation.10M.genes, intersect_age.csc)
# write.table(age.related.stable, paste(newdir,"/TSS_", genelist.name1, "_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.
# csc.related.stable <- setdiff(csc.hypomethylation.10M.genes, intersect_age.csc)
# write.table(csc.related.stable, paste(newdir,"/TSS_", genelist.name2, "_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.
# 
# #####
# # 15M HYPOMETHYLATION
# #####
# print("checkpoint4")
# genelist.name1="age-related.2.hypomethylation_15M_rep1_sb"
# age.hypomethylation.15M <- tss[which((tss$c15-tss$c1)<=-.4 & abs(tss$csc15-tss$csc1)<=.2),]
# age.hypomethylation.15M.genes <- get.genelist(age.hypomethylation.15M, genelist.nam=genelist.name1)
# head(age.hypomethylation.15M)
# genelist.name2="csc-related.2.hypomethylation_15M_rep1_sb"
# csc.hypomethylation.15M <- tss[which(abs(tss$c15-tss$c1)<=.2 & (tss$csc15-tss$csc1)<=-.4),]
# csc.hypomethylation.15M.genes <- get.genelist(csc.hypomethylation.15M, genelist.nam=genelist.name2)
# head(csc.hypomethylation.15M)
# 
# length(age.hypomethylation.15M.genes)
# length(csc.hypomethylation.15M.genes)
# intersect_age.csc <- intersect(age.hypomethylation.15M.genes, csc.hypomethylation.15M.genes)
# age.related.stable <- setdiff(age.hypomethylation.15M.genes, intersect_age.csc)
# write.table(age.related.stable, paste(newdir,"/TSS_", genelist.name1, "_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.
# csc.related.stable <- setdiff(csc.hypomethylation.15M.genes, intersect_age.csc)
# write.table(csc.related.stable, paste(newdir,"/TSS_", genelist.name2, "_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.
