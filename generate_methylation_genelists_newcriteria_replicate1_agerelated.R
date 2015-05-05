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
print("CORRECTED")
# To generate genelists from raw idat files
# This particular script is specific to Michelle's untreated vs CSC-treated lung cells at 4 timepoints (1M, 6M, 10M, 15M)

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
#Load images created by Ash
load(paste0( system.dir, "/Michelle/ash/jhu-methylation-scripts/R_image/all.probes.hg19.RData"))
load(paste0( system.dir, "/Michelle/ash/jhu-methylation-scripts/R_image/cgi.xy.filter.RData"))
load(paste0( system.dir, "/Michelle/ash/jhu-methylation-scripts/R_image/probesWithin5000TSS.RData"))
load(paste0( system.dir, "/Michelle/ash/jhu-methylation-scripts/R_image/probe.to.gene.RData"))
load(paste0( system.dir, "/Michelle/ash/jhu-methylation-scripts/R_image/annotation.RData"))
#Load object, new beta tables generated with new data
load(file=paste0(system.dir, "/Michelle/Robjects/michelle_new-annotation.RData"))
#load functions
source(paste0( system.dir, "/Michelle/Rscripts/methylation450k_functions.R"))
date <- gsub(" \\d+\\:\\d+\\:\\d+.*", "", Sys.time())
newdir <- paste0(system.dir, "Michelle/MethylationData/methylatedGeneLists/", date, "/rep1_agerelated/")
newdir
# Create new working directory
dir.create(newdir)
setwd(newdir)

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
# #get +/- 1500 TSS probes
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
#   write.csv(beta.tab_CpGI.PlusMinus1500_TSS,paste(folder, "tss_probes_old.csv",sep="") , quote=F, row.names=T) #row.names set to T to get the probe ids.
# }
# 
# #save image
# save(beta.tab_CpGI,beta.tab_CpGI.PlusMinus1500_TSS,file=paste0(system.dir, "/Michelle/Robjects/michelle_new-annotation_old.RData"))
load(file=paste0(system.dir, "/Michelle/Robjects/michelle_new-annotation_old.RData"))

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

######################################################################################################################################################################################
######################################################################################################################################################################################
# STABLE
######################################################################################################################################################################################
######################################################################################################################################################################################
#####
# 10M
#####
getwd()
#age.stable, treatment.specific methylation
print("checkpoint5")
genelist.name="age.stable.hypermethylation_10M_rep1"
# (tss$c10-tss$c1)>=.4 & (tss$csc10-tss$csc1)<=.05
age.hypermethylation.10M.true <- rownames(tss[which((tss$c10-tss$c1)>=.4 & (tss$csc10-tss$csc1)<=.05),])
age.hypermethylation.10M.true.genelist <- get.only.genelist(age.hypermethylation.10M.true)
age.hypermethylation.10M.false <- rownames(tss[which((tss$c10-tss$c1)>=.4 & (tss$csc10-tss$csc1)>=.05),])
age.hypermethylation.10M.false.genelist <- get.only.genelist(age.hypermethylation.10M.false)
age.hypermethylation.10M.genelist <- setdiff(age.hypermethylation.10M.true.genelist, age.hypermethylation.10M.false.genelist)
write.table(age.hypermethylation.10M.genelist, paste(newdir,"/TSS_", genelist.name, "_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.
length(age.hypermethylation.10M.true.genelist); length(age.hypermethylation.10M.false.genelist); length(age.hypermethylation.10M.genelist)
age.hypermethyation.10M.probes <- setdiff(age.hypermethylation.10M.true, age.hypermethylation.10M.false)
age.hypermethyation.10M.beta <- tss[age.hypermethyation.10M.probes,]
if(dim(age.hypermethyation.10M.beta)[1]>0){
     dat <- get.genelist(age.hypermethyation.10M.beta)
     length(unique(dat)); print(genelist.name)
}else{
     print("No data: ", genelist.name)
}

######################################################################################################################################################################################
######################################################################################################################################################################################
# INTERMEDIATE
######################################################################################################################################################################################
######################################################################################################################################################################################

#####
# 10M
#####
getwd()
#age.intermediate, treatment.specific methylation
print("checkpoint17")
genelist.name="age.intermediate.hypermethylation_10M_rep1"
age.hypermethylation.10M.true <- rownames(tss[which((tss$c10-tss$c1)>=.4 & (tss$csc10-tss$csc1)<=.2),])
age.hypermethylation.10M.true.genelist <- get.only.genelist(age.hypermethylation.10M.true)
age.hypermethylation.10M.false <- rownames(tss[which((tss$c10-tss$c1)>=.4 & (tss$csc10-tss$csc1)>=.2),])
age.hypermethylation.10M.false.genelist <- get.only.genelist(age.hypermethylation.10M.false)
age.hypermethylation.10M.genelist <- setdiff(age.hypermethylation.10M.true.genelist, age.hypermethylation.10M.false.genelist)
write.table(age.hypermethylation.10M.genelist, paste(newdir,"/TSS_", genelist.name, "_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.
length(age.hypermethylation.10M.true.genelist); length(age.hypermethylation.10M.false.genelist); length(age.hypermethylation.10M.genelist)
age.hypermethyation.10M.probes <- setdiff(age.hypermethylation.10M.true, age.hypermethylation.10M.false)
age.hypermethyation.10M.beta <- tss[age.hypermethyation.10M.probes,]
if(dim(age.hypermethyation.10M.beta)[1]>0){
     dat <- get.genelist(age.hypermethyation.10M.beta)
     length(unique(dat)); print(genelist.name)
}else{
     print("No data: ", genelist.name)
}
