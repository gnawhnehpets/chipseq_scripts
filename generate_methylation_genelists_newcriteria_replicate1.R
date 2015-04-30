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
newdir <- paste0(system.dir, "Michelle/MethylationData/methylatedGeneLists/", date)
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

#####
# 1M
#####

#Get CpG island probes within +/-1500 bp from TSS that get methylated/demethylated by more than 20% (age/proliferation dependent methylation), and then the probes that further get methylated/demethyalted upon CSC exposure.
#Get CpG island probes that get methylated specifically upon further exposure 
#abs() <= 0.05 will pick probes that do not undergo much changes in methylation upon culturing/ageing.
#treatment specific (de)methylation
genelist.name="age.hypermethylation_1M_rep1"
age.stable.1M.hypermethylation <- tss[which((tss$csc1-tss$c1)>=.2),]  #age.stable.1M.newdata
age.stable.1M.hypermethylation.genelist <- get.genelist(age.stable.1M.hypermethylation)
dim(age.stable.1M.hypermethylation)
genelist.name="age.hypomethylation_1M_rep1"
age.stable.1M.hypomethylation <- tss[which((tss$csc1-tss$c1)<=-.2),]  #age.stable.1M.newdata
age.stable.1M.hypomethylation.genelist <- get.genelist(age.stable.1M.hypomethylation)
dim(age.stable.1M.hypomethylation)

######################################################################################################################################################################################
######################################################################################################################################################################################
# STABLE
######################################################################################################################################################################################
######################################################################################################################################################################################

#####
# 6M
#####
getwd()
#age.stable, treatment.specific methylation
print("checkpoint1")
genelist.name="age.stable.hypermethylation_6M_rep1"
age.hypermethylation.6M.true <- rownames(tss[which((tss$c6-tss$c1)<=.05 & (tss$csc6-tss$csc1)>=.2),])
age.hypermethylation.6M.true.genelist <- get.only.genelist(age.hypermethylation.6M.true)
age.hypermethylation.6M.false <- rownames(tss[which((tss$c6-tss$c1)>=.05 & (tss$csc6-tss$csc1)>=.2),])
age.hypermethylation.6M.false.genelist <- get.only.genelist(age.hypermethylation.6M.false)
age.hypermethylation.6M.genelist <- setdiff(age.hypermethylation.6M.true.genelist, age.hypermethylation.6M.false.genelist)
write.table(age.hypermethylation.6M.genelist, paste(newdir,"/TSS_", genelist.name, "_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.
length(age.hypermethylation.6M.true.genelist); length(age.hypermethylation.6M.false.genelist); length(age.hypermethylation.6M.genelist)
age.hypermethyation.6M.probes <- setdiff(age.hypermethylation.6M.true, age.hypermethylation.6M.false)
age.hypermethyation.6M.beta <- tss[age.hypermethyation.6M.probes,]
if(dim(age.hypermethyation.6M.beta)[1]>0){
     dat <- get.genelist(age.hypermethyation.6M.beta)
     length(unique(dat)); print(genelist.name)
}else{
     print("No data: ", genelist.name)
}
print("checkpoint2")
genelist.name="age.stable.hypomethylation_6M_rep1"
age.hypomethylation.6M.true <- rownames(tss[which((tss$c6-tss$c1)<=.05 & (tss$csc6-tss$csc1)<=-.2),])
age.hypomethylation.6M.true.genelist <- get.only.genelist(age.hypomethylation.6M.true)
age.hypomethylation.6M.false <- rownames(tss[which((tss$c6-tss$c1)>=.05 & (tss$csc6-tss$csc1)<=-.2),])
age.hypomethylation.6M.false.genelist <- get.only.genelist(age.hypomethylation.6M.false)
age.hypomethylation.6M.genelist <- setdiff(age.hypomethylation.6M.true.genelist, age.hypomethylation.6M.false.genelist)
write.table(age.hypomethylation.6M.genelist, paste(newdir,"/TSS_", genelist.name, "_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.
length(age.hypomethylation.6M.true.genelist); length(age.hypomethylation.6M.false.genelist); length(age.hypomethylation.6M.genelist)
age.hypermethyation.6M.probes <- setdiff(age.hypomethylation.6M.true, age.hypomethylation.6M.false)
age.hypermethyation.6M.beta <- tss[age.hypermethyation.6M.probes,]
if(dim(age.hypermethyation.6M.beta)[1]>0){
     dat <- get.genelist(age.hypermethyation.6M.beta)
     length(unique(dat)); print(genelist.name)
}else{
     print("No data: ", genelist.name)
}
print("checkpoint3")
#treatment.stable, age.specific methylation
genelist.name="treatment.stable.hypermethylation_6M_rep1"
treatment.hypermethylation.6M.true <- rownames(tss[which((tss$c6-tss$c1)>=.2 & (tss$csc6-tss$csc1)<=.05),])
treatment.hypermethylation.6M.true.genelist <- get.only.genelist(treatment.hypermethylation.6M.true)
treatment.hypermethylation.6M.false <- rownames(tss[which((tss$c6-tss$c1)>=.2 & (tss$csc6-tss$csc1)>=.05),])
treatment.hypermethylation.6M.false.genelist <- get.only.genelist(treatment.hypermethylation.6M.false)
treatment.hypermethylation.6M.genelist <- setdiff(treatment.hypermethylation.6M.true.genelist, treatment.hypermethylation.6M.false.genelist)
write.table(treatment.hypermethylation.6M.genelist, paste(newdir,"/TSS_", genelist.name, "_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.
length(treatment.hypermethylation.6M.true.genelist); length(treatment.hypermethylation.6M.false.genelist); length(treatment.hypermethylation.6M.genelist)
treatment.hypermethyation.6M.probes <- setdiff(treatment.hypermethylation.6M.true, treatment.hypermethylation.6M.false)
treatment.hypermethyation.6M.beta <- tss[treatment.hypermethyation.6M.probes,]
if(dim(treatment.hypermethyation.6M.beta)[1]>0){
     dat <- get.genelist(treatment.hypermethyation.6M.beta)
     length(unique(dat)); print(genelist.name)
}else{
     print("No data: ", genelist.name)
}
print("checkpoint4")
genelist.name="treatment.stable.hypomethylation_6M_rep1"
treatment.hypomethylation.6M.true <- rownames(tss[which((tss$c6-tss$c1)<=-.2 & (tss$csc6-tss$csc1)<=.05),])
treatment.hypomethylation.6M.true.genelist <- get.only.genelist(treatment.hypomethylation.6M.true)
treatment.hypomethylation.6M.false <- rownames(tss[which((tss$c6-tss$c1)<=-.2 & (tss$csc6-tss$csc1)>=.05),])
treatment.hypomethylation.6M.false.genelist <- get.only.genelist(treatment.hypomethylation.6M.false)
treatment.hypomethylation.6M.genelist <- setdiff(treatment.hypomethylation.6M.true.genelist, treatment.hypomethylation.6M.false.genelist)
write.table(treatment.hypomethylation.6M.genelist, paste(newdir,"/TSS_", genelist.name, "_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.
length(treatment.hypomethylation.6M.true.genelist); length(treatment.hypomethylation.6M.false.genelist); length(treatment.hypomethylation.6M.genelist)
treatment.hypermethyation.6M.probes <- setdiff(treatment.hypomethylation.6M.true, treatment.hypomethylation.6M.false)
treatment.hypermethyation.6M.beta <- tss[treatment.hypermethyation.6M.probes,]
if(dim(treatment.hypermethyation.6M.beta)[1]>0){
     dat <- get.genelist(treatment.hypermethyation.6M.beta)
     length(unique(dat)); print(genelist.name)
}else{
     print("No data: ", genelist.name)
}

#####
# 10M
#####
getwd()
#age.stable, treatment.specific methylation
print("checkpoint5")
genelist.name="age.stable.hypermethylation_10M_rep1"
age.hypermethylation.10M.true <- rownames(tss[which((tss$c10-tss$c1)<=.05 & (tss$csc10-tss$csc1)>=.2),])
age.hypermethylation.10M.true.genelist <- get.only.genelist(age.hypermethylation.10M.true)
age.hypermethylation.10M.false <- rownames(tss[which((tss$c10-tss$c1)>=.05 & (tss$csc10-tss$csc1)>=.2),])
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
print("checkpoint6")
genelist.name="age.stable.hypomethylation_10M_rep1"
age.hypomethylation.10M.true <- rownames(tss[which((tss$c10-tss$c1)<=.05 & (tss$csc10-tss$csc1)<=-.2),])
age.hypomethylation.10M.true.genelist <- get.only.genelist(age.hypomethylation.10M.true)
age.hypomethylation.10M.false <- rownames(tss[which((tss$c10-tss$c1)>=.05 & (tss$csc10-tss$csc1)<=-.2),])
age.hypomethylation.10M.false.genelist <- get.only.genelist(age.hypomethylation.10M.false)
age.hypomethylation.10M.genelist <- setdiff(age.hypomethylation.10M.true.genelist, age.hypomethylation.10M.false.genelist)
write.table(age.hypomethylation.10M.genelist, paste(newdir,"/TSS_", genelist.name, "_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.
length(age.hypomethylation.10M.true.genelist); length(age.hypomethylation.10M.false.genelist); length(age.hypomethylation.10M.genelist)
age.hypermethyation.10M.probes <- setdiff(age.hypomethylation.10M.true, age.hypomethylation.10M.false)
age.hypermethyation.10M.beta <- tss[age.hypermethyation.10M.probes,]
if(dim(age.hypermethyation.10M.beta)[1]>0){
     dat <- get.genelist(age.hypermethyation.10M.beta)
     length(unique(dat)); print(genelist.name)
}else{
     print("No data: ", genelist.name)
}

#treatment.stable, age.specific methylation
print("checkpoint7")
genelist.name="treatment.stable.hypermethylation_10M_rep1"
treatment.hypermethylation.10M.true <- rownames(tss[which((tss$c10-tss$c1)>=.2 & (tss$csc10-tss$csc1)<=.05),])
treatment.hypermethylation.10M.true.genelist <- get.only.genelist(treatment.hypermethylation.10M.true)
treatment.hypermethylation.10M.false <- rownames(tss[which((tss$c10-tss$c1)>=.2 & (tss$csc10-tss$csc1)>=.05),])
treatment.hypermethylation.10M.false.genelist <- get.only.genelist(treatment.hypermethylation.10M.false)
treatment.hypermethylation.10M.genelist <- setdiff(treatment.hypermethylation.10M.true.genelist, treatment.hypermethylation.10M.false.genelist)
write.table(treatment.hypermethylation.10M.genelist, paste(newdir,"/TSS_", genelist.name, "_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.
length(treatment.hypermethylation.10M.true.genelist); length(treatment.hypermethylation.10M.false.genelist); length(treatment.hypermethylation.10M.genelist)
treatment.hypermethyation.10M.probes <- setdiff(treatment.hypermethylation.10M.true, treatment.hypermethylation.10M.false)
treatment.hypermethyation.10M.beta <- tss[treatment.hypermethyation.10M.probes,]
if(dim(treatment.hypermethyation.10M.beta)[1]>0){
     dat <- get.genelist(treatment.hypermethyation.10M.beta)
     length(unique(dat)); print(genelist.name)
}else{
     print("No data: ", genelist.name)
}
print("checkpoint8")
genelist.name="treatment.stable.hypomethylation_10M_rep1"
treatment.hypomethylation.10M.true <- rownames(tss[which((tss$c10-tss$c1)<=-.2 & (tss$csc10-tss$csc1)<=.05),])
treatment.hypomethylation.10M.true.genelist <- get.only.genelist(treatment.hypomethylation.10M.true)
treatment.hypomethylation.10M.false <- rownames(tss[which((tss$c10-tss$c1)<=-.2 & (tss$csc10-tss$csc1)>=.05),])
treatment.hypomethylation.10M.false.genelist <- get.only.genelist(treatment.hypomethylation.10M.false)
treatment.hypomethylation.10M.genelist <- setdiff(treatment.hypomethylation.10M.true.genelist, treatment.hypomethylation.10M.false.genelist)
write.table(treatment.hypomethylation.10M.genelist, paste(newdir,"/TSS_", genelist.name, "_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.
length(treatment.hypomethylation.10M.true.genelist); length(treatment.hypomethylation.10M.false.genelist); length(treatment.hypomethylation.10M.genelist)
treatment.hypermethyation.10M.probes <- setdiff(treatment.hypomethylation.10M.true, treatment.hypomethylation.10M.false)
treatment.hypermethyation.10M.beta <- tss[treatment.hypermethyation.10M.probes,]
if(dim(treatment.hypermethyation.10M.beta)[1]>0){
     dat <- get.genelist(treatment.hypermethyation.10M.beta)
     length(unique(dat)); print(genelist.name)
}else{
     print("No data: ", genelist.name)
}

#####
# 15M
#####
getwd()
#age.stable, treatment.specific methylation
print("checkpoint9")
genelist.name="age.stable.hypermethylation_15M_rep1"
age.hypermethylation.15M.true <- rownames(tss[which((tss$c15-tss$c1)<=.05 & (tss$csc15-tss$csc1)>=.2),])
age.hypermethylation.15M.true.genelist <- get.only.genelist(age.hypermethylation.15M.true)
age.hypermethylation.15M.false <- rownames(tss[which((tss$c15-tss$c1)>=.05 & (tss$csc15-tss$csc1)>=.2),])
age.hypermethylation.15M.false.genelist <- get.only.genelist(age.hypermethylation.15M.false)
age.hypermethylation.15M.genelist <- setdiff(age.hypermethylation.15M.true.genelist, age.hypermethylation.15M.false.genelist)
write.table(age.hypermethylation.15M.genelist, paste(newdir,"/TSS_", genelist.name, "_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.
length(age.hypermethylation.15M.true.genelist); length(age.hypermethylation.15M.false.genelist); length(age.hypermethylation.15M.genelist)
age.hypermethyation.15M.probes <- setdiff(age.hypermethylation.15M.true, age.hypermethylation.15M.false)
age.hypermethyation.15M.beta <- tss[age.hypermethyation.15M.probes,]
if(dim(age.hypermethyation.15M.beta)[1]>0){
     dat <- get.genelist(age.hypermethyation.15M.beta)
     length(unique(dat)); print(genelist.name)
}else{
     print("No data: ", genelist.name)
}
print("checkpoint10")
genelist.name="age.stable.hypomethylation_15M_rep1"
age.hypomethylation.15M.true <- rownames(tss[which((tss$c15-tss$c1)<=.05 & (tss$csc15-tss$csc1)<=-.2),])
age.hypomethylation.15M.true.genelist <- get.only.genelist(age.hypomethylation.15M.true)
age.hypomethylation.15M.false <- rownames(tss[which((tss$c15-tss$c1)>=.05 & (tss$csc15-tss$csc1)<=-.2),])
age.hypomethylation.15M.false.genelist <- get.only.genelist(age.hypomethylation.15M.false)
age.hypomethylation.15M.genelist <- setdiff(age.hypomethylation.15M.true.genelist, age.hypomethylation.15M.false.genelist)
write.table(age.hypomethylation.15M.genelist, paste(newdir,"/TSS_", genelist.name, "_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.
length(age.hypomethylation.15M.true.genelist); length(age.hypomethylation.15M.false.genelist); length(age.hypomethylation.15M.genelist)
age.hypermethyation.15M.probes <- setdiff(age.hypomethylation.15M.true, age.hypomethylation.15M.false)
age.hypermethyation.15M.beta <- tss[age.hypermethyation.15M.probes,]
if(dim(age.hypermethyation.15M.beta)[1]>0){
     dat <- get.genelist(age.hypermethyation.15M.beta)
     length(unique(dat)); print(genelist.name)
}else{
     print("No data: ", genelist.name)
}

#treatment.stable, age.specific methylation
print("checkpoint11")
genelist.name="treatment.stable.hypermethylation_15M_rep1"
treatment.hypermethylation.15M.true <- rownames(tss[which((tss$c15-tss$c1)>=.2 & (tss$csc15-tss$csc1)<=.05),])
treatment.hypermethylation.15M.true.genelist <- get.only.genelist(treatment.hypermethylation.15M.true)
treatment.hypermethylation.15M.false <- rownames(tss[which((tss$c15-tss$c1)>=.2 & (tss$csc15-tss$csc1)>=.05),])
treatment.hypermethylation.15M.false.genelist <- get.only.genelist(treatment.hypermethylation.15M.false)
treatment.hypermethylation.15M.genelist <- setdiff(treatment.hypermethylation.15M.true.genelist, treatment.hypermethylation.15M.false.genelist)
write.table(treatment.hypermethylation.15M.genelist, paste(newdir,"/TSS_", genelist.name, "_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.
length(treatment.hypermethylation.15M.true.genelist); length(treatment.hypermethylation.15M.false.genelist); length(treatment.hypermethylation.15M.genelist)
treatment.hypermethyation.15M.probes <- setdiff(treatment.hypermethylation.15M.true, treatment.hypermethylation.15M.false)
treatment.hypermethyation.15M.beta <- tss[treatment.hypermethyation.15M.probes,]
if(dim(treatment.hypermethyation.15M.beta)[1]>0){
     dat <- get.genelist(treatment.hypermethyation.15M.beta)
     length(unique(dat)); print(genelist.name)
}else{
     print("No data: ", genelist.name)
}

print("checkpoint12")
genelist.name="treatment.stable.hypomethylation_15M_rep1"
treatment.hypomethylation.15M.true <- rownames(tss[which((tss$c15-tss$c1)<=-.2 & (tss$csc15-tss$csc1)<=.05),])
treatment.hypomethylation.15M.true.genelist <- get.only.genelist(treatment.hypomethylation.15M.true)
treatment.hypomethylation.15M.false <- rownames(tss[which((tss$c15-tss$c1)<=-.2 & (tss$csc15-tss$csc1)>=.05),])
treatment.hypomethylation.15M.false.genelist <- get.only.genelist(treatment.hypomethylation.15M.false)
treatment.hypomethylation.15M.genelist <- setdiff(treatment.hypomethylation.15M.true.genelist, treatment.hypomethylation.15M.false.genelist)
write.table(treatment.hypomethylation.15M.genelist, paste(newdir,"/TSS_", genelist.name, "_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.
length(treatment.hypomethylation.15M.true.genelist); length(treatment.hypomethylation.15M.false.genelist); length(treatment.hypomethylation.15M.genelist)
treatment.hypermethyation.15M.probes <- setdiff(treatment.hypomethylation.15M.true, treatment.hypomethylation.15M.false)
treatment.hypermethyation.15M.beta <- tss[treatment.hypermethyation.15M.probes,]
if(dim(treatment.hypermethyation.15M.beta)[1]>0){
     dat <- get.genelist(treatment.hypermethyation.15M.beta)
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
# 6M
#####
getwd()
#age.intermediate, treatment.specific methylation
print("checkpoint13")
genelist.name="age.intermediate.hypermethylation_6M_rep1"
age.hypermethylation.6M.true <- rownames(tss[which((tss$c6-tss$c1)<=.2 & (tss$csc6-tss$csc1)>=.2),])
age.hypermethylation.6M.true.genelist <- get.only.genelist(age.hypermethylation.6M.true)
age.hypermethylation.6M.false <- rownames(tss[which((tss$c6-tss$c1)>=.2 & (tss$csc6-tss$csc1)>=.2),])
age.hypermethylation.6M.false.genelist <- get.only.genelist(age.hypermethylation.6M.false)
age.hypermethylation.6M.genelist <- setdiff(age.hypermethylation.6M.true.genelist, age.hypermethylation.6M.false.genelist)
write.table(age.hypermethylation.6M.genelist, paste(newdir,"/TSS_", genelist.name, "_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.
length(age.hypermethylation.6M.true.genelist); length(age.hypermethylation.6M.false.genelist); length(age.hypermethylation.6M.genelist)
age.hypermethyation.6M.probes <- setdiff(age.hypermethylation.6M.true, age.hypermethylation.6M.false)
age.hypermethyation.6M.beta <- tss[age.hypermethyation.6M.probes,]
if(dim(age.hypermethyation.6M.beta)[1]>0){
     dat <- get.genelist(age.hypermethyation.6M.beta)
     length(unique(dat)); print(genelist.name)
}else{
     print("No data: ", genelist.name)
}
print("checkpoint14")
genelist.name="age.intermediate.hypomethylation_6M_rep1"
age.hypomethylation.6M.true <- rownames(tss[which((tss$c6-tss$c1)<=.2 & (tss$csc6-tss$csc1)<=-.2),])
age.hypomethylation.6M.true.genelist <- get.only.genelist(age.hypomethylation.6M.true)
age.hypomethylation.6M.false <- rownames(tss[which((tss$c6-tss$c1)>=.2 & (tss$csc6-tss$csc1)<=-.2),])
age.hypomethylation.6M.false.genelist <- get.only.genelist(age.hypomethylation.6M.false)
age.hypomethylation.6M.genelist <- setdiff(age.hypomethylation.6M.true.genelist, age.hypomethylation.6M.false.genelist)
write.table(age.hypomethylation.6M.genelist, paste(newdir,"/TSS_", genelist.name, "_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.
length(age.hypomethylation.6M.true.genelist); length(age.hypomethylation.6M.false.genelist); length(age.hypomethylation.6M.genelist)
age.hypermethyation.6M.probes <- setdiff(age.hypomethylation.6M.true, age.hypomethylation.6M.false)
age.hypermethyation.6M.beta <- tss[age.hypermethyation.6M.probes,]
if(dim(age.hypermethyation.6M.beta)[1]>0){
     dat <- get.genelist(age.hypermethyation.6M.beta)
     length(unique(dat)); print(genelist.name)
}else{
     print("No data: ", genelist.name)
}

#treatment.intermediate, age.specific methylation
print("checkpoint15")
genelist.name="treatment.intermediate.hypermethylation_6M_rep1"
treatment.hypermethylation.6M.true <- rownames(tss[which((tss$c6-tss$c1)>=.2 & (tss$csc6-tss$csc1)<=.2),])
treatment.hypermethylation.6M.true.genelist <- get.only.genelist(treatment.hypermethylation.6M.true)
treatment.hypermethylation.6M.false <- rownames(tss[which((tss$c6-tss$c1)>=.2 & (tss$csc6-tss$csc1)>=.2),])
treatment.hypermethylation.6M.false.genelist <- get.only.genelist(treatment.hypermethylation.6M.false)
treatment.hypermethylation.6M.genelist <- setdiff(treatment.hypermethylation.6M.true.genelist, treatment.hypermethylation.6M.false.genelist)
write.table(treatment.hypermethylation.6M.genelist, paste(newdir,"/TSS_", genelist.name, "_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.
length(treatment.hypermethylation.6M.true.genelist); length(treatment.hypermethylation.6M.false.genelist); length(treatment.hypermethylation.6M.genelist)
treatment.hypermethyation.6M.probes <- setdiff(treatment.hypermethylation.6M.true, treatment.hypermethylation.6M.false)
treatment.hypermethyation.6M.beta <- tss[treatment.hypermethyation.6M.probes,]
if(dim(treatment.hypermethyation.6M.beta)[1]>0){
     dat <- get.genelist(treatment.hypermethyation.6M.beta)
     length(unique(dat)); print(genelist.name)
}else{
     print("No data: ", genelist.name)
}
print("checkpoint16")
genelist.name="treatment.intermediate.hypomethylation_6M_rep1"
treatment.hypomethylation.6M.true <- rownames(tss[which((tss$c6-tss$c1)<=-.2 & (tss$csc6-tss$csc1)<=.2),])
treatment.hypomethylation.6M.true.genelist <- get.only.genelist(treatment.hypomethylation.6M.true)
treatment.hypomethylation.6M.false <- rownames(tss[which((tss$c6-tss$c1)<=-.2 & (tss$csc6-tss$csc1)>=.2),])
treatment.hypomethylation.6M.false.genelist <- get.only.genelist(treatment.hypomethylation.6M.false)
treatment.hypomethylation.6M.genelist <- setdiff(treatment.hypomethylation.6M.true.genelist, treatment.hypomethylation.6M.false.genelist)
write.table(treatment.hypomethylation.6M.genelist, paste(newdir,"/TSS_", genelist.name, "_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.
length(treatment.hypomethylation.6M.true.genelist); length(treatment.hypomethylation.6M.false.genelist); length(treatment.hypomethylation.6M.genelist)
treatment.hypermethyation.6M.probes <- setdiff(treatment.hypomethylation.6M.true, treatment.hypomethylation.6M.false)
treatment.hypermethyation.6M.beta <- tss[treatment.hypermethyation.6M.probes,]
if(dim(treatment.hypermethyation.6M.beta)[1]>0){
     dat <- get.genelist(treatment.hypermethyation.6M.beta)
     length(unique(dat)); print(genelist.name)
}else{
     print("No data: ", genelist.name)
}

#####
# 10M
#####
getwd()
#age.intermediate, treatment.specific methylation
print("checkpoint17")
genelist.name="age.intermediate.hypermethylation_10M_rep1"
age.hypermethylation.10M.true <- rownames(tss[which((tss$c10-tss$c1)<=.2 & (tss$csc10-tss$csc1)>=.2),])
age.hypermethylation.10M.true.genelist <- get.only.genelist(age.hypermethylation.10M.true)
age.hypermethylation.10M.false <- rownames(tss[which((tss$c10-tss$c1)>=.2 & (tss$csc10-tss$csc1)>=.2),])
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

genelist.name="age.intermediate.hypomethylation_10M_rep1"
print("checkpoint18")
age.hypomethylation.10M.true <- rownames(tss[which((tss$c10-tss$c1)<=.2 & (tss$csc10-tss$csc1)<=-.2),])
age.hypomethylation.10M.true.genelist <- get.only.genelist(age.hypomethylation.10M.true)
age.hypomethylation.10M.false <- rownames(tss[which((tss$c10-tss$c1)>=.2 & (tss$csc10-tss$csc1)<=-.2),])
age.hypomethylation.10M.false.genelist <- get.only.genelist(age.hypomethylation.10M.false)
age.hypomethylation.10M.genelist <- setdiff(age.hypomethylation.10M.true.genelist, age.hypomethylation.10M.false.genelist)
write.table(age.hypomethylation.10M.genelist, paste(newdir,"/TSS_", genelist.name, "_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.
length(age.hypomethylation.10M.true.genelist); length(age.hypomethylation.10M.false.genelist); length(age.hypomethylation.10M.genelist)
age.hypermethyation.10M.probes <- setdiff(age.hypomethylation.10M.true, age.hypomethylation.10M.false)
age.hypermethyation.10M.beta <- tss[age.hypermethyation.10M.probes,]
if(dim(age.hypermethyation.10M.beta)[1]>0){
     dat <- get.genelist(age.hypermethyation.10M.beta)
     length(unique(dat)); print(genelist.name)
}else{
     print("No data: ", genelist.name)
}

#treatment.intermediate, age.specific methylation
genelist.name="treatment.intermediate.hypermethylation_10M_rep1"
print("checkpoint19")
treatment.hypermethylation.10M.true <- rownames(tss[which((tss$c10-tss$c1)>=.2 & (tss$csc10-tss$csc1)<=.2),])
treatment.hypermethylation.10M.true.genelist <- get.only.genelist(treatment.hypermethylation.10M.true)
treatment.hypermethylation.10M.false <- rownames(tss[which((tss$c10-tss$c1)>=.2 & (tss$csc10-tss$csc1)>=.2),])
treatment.hypermethylation.10M.false.genelist <- get.only.genelist(treatment.hypermethylation.10M.false)
treatment.hypermethylation.10M.genelist <- setdiff(treatment.hypermethylation.10M.true.genelist, treatment.hypermethylation.10M.false.genelist)
write.table(treatment.hypermethylation.10M.genelist, paste(newdir,"/TSS_", genelist.name, "_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.
length(treatment.hypermethylation.10M.true.genelist); length(treatment.hypermethylation.10M.false.genelist); length(treatment.hypermethylation.10M.genelist)
treatment.hypermethyation.10M.probes <- setdiff(treatment.hypermethylation.10M.true, treatment.hypermethylation.10M.false)
treatment.hypermethyation.10M.beta <- tss[treatment.hypermethyation.10M.probes,]
if(dim(treatment.hypermethyation.10M.beta)[1]>0){
     dat <- get.genelist(treatment.hypermethyation.10M.beta)
     length(unique(dat)); print(genelist.name)
}else{
     print("No data: ", genelist.name)
}

genelist.name="treatment.intermediate.hypomethylation_10M_rep1"
print("checkpoint20")
treatment.hypomethylation.10M.true <- rownames(tss[which((tss$c10-tss$c1)<=-.2 & (tss$csc10-tss$csc1)<=.2),])
treatment.hypomethylation.10M.true.genelist <- get.only.genelist(treatment.hypomethylation.10M.true)
treatment.hypomethylation.10M.false <- rownames(tss[which((tss$c10-tss$c1)<=-.2 & (tss$csc10-tss$csc1)>=.2),])
treatment.hypomethylation.10M.false.genelist <- get.only.genelist(treatment.hypomethylation.10M.false)
treatment.hypomethylation.10M.genelist <- setdiff(treatment.hypomethylation.10M.true.genelist, treatment.hypomethylation.10M.false.genelist)
write.table(treatment.hypomethylation.10M.genelist, paste(newdir,"/TSS_", genelist.name, "_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.
length(treatment.hypomethylation.10M.true.genelist); length(treatment.hypomethylation.10M.false.genelist); length(treatment.hypomethylation.10M.genelist)
treatment.hypermethyation.10M.probes <- setdiff(treatment.hypomethylation.10M.true, treatment.hypomethylation.10M.false)
treatment.hypermethyation.10M.beta <- tss[treatment.hypermethyation.10M.probes,]
if(dim(treatment.hypermethyation.10M.beta)[1]>0){
     dat <- get.genelist(treatment.hypermethyation.10M.beta)
     length(unique(dat)); print(genelist.name)
}else{
     print("No data: ", genelist.name)
}

#####
# 15M
#####
getwd()
#age.intermediate, treatment.specific methylation
genelist.name="age.intermediate.hypermethylation_15M_rep1"
print("checkpoint21")
age.hypermethylation.15M.true <- rownames(tss[which((tss$c15-tss$c1)<=.2 & (tss$csc15-tss$csc1)>=.2),])
age.hypermethylation.15M.true.genelist <- get.only.genelist(age.hypermethylation.15M.true)
age.hypermethylation.15M.false <- rownames(tss[which((tss$c15-tss$c1)>=.2 & (tss$csc15-tss$csc1)>=.2),])
age.hypermethylation.15M.false.genelist <- get.only.genelist(age.hypermethylation.15M.false)
age.hypermethylation.15M.genelist <- setdiff(age.hypermethylation.15M.true.genelist, age.hypermethylation.15M.false.genelist)
write.table(age.hypermethylation.15M.genelist, paste(newdir,"/TSS_", genelist.name, "_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.
length(age.hypermethylation.15M.true.genelist); length(age.hypermethylation.15M.false.genelist); length(age.hypermethylation.15M.genelist)
age.hypermethyation.15M.probes <- setdiff(age.hypermethylation.15M.true, age.hypermethylation.15M.false)
age.hypermethyation.15M.beta <- tss[age.hypermethyation.15M.probes,]
if(dim(age.hypermethyation.15M.beta)[1]>0){
     dat <- get.genelist(age.hypermethyation.15M.beta)
     length(unique(dat)); print(genelist.name)
}else{
     print("No data: ", genelist.name)
}

genelist.name="age.intermediate.hypomethylation_15M_rep1"
print("checkpoint22")
age.hypomethylation.15M.true <- rownames(tss[which((tss$c15-tss$c1)<=.2 & (tss$csc15-tss$csc1)<=-.2),])
age.hypomethylation.15M.true.genelist <- get.only.genelist(age.hypomethylation.15M.true)
age.hypomethylation.15M.false <- rownames(tss[which((tss$c15-tss$c1)>=.2 & (tss$csc15-tss$csc1)<=-.2),])
age.hypomethylation.15M.false.genelist <- get.only.genelist(age.hypomethylation.15M.false)
age.hypomethylation.15M.genelist <- setdiff(age.hypomethylation.15M.true.genelist, age.hypomethylation.15M.false.genelist)
write.table(age.hypomethylation.15M.genelist, paste(newdir,"/TSS_", genelist.name, "_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.
length(age.hypomethylation.15M.true.genelist); length(age.hypomethylation.15M.false.genelist); length(age.hypomethylation.15M.genelist)
age.hypermethyation.15M.probes <- setdiff(age.hypomethylation.15M.true, age.hypomethylation.15M.false)
age.hypermethyation.15M.beta <- tss[age.hypermethyation.15M.probes,]
if(dim(age.hypermethyation.15M.beta)[1]>0){
     dat <- get.genelist(age.hypermethyation.15M.beta)
     length(unique(dat)); print(genelist.name)
}else{
     print("No data: ", genelist.name)
}

#treatment.intermediate, age.specific methylation
genelist.name="treatment.intermediate.hypermethylation_15M_rep1"
print("checkpoint23")
treatment.hypermethylation.15M.true <- rownames(tss[which((tss$c15-tss$c1)>=.2 & (tss$csc15-tss$csc1)<=.2),])
treatment.hypermethylation.15M.true.genelist <- get.only.genelist(treatment.hypermethylation.15M.true)
treatment.hypermethylation.15M.false <- rownames(tss[which((tss$c15-tss$c1)>=.2 & (tss$csc15-tss$csc1)>=.2),])
treatment.hypermethylation.15M.false.genelist <- get.only.genelist(treatment.hypermethylation.15M.false)
treatment.hypermethylation.15M.genelist <- setdiff(treatment.hypermethylation.15M.true.genelist, treatment.hypermethylation.15M.false.genelist)
write.table(treatment.hypermethylation.15M.genelist, paste(newdir,"/TSS_", genelist.name, "_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.
length(treatment.hypermethylation.15M.true.genelist); length(treatment.hypermethylation.15M.false.genelist); length(treatment.hypermethylation.15M.genelist)
treatment.hypermethyation.15M.probes <- setdiff(treatment.hypermethylation.15M.true, treatment.hypermethylation.15M.false)
treatment.hypermethyation.15M.beta <- tss[treatment.hypermethyation.15M.probes,]
if(dim(treatment.hypermethyation.15M.beta)[1]>0){
     dat <- get.genelist(treatment.hypermethyation.15M.beta)
     length(unique(dat)); print(genelist.name)
}else{
     print("No data: ", genelist.name)
}

genelist.name="treatment.intermediate.hypomethylation_15M_rep1"
print("checkpoint24")
treatment.hypomethylation.15M.true <- rownames(tss[which((tss$c15-tss$c1)<=-.2 & (tss$csc15-tss$csc1)<=.2),])
treatment.hypomethylation.15M.true.genelist <- get.only.genelist(treatment.hypomethylation.15M.true)
treatment.hypomethylation.15M.false <- rownames(tss[which((tss$c15-tss$c1)<=-.2 & (tss$csc15-tss$csc1)>=.2),])
treatment.hypomethylation.15M.false.genelist <- get.only.genelist(treatment.hypomethylation.15M.false)
treatment.hypomethylation.15M.genelist <- setdiff(treatment.hypomethylation.15M.true.genelist, treatment.hypomethylation.15M.false.genelist)
write.table(treatment.hypomethylation.15M.genelist, paste(newdir,"/TSS_", genelist.name, "_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.
length(treatment.hypomethylation.15M.true.genelist); length(treatment.hypomethylation.15M.false.genelist); length(treatment.hypomethylation.15M.genelist)
treatment.hypermethyation.15M.probes <- setdiff(treatment.hypomethylation.15M.true, treatment.hypomethylation.15M.false)
treatment.hypermethyation.15M.beta <- tss[treatment.hypermethyation.15M.probes,]
if(dim(treatment.hypermethyation.15M.beta)[1]>0){
     dat <- get.genelist(treatment.hypermethyation.15M.beta)
     length(unique(dat)); print(genelist.name)
}else{
     print("No data: ", genelist.name)
}
