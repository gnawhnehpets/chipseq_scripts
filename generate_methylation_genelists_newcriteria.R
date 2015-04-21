##########
# Purpose:
##########

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
# idat <- paste0( system.dir, "/users/shwang26/Michelle/MethylationData/idat")
# sheet_genelist <- paste0( system.dir, "/users/shwang26/Michelle/ash/sheet_genelist.csv")
# sheet_deltaplot <- paste0( system.dir, "/users/shwang26/Michelle/ash/sheet_deltaplot.csv")
# folder <- paste0( system.dir, "/users/shwang26/Michelle/ash/new_iteration/")
# 
# #Fuction for reading and normalizing data
# beta.tab <- data_methylation(idat, folder)
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
#   write.csv(beta.tab_CpGI.PlusMinus1500_TSS,paste(folder, "tss_probes.csv",sep="") , quote=F, row.names=T) #row.names set to T to get the probe ids.
# }
# 
# #save image
# save(beta.tab_CpGI,beta.tab_CpGI.PlusMinus1500_TSS,file=paste0(system.dir, "/Michelle/Robjects/michelle_new-annotation.RData"))
# load(file=paste0(system.dir, "/Michelle/Robjects/michelle_new-annotation.RData"))

# transition
subdir <- "combined"
timepoint="combined"
type <- "TSS"

head(beta.tab)
dim(beta.tab)
head(beta.tab_CpGI.PlusMinus1500_TSS)

tss <- beta.tab_CpGI.PlusMinus1500_TSS
tss_reorder <- cbind.data.frame(tss$C.1M, tss$CSC.1M, tss$C.6M, tss$CSC.6M, tss$C.10M, tss$CSC.10M, tss$C.15M, tss$CSC.15M)
colnames(tss_reorder) <- c("c1", "csc1", "c6", "csc6", "c10", "csc10", "c15", "csc15")
tss <- tss_reorder
colnames(tss)
rownames(tss) <- rownames(beta.tab_CpGI.PlusMinus1500_TSS)
rownames(tss)
cbind(tss$c1, tss$csc1)
#####
# 1M
#####

#Get CpG island probes within +/-1500 bp from TSS that get methylated/demethylated by more than 20% (age/proliferation dependent methylation), and then the probes that further get methylated/demethyalted upon CSC exposure.
#Get CpG island probes that get methylated specifically upon further exposure 
#abs() <= 0.05 will pick probes that do not undergo much changes in methylation upon culturing/ageing.
#treatment specific (de)methylation
genelist.name="treatment.specific.hypermethylation_1M"
stable.1M.treatment.hypermethylation <- tss[which((tss$csc1-tss$c1)>=.2),]  #stable.1M.newdata
stable.1M.treatment.hypermethylation.genelist <- get.genelist(stable.1M.treatment.hypermethylation)
dim(stable.1M.treatment.hypermethylation)
genelist.name="treatment.specific.hypomethylation_1M"
stable.1M.treatment.hypomethylation <- tss[which((tss$csc1-tss$c1)<=-.2),]  #stable.1M.newdata 
stable.1M.treatment.hypomethylation.genelist <- get.genelist(stable.1M.treatment.hypomethylation)
dim(stable.1M.treatment.hypomethylation)

#####
# 6M
#####
getwd()
#age specific (de)methylation
genelist.name="age.specific.hypermethylation_6M"
age.specific.hypermethylation.6M.true <- rownames(tss[which((tss$c6-tss$c1)>=.2 & (tss$csc6-tss$csc1)<=.2),])
age.specific.hypermethylation.6M.true.genelist <- get.only.genelist(age.specific.hypermethylation.6M.true)
age.specific.hypermethylation.6M.false <- rownames(tss[which((tss$c6-tss$c1)>=.2 & (tss$csc6-tss$csc1)>=.2),])
age.specific.hypermethylation.6M.false.genelist <- get.only.genelist(age.specific.hypermethylation.6M.false)
age.specific.hypermethylation.6M.genelist <- setdiff(age.specific.hypermethylation.6M.true.genelist, age.specific.hypermethylation.6M.false.genelist)
write.table(age.specific.hypermethylation.6M.genelist, paste(newdir,"/TSS_", genelist.name, "_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.
length(age.specific.hypermethylation.6M.true.genelist); length(age.specific.hypermethylation.6M.false.genelist); length(age.specific.hypermethylation.6M.genelist)
age.specific.hypermethyation.6M.probes <- setdiff(age.specific.hypermethylation.6M.true, age.specific.hypermethylation.6M.false)
age.specific.hypermethyation.6M.beta <- tss[age.specific.hypermethyation.6M.probes,]
dat <- get.genelist(age.specific.hypermethyation.6M.beta)
length(unique(dat))
head(age.specific.hypermethyation.6M.beta)


genelist.name="age.specific.hypomethylation_6M"
age.specific.hypomethylation.6M.true <- rownames(tss[which((tss$c6-tss$c1)<=-.2 & (tss$csc6-tss$csc1)<=.2),])
age.specific.hypomethylation.6M.true.genelist <- get.only.genelist(age.specific.hypomethylation.6M.true)
age.specific.hypomethylation.6M.false <- rownames(tss[which((tss$c6-tss$c1)<=-.2 & (tss$csc6-tss$csc1)>=.2),])
age.specific.hypomethylation.6M.false.genelist <- get.only.genelist(age.specific.hypomethylation.6M.false)
age.specific.hypomethylation.6M.genelist <- setdiff(age.specific.hypomethylation.6M.true.genelist, age.specific.hypomethylation.6M.false.genelist)
write.table(age.specific.hypomethylation.6M.genelist, paste(newdir,"/TSS_", genelist.name, "_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.
length(age.specific.hypomethylation.6M.true.genelist); length(age.specific.hypomethylation.6M.false.genelist); length(age.specific.hypomethylation.6M.genelist)
age.specific.hypomethyation.6M.probes <- setdiff(age.specific.hypomethylation.6M.true, age.specific.hypomethylation.6M.false)
age.specific.hypomethyation.6M.beta <- tss[age.specific.hypomethyation.6M.probes,]
dat <- get.genelist(age.specific.hypomethyation.6M.beta)
length(unique(dat))
head(age.specific.hypomethyation.6M.beta)

#treatment specific (de)methylation
genelist.name="treatment.specific.hypermethylation_6M"
treatment.specific.hypermethylation.6M.true <- rownames(tss[which((tss$c6-tss$c1)<=.2 & (tss$csc6-tss$csc1)>=.2),])
treatment.specific.hypermethylation.6M.true.genelist <- get.only.genelist(treatment.specific.hypermethylation.6M.true)
treatment.specific.hypermethylation.6M.false <- rownames(tss[which((tss$c6-tss$c1)>=.2 & (tss$csc6-tss$csc1)>=.2),])
treatment.specific.hypermethylation.6M.false.genelist <- get.only.genelist(treatment.specific.hypermethylation.6M.false)
treatment.specific.hypermethylation.6M.genelist <- setdiff(treatment.specific.hypermethylation.6M.true.genelist, treatment.specific.hypermethylation.6M.false.genelist)
write.table(treatment.specific.hypermethylation.6M.genelist, paste(newdir,"/TSS_", genelist.name, "_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.
length(treatment.specific.hypermethylation.6M.true.genelist); length(treatment.specific.hypermethylation.6M.false.genelist); length(treatment.specific.hypermethylation.6M.genelist)
treatment.specific.hypermethyation.6M.probes <- setdiff(treatment.specific.hypermethylation.6M.true, treatment.specific.hypermethylation.6M.false)
treatment.specific.hypermethyation.6M.beta <- tss[treatment.specific.hypermethyation.6M.probes,]
dat <- get.genelist(treatment.specific.hypermethyation.6M.beta)
length(unique(dat))
head(treatment.specific.hypermethyation.6M.beta)

genelist.name="treatment.specific.hypomethylation_6M"
treatment.specific.hypomethylation.6M.true <- rownames(tss[which((tss$c6-tss$c1)<=.2 & (tss$csc6-tss$csc1)<=-.2),])
treatment.specific.hypomethylation.6M.true.genelist <- get.only.genelist(treatment.specific.hypomethylation.6M.true)
treatment.specific.hypomethylation.6M.false <- rownames(tss[which((tss$c6-tss$c1)>=.2 & (tss$csc6-tss$csc1)<=-.2),])
treatment.specific.hypomethylation.6M.false.genelist <- get.only.genelist(treatment.specific.hypomethylation.6M.false)
treatment.specific.hypomethylation.6M.genelist <- setdiff(treatment.specific.hypomethylation.6M.true.genelist, treatment.specific.hypomethylation.6M.false.genelist)
write.table(treatment.specific.hypomethylation.6M.genelist, paste(newdir,"/TSS_", genelist.name, "_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.
length(treatment.specific.hypomethylation.6M.true.genelist); length(treatment.specific.hypomethylation.6M.false.genelist); length(treatment.specific.hypomethylation.6M.genelist)
treatment.specific.hypomethyation.6M.probes <- setdiff(treatment.specific.hypomethylation.6M.true, treatment.specific.hypomethylation.6M.false)
treatment.specific.hypomethyation.6M.beta <- tss[treatment.specific.hypomethyation.6M.probes,]
dat <- get.genelist(treatment.specific.hypomethyation.6M.beta)
length(unique(dat))
head(treatment.specific.hypomethyation.6M.beta)

#####
# 10M
#####
#age specific (de)methylation
genelist.name="age.specific.hypermethylation_10M"
age.specific.hypermethylation.10M.true <- rownames(tss[which((tss$c10-tss$c1)>=.2 & (tss$csc10-tss$csc1)<=.2),])
age.specific.hypermethylation.10M.true.genelist <- get.only.genelist(age.specific.hypermethylation.10M.true)
age.specific.hypermethylation.10M.false <- rownames(tss[which((tss$c10-tss$c1)>=.2 & (tss$csc10-tss$csc1)>=.2),])
age.specific.hypermethylation.10M.false.genelist <- get.only.genelist(age.specific.hypermethylation.10M.false)
age.specific.hypermethylation.10M.genelist <- setdiff(age.specific.hypermethylation.10M.true.genelist, age.specific.hypermethylation.10M.false.genelist)
write.table(age.specific.hypermethylation.10M.genelist, paste(newdir,"/TSS_", genelist.name, "_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.
length(age.specific.hypermethylation.10M.true.genelist); length(age.specific.hypermethylation.10M.false.genelist); length(age.specific.hypermethylation.10M.genelist)
age.specific.hypermethyation.10M.probes <- setdiff(age.specific.hypermethylation.10M.true, age.specific.hypermethylation.10M.false)
age.specific.hypermethyation.10M.beta <- tss[age.specific.hypermethyation.10M.probes,]
dat <- get.genelist(age.specific.hypermethyation.10M.beta)
length(unique(dat))
head(age.specific.hypermethyation.10M.beta)


genelist.name="age.specific.hypomethylation_10M"
age.specific.hypomethylation.10M.true <- rownames(tss[which((tss$c10-tss$c1)<=-.2 & (tss$csc10-tss$csc1)<=.2),])
age.specific.hypomethylation.10M.true.genelist <- get.only.genelist(age.specific.hypomethylation.10M.true)
age.specific.hypomethylation.10M.false <- rownames(tss[which((tss$c10-tss$c1)<=-.2 & (tss$csc10-tss$csc1)>=.2),])
age.specific.hypomethylation.10M.false.genelist <- get.only.genelist(age.specific.hypomethylation.10M.false)
age.specific.hypomethylation.10M.genelist <- setdiff(age.specific.hypomethylation.10M.true.genelist, age.specific.hypomethylation.10M.false.genelist)
write.table(age.specific.hypomethylation.10M.genelist, paste(newdir,"/TSS_", genelist.name, "_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.
length(age.specific.hypomethylation.10M.true.genelist); length(age.specific.hypomethylation.10M.false.genelist); length(age.specific.hypomethylation.10M.genelist)
age.specific.hypomethyation.10M.probes <- setdiff(age.specific.hypomethylation.10M.true, age.specific.hypomethylation.10M.false)
age.specific.hypomethyation.10M.beta <- tss[age.specific.hypomethyation.10M.probes,]
dat <- get.genelist(age.specific.hypomethyation.10M.beta)
length(unique(dat))
head(age.specific.hypomethyation.10M.beta)

#treatment specific (de)methylation
genelist.name="treatment.specific.hypermethylation_10M"
treatment.specific.hypermethylation.10M.true <- rownames(tss[which((tss$c10-tss$c1)<=.2 & (tss$csc10-tss$csc1)>=.2),])
treatment.specific.hypermethylation.10M.true.genelist <- get.only.genelist(treatment.specific.hypermethylation.10M.true)
treatment.specific.hypermethylation.10M.false <- rownames(tss[which((tss$c10-tss$c1)>=.2 & (tss$csc10-tss$csc1)>=.2),])
treatment.specific.hypermethylation.10M.false.genelist <- get.only.genelist(treatment.specific.hypermethylation.10M.false)
treatment.specific.hypermethylation.10M.genelist <- setdiff(treatment.specific.hypermethylation.10M.true.genelist, treatment.specific.hypermethylation.10M.false.genelist)
write.table(treatment.specific.hypermethylation.10M.genelist, paste(newdir,"/TSS_", genelist.name, "_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.
length(treatment.specific.hypermethylation.10M.true.genelist); length(treatment.specific.hypermethylation.10M.false.genelist); length(treatment.specific.hypermethylation.10M.genelist)
treatment.specific.hypermethyation.10M.probes <- setdiff(treatment.specific.hypermethylation.10M.true, treatment.specific.hypermethylation.10M.false)
treatment.specific.hypermethyation.10M.beta <- tss[treatment.specific.hypermethyation.10M.probes,]
dat <- get.genelist(treatment.specific.hypermethyation.10M.beta)
length(unique(dat))
head(treatment.specific.hypermethyation.10M.beta)

genelist.name="treatment.specific.hypomethylation_10M"
treatment.specific.hypomethylation.10M.true <- rownames(tss[which((tss$c10-tss$c1)<=.2 & (tss$csc10-tss$csc1)<=-.2),])
treatment.specific.hypomethylation.10M.true.genelist <- get.only.genelist(treatment.specific.hypomethylation.10M.true)
treatment.specific.hypomethylation.10M.false <- rownames(tss[which((tss$c10-tss$c1)>=.2 & (tss$csc10-tss$csc1)<=-.2),])
treatment.specific.hypomethylation.10M.false.genelist <- get.only.genelist(treatment.specific.hypomethylation.10M.false)
treatment.specific.hypomethylation.10M.genelist <- setdiff(treatment.specific.hypomethylation.10M.true.genelist, treatment.specific.hypomethylation.10M.false.genelist)
write.table(treatment.specific.hypomethylation.10M.genelist, paste(newdir,"/TSS_", genelist.name, "_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.
length(treatment.specific.hypomethylation.10M.true.genelist); length(treatment.specific.hypomethylation.10M.false.genelist); length(treatment.specific.hypomethylation.10M.genelist)
treatment.specific.hypomethyation.10M.probes <- setdiff(treatment.specific.hypomethylation.10M.true, treatment.specific.hypomethylation.10M.false)
treatment.specific.hypomethyation.10M.beta <- tss[treatment.specific.hypomethyation.10M.probes,]
dat <- get.genelist(treatment.specific.hypomethyation.10M.beta)
length(unique(dat))
head(treatment.specific.hypomethyation.10M.beta)

#####
# 15M
#####

#age specific (de)methylation
genelist.name="age.specific.hypermethylation_15M"
age.specific.hypermethylation.15M.true <- rownames(tss[which((tss$c15-tss$c1)>=.2 & (tss$csc15-tss$csc1)<=.2),])
age.specific.hypermethylation.15M.true.genelist <- get.only.genelist(age.specific.hypermethylation.15M.true)
age.specific.hypermethylation.15M.false <- rownames(tss[which((tss$c15-tss$c1)>=.2 & (tss$csc15-tss$csc1)>=.2),])
age.specific.hypermethylation.15M.false.genelist <- get.only.genelist(age.specific.hypermethylation.15M.false)
age.specific.hypermethylation.15M.genelist <- setdiff(age.specific.hypermethylation.15M.true.genelist, age.specific.hypermethylation.15M.false.genelist)
write.table(age.specific.hypermethylation.15M.genelist, paste(newdir,"/TSS_", genelist.name, "_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.
length(age.specific.hypermethylation.15M.true.genelist); length(age.specific.hypermethylation.15M.false.genelist); length(age.specific.hypermethylation.15M.genelist)
age.specific.hypermethyation.15M.probes <- setdiff(age.specific.hypermethylation.15M.true, age.specific.hypermethylation.15M.false)
age.specific.hypermethyation.15M.beta <- tss[age.specific.hypermethyation.15M.probes,]
dat <- get.genelist(age.specific.hypermethyation.15M.beta)
length(unique(dat))
head(age.specific.hypermethyation.15M.beta)


genelist.name="age.specific.hypomethylation_15M"
age.specific.hypomethylation.15M.true <- rownames(tss[which((tss$c15-tss$c1)<=-.2 & (tss$csc15-tss$csc1)<=.2),])
age.specific.hypomethylation.15M.true.genelist <- get.only.genelist(age.specific.hypomethylation.15M.true)
age.specific.hypomethylation.15M.false <- rownames(tss[which((tss$c15-tss$c1)<=-.2 & (tss$csc15-tss$csc1)>=.2),])
age.specific.hypomethylation.15M.false.genelist <- get.only.genelist(age.specific.hypomethylation.15M.false)
age.specific.hypomethylation.15M.genelist <- setdiff(age.specific.hypomethylation.15M.true.genelist, age.specific.hypomethylation.15M.false.genelist)
write.table(age.specific.hypomethylation.15M.genelist, paste(newdir,"/TSS_", genelist.name, "_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.
length(age.specific.hypomethylation.15M.true.genelist); length(age.specific.hypomethylation.15M.false.genelist); length(age.specific.hypomethylation.15M.genelist)
age.specific.hypomethyation.15M.probes <- setdiff(age.specific.hypomethylation.15M.true, age.specific.hypomethylation.15M.false)
age.specific.hypomethyation.15M.beta <- tss[age.specific.hypomethyation.15M.probes,]
dat <- get.genelist(age.specific.hypomethyation.15M.beta)
length(unique(dat))
head(age.specific.hypomethyation.15M.beta)

#treatment specific (de)methylation
genelist.name="treatment.specific.hypermethylation_15M"
treatment.specific.hypermethylation.15M.true <- rownames(tss[which((tss$c15-tss$c1)<=.2 & (tss$csc15-tss$csc1)>=.2),])
treatment.specific.hypermethylation.15M.true.genelist <- get.only.genelist(treatment.specific.hypermethylation.15M.true)
treatment.specific.hypermethylation.15M.false <- rownames(tss[which((tss$c15-tss$c1)>=.2 & (tss$csc15-tss$csc1)>=.2),])
treatment.specific.hypermethylation.15M.false.genelist <- get.only.genelist(treatment.specific.hypermethylation.15M.false)
treatment.specific.hypermethylation.15M.genelist <- setdiff(treatment.specific.hypermethylation.15M.true.genelist, treatment.specific.hypermethylation.15M.false.genelist)
write.table(treatment.specific.hypermethylation.15M.genelist, paste(newdir,"/TSS_", genelist.name, "_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.
length(treatment.specific.hypermethylation.15M.true.genelist); length(treatment.specific.hypermethylation.15M.false.genelist); length(treatment.specific.hypermethylation.15M.genelist)
treatment.specific.hypermethyation.15M.probes <- setdiff(treatment.specific.hypermethylation.15M.true, treatment.specific.hypermethylation.15M.false)
treatment.specific.hypermethyation.15M.beta <- tss[treatment.specific.hypermethyation.15M.probes,]
dat <- get.genelist(treatment.specific.hypermethyation.15M.beta)
length(unique(dat))
head(treatment.specific.hypermethyation.15M.beta)

genelist.name="treatment.specific.hypomethylation_15M"
treatment.specific.hypomethylation.15M.true <- rownames(tss[which((tss$c15-tss$c1)<=.2 & (tss$csc15-tss$csc1)<=-.2),])
treatment.specific.hypomethylation.15M.true.genelist <- get.only.genelist(treatment.specific.hypomethylation.15M.true)
treatment.specific.hypomethylation.15M.false <- rownames(tss[which((tss$c15-tss$c1)>=.2 & (tss$csc15-tss$csc1)<=-.2),])
treatment.specific.hypomethylation.15M.false.genelist <- get.only.genelist(treatment.specific.hypomethylation.15M.false)
treatment.specific.hypomethylation.15M.genelist <- setdiff(treatment.specific.hypomethylation.15M.true.genelist, treatment.specific.hypomethylation.15M.false.genelist)
write.table(treatment.specific.hypomethylation.15M.genelist, paste(newdir,"/TSS_", genelist.name, "_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.
length(treatment.specific.hypomethylation.15M.true.genelist); length(treatment.specific.hypomethylation.15M.false.genelist); length(treatment.specific.hypomethylation.15M.genelist)
treatment.specific.hypomethyation.15M.probes <- setdiff(treatment.specific.hypomethylation.15M.true, treatment.specific.hypomethylation.15M.false)
treatment.specific.hypomethyation.15M.beta <- tss[treatment.specific.hypomethyation.15M.probes,]
dat <- get.genelist(treatment.specific.hypomethyation.15M.beta)
length(unique(dat))
head(treatment.specific.hypomethyation.15M.beta)

#filtered beta.tab matrix 
get.genelist <- function(x, system.di=system.dir, probe.to.gen=probe.to.gene, newdi=newdir, typ=type, genelist.nam=genelist.name){
          if (length(rownames(x)) > 0) {
               x <- cbind.data.frame(x, "Gene"=rep(NA, times=nrow(x)))
               for(j in 1:nrow(x)) {
                    probe.id <- rownames(x)[j]
                    gene <- probe.to.gen$gene[probe.to.gen$probe==probe.id]
                    if(length(gene) == 1) { x$Gene[j] <- as.character(gene) }
               }
               write.table(x, paste(newdi,"/", typ, "_", genelist.nam, "_beta.values.txt",sep=""), quote=F, row.names=T, col.names=NA, sep="\t") #row.names set to T to get the probe ids.
          }     
          print("beta values written!")
          tmp <- as.matrix(sort(unique(unlist(strsplit(x$Gene, ";")))))
          write.table(tmp, paste(newdi,"/", typ, "_", genelist.nam, "_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.
     return(tmp)
}

#probes aka rownames() only
get.only.genelist <- function(x, probe.to.gen=probe.to.gene, newdi=newdir, typ=type, genelist.nam=genelist.name){
     list.of.genes <- vector()
     for(i in 1:length(x)) {
          probe.id <- x[i]
          gene <- unlist(strsplit(as.character(probe.to.gen$gene[probe.to.gen$probe==probe.id]), ";"))
          if(length(gene) == 1) { list.of.genes <- c(list.of.genes, gene) }
     }
     list.of.genes <- sort(unique(list.of.genes))
#      write.table(list.of.genes, paste(newdi,"/", typ, "_", genelist.nam, "_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.
     return(list.of.genes)
}