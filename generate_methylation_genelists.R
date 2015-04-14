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
genelist.name="treatment.specific.hypermethylation_stable.1M.newdata"
stable.1M.treatment.hypermethylation <- tss[which((tss$csc1-tss$c1)>=.2),]  #stable.1M.newdata
stable.1M.treatment.hypermethylation.genelist <- get.genelist(stable.1M.treatment.hypermethylation)
dim(stable.1M.treatment.hypermethylation)
genelist.name="treatment.specific.hypomethylation_stable.1M.newdata"
stable.1M.treatment.hypomethylation <- tss[which((tss$csc1-tss$c1)<=-.2),]  #stable.1M.newdata 
stable.1M.treatment.hypomethylation.genelist <- get.genelist(stable.1M.treatment.hypomethylation)
dim(stable.1M.treatment.hypomethylation)


#####
# 6M
#####

#Get CpG island probes within +/-1500 bp from TSS that get methylated/demethylated by more than 20% (age/proliferation dependent methylation), and then the probes that further get methylated/demethyalted upon CSC exposure.
#Get CpG island probes that get methylated specifically upon further exposure 
#abs() <= 0.05 will pick probes that do not undergo much changes in methylation upon culturing/ageing.
#treatment specific (de)methylation
genelist.name="treatment.specific.hypermethylation_stable.6M.newdata"
stable.6M.treatment.hypermethylation <- tss[which(abs(tss$c6-tss$c1)<=.05 & (tss$csc6-tss$c6)>=.2),]  #stable.6M.newdata
stable.6M.treatment.hypermethylation.genelist <- get.genelist(stable.6M.treatment.hypermethylation)
dim(stable.6M.treatment.hypermethylation)
genelist.name="treatment.specific.hypomethylation_stable.6M.newdata"
stable.6M.treatment.hypomethylation <- tss[which(abs(tss$c6-tss$c1)<=.05 & (tss$csc6-tss$c6)<=-.2),]  #stable.6M.newdata 
stable.6M.treatment.hypomethylation.genelist <- get.genelist(stable.6M.treatment.hypomethylation)
dim(stable.6M.treatment.hypomethylation)

#Get probes that range in methylation between +/-0.05 to +/-0.2 during aging/proliferation but are methylated/demethylated >= 0.2 
genelist.name="treatment.specific.hypermethylation_intermediate.6M.newdata"
intermediate.6M.treatment.hypermethylation <- tss[which(abs(tss$c6-tss$c1)>=.05 & abs(tss$c6-tss$c1)<=.2 & (tss$csc6-tss$c6)>=.2),]  #stable.6M.newdata
intermediate.6M.treatment.hypermethylation.genelist <- get.genelist(intermediate.6M.treatment.hypermethylation)
dim(intermediate.6M.treatment.hypermethylation.genelist)
genelist.name="treatment.specific.hypomethylation_intermediate.6M.newdata"
intermediate.6M.treatment.hypomethylation <- tss[which(abs(tss$c6-tss$c1)>=.05 & abs(tss$c6-tss$c1)<=.2 & (tss$csc6-tss$c6)<=-.2),]  #stable.6M.newdata
intermediate.6M.treatment.hypomethylation.genelist <- get.genelist(intermediate.6M.treatment.hypomethylation)
dim(intermediate.6M.treatment.hypomethylation)

#Get probes that are methylated/demethylated by age and then further methylated/demethylated by further treatment
#probes that show age specific (de)methylation, no treatment specific changes, but further (de)methylation in further treatment with no changes in age (de)methylayion
genelist.name="age.csc.further.hypermethylation_stable.6M.newdata"
stable.6M.further.hypermethylation <- tss[which((tss$c6-tss$c1)>=.2 & abs(tss$csc6-tss$csc1)<=.05 & (tss$csc10-tss$csc6)>=.2),] #stable.6M.newdata.further.meth
stable.6M.further.hypermethylation.genelist <- get.genelist(stable.6M.further.hypermethylation)
dim(stable.6M.further.hypermethylation)
genelist.name="age.csc.further.hypomethylation_stable.6M.newdata"
stable.6M.further.hypomethylation <- tss[which((tss$c6-tss$c1)<=-.2 & abs(tss$csc6-tss$csc1)<=.05 & (tss$csc10-tss$csc6)<=-.2),] #stable.6M.newdata.further.meth
stable.6M.further.hypomethylation.genelist <- get.genelist(stable.6M.further.hypomethylation)
dim(stable.6M.further.hypomethylation)
genelist.name="age.csc.further.hypermethylation_intermediate.6M.newdata"
intermediate.6M.further.hypermethylation <- tss[which((tss$c6-tss$c1)>=.2 & abs(tss$csc6-tss$csc1)>=.05 & abs(tss$csc6-tss$csc1)<=.2 & (tss$csc10-tss$csc6)>=.2),] #intermediate.6M.newdata.further.meth
intermediate.6M.further.hypermethylation.genelist <- get.genelist(intermediate.6M.further.hypermethylation)
dim(intermediate.6M.further.hypermethylation)
genelist.name="age.csc.further.hypomethylation_intermediate.6M.newdata"
intermediate.6M.further.hypomethylation <- tss[which((tss$c6-tss$c1)<=-.2 & abs(tss$csc6-tss$csc1)>=.05 & abs(tss$csc6-tss$csc1)<=.2 & (tss$csc10-tss$csc6)<=-.2),] #intermediate.6M.newdata.further.meth
intermediate.6M.further.hypomethylation.genelist <- get.genelist(intermediate.6M.further.hypomethylation)
dim(intermediate.6M.further.hypomethylation)

#####
# 10M
#####

#Get CpG island probes within +/-1500 bp from TSS that get methylated/demethylated by more than 20% (age/proliferation dependent methylation), and then the probes that further get methylated/demethyalted upon CSC exposure.
#Get CpG island probes that get methylated specifically upon further exposure 
#abs() <= 0.05 will pick probes that do not undergo much changes in methylation upon culturing/ageing.
#treatment specific (de)methylation
genelist.name="treatment.specific.hypermethylation_stable.10M.newdata"
stable.10M.treatment.hypermethylation <- tss[which(abs(tss$c10-tss$c1)<=.05 & (tss$csc10-tss$c10)>=.2),]  #stable.10M.newdata
stable.10M.treatment.hypermethylation.genelist <- get.genelist(stable.10M.treatment.hypermethylation)
dim(stable.10M.treatment.hypermethylation)
genelist.name="treatment.specific.hypomethylation_stable.10M.newdata"
stable.10M.treatment.hypomethylation <- tss[which(abs(tss$c10-tss$c1)<=.05 & (tss$csc10-tss$c10)<=-.2),]  #stable.10M.newdata 
stable.10M.treatment.hypomethylation.genelist <- get.genelist(stable.10M.treatment.hypomethylation)
dim(stable.10M.treatment.hypomethylation)

#Get probes that range in methylation between +/-0.05 to +/-0.2 during aging/proliferation but are methylated/demethylated >= 0.2 
genelist.name="treatment.specific.hypermethylation_intermediate.10M.newdata"
intermediate.10M.treatment.hypermethylation <- tss[which(abs(tss$c10-tss$c1)>=.05 & abs(tss$c10-tss$c1)<=.2 & (tss$csc10-tss$c10)>=.2),]  #stable.10M.newdata
intermediate.10M.treatment.hypermethylation.genelist <- get.genelist(intermediate.10M.treatment.hypermethylation)
dim(intermediate.10M.treatment.hypermethylation.genelist)
genelist.name="treatment.specific.hypomethylation_intermediate.10M.newdata"
intermediate.10M.treatment.hypomethylation <- tss[which(abs(tss$c10-tss$c1)>=.05 & abs(tss$c10-tss$c1)<=.2 & (tss$csc10-tss$c10)<=-.2),]  #stable.10M.newdata
intermediate.10M.treatment.hypomethylation.genelist <- get.genelist(intermediate.10M.treatment.hypomethylation)
dim(intermediate.10M.treatment.hypomethylation)

#Get probes that are methylated/demethylated by age and then further methylated/demethylated by further treatment
#probes that show age specific (de)methylation, no treatment specific changes, but further (de)methylation in further treatment with no changes in age (de)methylayion
genelist.name="age.csc.further.hypermethylation_stable.10M.newdata"
stable.10M.further.hypermethylation <- tss[which((tss$c10-tss$c1)>=.2 & abs(tss$csc10-tss$csc1)<=.05 & (tss$csc15-tss$csc10)>=.2),] #stable.10M.newdata.further.meth
stable.10M.further.hypermethylation.genelist <- get.genelist(stable.10M.further.hypermethylation)
dim(stable.10M.further.hypermethylation)
genelist.name="age.csc.further.hypomethylation_stable.10M.newdata"
stable.10M.further.hypomethylation <- tss[which((tss$c10-tss$c1)<=-.2 & abs(tss$csc10-tss$csc1)<=.05 & (tss$csc15-tss$csc10)<=-.2),] #stable.10M.newdata.further.meth
stable.10M.further.hypomethylation.genelist <- get.genelist(stable.10M.further.hypomethylation)
dim(stable.10M.further.hypomethylation)
genelist.name="age.csc.further.hypermethylation_intermediate.10M.newdata"
intermediate.10M.further.hypermethylation <- tss[which((tss$c10-tss$c1)>=.2 & abs(tss$csc10-tss$csc1)>=.05 & abs(tss$csc10-tss$csc1)<=.2 & (tss$csc15-tss$csc10)>=.2),] #intermediate.10M.newdata.further.meth
intermediate.10M.further.hypermethylation.genelist <- get.genelist(intermediate.10M.further.hypermethylation)
dim(intermediate.10M.further.hypermethylation)
genelist.name="age.csc.further.hypomethylation_intermediate.10M.newdata"
intermediate.10M.further.hypomethylation <- tss[which((tss$c10-tss$c1)<=-.2 & abs(tss$csc10-tss$csc1)>=.05 & abs(tss$csc10-tss$csc1)<=.2 & (tss$csc15-tss$csc10)<=-.2),] #intermediate.10M.newdata.further.meth
intermediate.10M.further.hypomethylation.genelist <- get.genelist(intermediate.10M.further.hypomethylation)
dim(intermediate.10M.further.hypomethylation)


#####
# 15M
#####

#Get CpG island probes within +/-1500 bp from TSS that get methylated/demethylated by more than 20% (age/proliferation dependent methylation), and then the probes that further get methylated/demethyalted upon CSC exposure.
#Get CpG island probes that get methylated specifically upon further exposure 
#abs() <= 0.05 will pick probes that do not undergo much changes in methylation upon culturing/ageing.
#treatment specific (de)methylation
genelist.name="treatment.specific.hypermethylation_stable.15M.newdata"
stable.15M.treatment.hypermethylation <- tss[which(abs(tss$c15-tss$c1)<=.05 & (tss$csc15-tss$c15)>=.2),]  #stable.15M.newdata
stable.15M.treatment.hypermethylation.genelist <- get.genelist(stable.15M.treatment.hypermethylation)
dim(stable.15M.treatment.hypermethylation)
genelist.name="treatment.specific.hypomethylation_stable.15M.newdata"
stable.15M.treatment.hypomethylation <- tss[which(abs(tss$c15-tss$c1)<=.05 & (tss$csc15-tss$c15)<=-.2),]  #stable.15M.newdata 
stable.15M.treatment.hypomethylation.genelist <- get.genelist(stable.15M.treatment.hypomethylation)
dim(stable.15M.treatment.hypomethylation)

#Get probes that range in methylation between +/-0.05 to +/-0.2 during aging/proliferation but are methylated/demethylated >= 0.2 
genelist.name="treatment.specific.hypermethylation_intermediate.15M.newdata"
intermediate.15M.treatment.hypermethylation <- tss[which(abs(tss$c15-tss$c1)>=.05 & abs(tss$c15-tss$c1)<=.2 & (tss$csc15-tss$c15)>=.2),]  #stable.15M.newdata
intermediate.15M.treatment.hypermethylation.genelist <- get.genelist(intermediate.15M.treatment.hypermethylation)
dim(intermediate.15M.treatment.hypermethylation.genelist)
genelist.name="treatment.specific.hypomethylation_intermediate.15M.newdata"
intermediate.15M.treatment.hypomethylation <- tss[which(abs(tss$c15-tss$c1)>=.05 & abs(tss$c15-tss$c1)<=.2 & (tss$csc15-tss$c15)<=-.2),]  #stable.15M.newdata
intermediate.15M.treatment.hypomethylation.genelist <- get.genelist(intermediate.15M.treatment.hypomethylation)
dim(intermediate.15M.treatment.hypomethylation)

dir.create(paste0(system.dir, "/Michelle/MethylationData/methylatedGeneLists/",subdir,"/"))

get.genelist <- function(x, system.di=system.dir, probe.to.gen=probe.to.gene, subdi=subdir, typ=type, genelist.nam=genelist.name){
     if (length(rownames(x)) > 0) {
          x <- cbind.data.frame(x, "Gene"=rep(NA, times=nrow(x)))
          for(j in 1:nrow(x)) {
               probe.id <- rownames(x)[j]
               gene <- probe.to.gen$gene[probe.to.gen$probe==probe.id]
               if(length(gene) == 1) { x$Gene[j] <- as.character(gene) }
          }
          write.table(x, paste(system.di, "/Michelle/MethylationData/methylatedGeneLists/", subdi,"/", typ, "_", genelist.nam, ".csv",sep=""), quote=F, row.names=T, col.names=NA, sep="\t") #row.names set to T to get the probe ids.
     }     
     tmp <- as.matrix(sort(unique(unlist(strsplit(x$Gene, ";")))))
     write.table(tmp, paste(system.di, "/Michelle/MethylationData/methylatedGeneLists/", subdi,"/", typ, "_", genelist.nam, "_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.
     return(tmp)
}
# # change dat
# stable.dat <- stable.10M.treatment.hypermethylation
# intermediate.dat <- intermediate.10M.treatment.hypermethylation
# # create directory to which genelist files will be saved to
# 
# # generate beta table with probe and gene annotation
# if (length(rownames(stable.dat)) > 0) {
#   stable.dat <- cbind.data.frame(stable.dat, "Gene"=rep(NA, times=nrow(stable.dat)))
#   for(j in 1:nrow(stable.dat)) {
#     probe.id <- rownames(stable.dat)[j]
#     gene <- probe.to.gene$gene[probe.to.gene$probe==probe.id]
#     if(length(gene) == 1) { stable.dat$Gene[j] <- as.character(gene) }
#   }
#   write.table(stable.dat, paste(system.dir, "/Michelle/MethylationData/methylatedGeneLists/", subdir,"/", type, "_", genelist.name, ".csv",sep=""), quote=F, row.names=T) #row.names set to T to get the probe ids.
# }
# 
# if (length(rownames(intermediate.dat)) > 0) {
#   intermediate.dat <- cbind.data.frame(intermediate.dat, "Gene"=rep(NA, times=nrow(intermediate.dat)))
#   for(j in 1:nrow(intermediate.dat)) {
#     probe.id <- rownames(intermediate.dat)[j]
#     gene <- probe.to.gene$gene[probe.to.gene$probe==probe.id]
#     if(length(gene) == 1) { intermediate.dat$Gene[j] <- as.character(gene) }
#   }
#   write.table(intermediate.dat, paste(system.dir, "/Michelle/MethylationData/methylatedGeneLists/", subdir,"/", type, "_", timepoint, "_age.intermediate.dat.csv",sep=""), quote=F, row.names=T) #row.names set to T to get the probe ids.
# }

# original stable.10M genelist
ten.stable.original <- as.matrix(read.table(paste0( system.dir, "/Michelle/MethylationData/methylatedGeneLists/new/age_stable_methylated_genes_at_10months_new_annotation.txt")))
ten.stable.original <- sort(unique(unlist(strsplit(ten.stable.original, ";"))))
ten.stable.original
# original intermediate.10M genelist
ten.intermediate.original <- as.matrix(read.table(paste0( system.dir, "/Michelle/MethylationData/methylatedGeneLists/new/age_intermediate_methylated_genes_at_10months_new_annotation.txt")))
ten.intermediate.original <- sort(unique(unlist(strsplit(ten.intermediate.original, ";"))))
ten.intermediate.original

# Find common genes between old and new stable genelists - 10month only
ten.stable.new <- stable.dat$Gene
ten.stable.new <- sort(unique(unlist(strsplit(ten.stable.new, ";"))))
# save the new stable genelist
write.table(ten.stable.new, file=paste(system.dir, "/Michelle/MethylationData/methylatedGeneLists/", subdir,"/", type, "_", timepoint, "_stable.genelist_new.data.txt",sep=""), quote=F, row.names=F, col.names=FALSE) #row.names set to T to get the probe ids.
stable.csc15csc1 <- as.matrix(intersect(ten.stable.original, ten.stable.new))
colnames(stable.csc10csc1) <- "# common between old and new"
write.table(stable.csc10csc1, file=paste(system.dir, "/Michelle/MethylationData/methylatedGeneLists/", subdir,"/", type, "_", timepoint, "_common.stable.txt",sep=""), quote=F, row.names=F, col.names=FALSE) #row.names set to T to get the probe ids.
length(stable.csc10csc1)
length(ten.stable.original)
length(ten.stable.new)
# Find genes unique to old and unique to new stable genelists
unique.stable.old <- setdiff(ten.stable.original, ten.stable.new)
write.table(unique.stable.old, paste(system.dir, "/Michelle/MethylationData/methylatedGeneLists/", subdir,"/", type, "_", timepoint, "_unique.stable.old.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.
unique.stable.new <- setdiff(ten.stable.new, ten.stable.original)
write.table(unique.stable.new, paste(system.dir, "/Michelle/MethylationData/methylatedGeneLists/", subdir,"/", type, "_", timepoint, "_unique.stable.new.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.

# Find common genes between old and new intermediate genelists
ten.intermediate.new <- intermediate.dat$Gene
ten.intermediate.new <- sort(unique(unlist(strsplit(ten.intermediate.new, ";"))))
# save the new intermediate genelist
write.table(ten.intermediate.new, file=paste(system.dir, "/Michelle/MethylationData/methylatedGeneLists/", subdir,"/", type, "_", timepoint, "_intermediate.genelist_new.data.txt",sep=""), quote=F, row.names=F, col.names=FALSE) #row.names set to T to get the probe ids.
intermediate.csc10csc1 <- intersect(ten.intermediate.original, ten.intermediate.new)
write.table(intermediate.csc10csc1, paste(system.dir, "/Michelle/MethylationData/methylatedGeneLists/", subdir,"/", type, "_", timepoint, "_common.intermediate.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.
length(intermediate.csc10csc1)
length(ten.intermediate.original)
length(ten.intermediate.new)
# Find genes unique to old and unique to new intermediate genelists
unique.intermediate.old <- setdiff(ten.intermediate.original, ten.intermediate.new)
write.table(unique.intermediate.old, paste(system.dir, "/Michelle/MethylationData/methylatedGeneLists/", subdir,"/", type, "_", timepoint, "_unique.intermediate.old.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.
unique.intermediate.new <- setdiff(ten.intermediate.new, ten.intermediate.original)
write.table(unique.intermediate.new, paste(system.dir, "/Michelle/MethylationData/methylatedGeneLists/", subdir,"/", type, "_", timepoint, "_unique.intermediate.new.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.


# old.intermediate <- read.table("/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/Michelle/MethylationData/methylatedGeneLists/new/age_intermediate_methylated_genes_at_10months_new_annotation.txt")
# old.stable <- read.table("/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/Michelle/MethylationData/methylatedGeneLists/new/age_stable_methylated_genes_at_10months_new_annotation.txt")
# 
# old.intermediate <- sort(unique(unlist(strsplit(as.matrix(old.intermediate), ";"))))
# old.stable <- sort(unique(unlist(strsplit(as.matrix(old.stable), ";"))))
# intersect(old.intermediate, ten.stable.new)
# intersect(old.stable, ten.stable.new)
# 
# intersect(intersect(old.stable, ten.stable.new), intersect(old.intermediate, ten.stable.new))
# 
# setdiff(ten.stable.new, old.stable)
# length(ten.stable.new)
# setdiff(ten.stable.new, old.intermediate)
# length(ten.stable.new)
