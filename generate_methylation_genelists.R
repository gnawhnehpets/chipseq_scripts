#Load libraries
# source("http://www.bioconductor.org/biocLite.R")
# biocLite("IlluminaHumanMethylation450kmanifest")
library(limma)
library(gplots)
library(heatmap.plus)
library(IlluminaHumanMethylation450k.db)
library(minfi)
library(IlluminaHumanMethylation450kmanifest)
system.dir = "/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/"
#Load images
load(paste0( system.dir, "/users/shwang26/Michelle/ash/jhu-methylation-scripts/R_image/all.probes.hg19.RData"))
load(paste0( system.dir, "/users/shwang26/Michelle/ash/jhu-methylation-scripts/R_image/cgi.xy.filter.RData"))
load(paste0( system.dir, "/users/shwang26/Michelle/ash/jhu-methylation-scripts/R_image/probesWithin5000TSS.RData"))
load(paste0( system.dir, "/users/shwang26/Michelle/ash/jhu-methylation-scripts/R_image/probe.to.gene.RData"))
load(paste0( system.dir, "/users/shwang26/Michelle/ash/jhu-methylation-scripts/R_image/annotation.RData"))

#load functions
source(paste0( system.dir, "/users/shwang26/Michelle/Rscripts/methylation450k_functions.R"))

#Get paths to input files and folders
idat <- paste0( system.dir, "/users/shwang26/Michelle/MethylationData/idat")
sheet_genelist <- paste0( system.dir, "/users/shwang26/Michelle/ash/sheet_genelist.csv")
sheet_deltaplot <- paste0( system.dir, "/users/shwang26/Michelle/ash/sheet_deltaplot.csv")
folder <- paste0( system.dir, "/users/shwang26/Michelle/ash/new_iteration/")

#Fuction for reading and normalizing data
beta.tab <- data_methylation(idat, folder)


#Retain the probes which are +/-1500 from TSS (loaded from annotation.Rdata image (for file annotation.R))
#Subset the beta.tab table to retain CpG-island probes.
dim(beta.tab)
head(beta.tab)
beta.tab_CpGI <- beta.tab[rownames(beta.tab) %in% cgi.xy.filter, ]
beta.tab_CpGI <- na.omit(beta.tab_CpGI)
dim(beta.tab_CpGI)

########################################################
#Subset the beta.tab table to retain gene body probes.
body <- beta.tab[names(body.probes), ]
body <- na.omit(body)
dim(body)

#Subset the beta.tab table to retain shelf probes.
shelf <- beta.tab[names(shelf.probes), ]
shelf <- na.omit(shelf)
dim(shelf)

#Subset the beta.tab table to retain shore probes.
shore <- beta.tab[names(shore.probes), ]
shore <- na.omit(shore)
dim(shore)
##########################################################

#get +/- 1500 TSS probes
head(probesWithin5000TSS)
tss.1500.probes <- as.character(probesWithin5000TSS[abs(probesWithin5000TSS$Dist.to.TSS) <= 1500, "ProbeName"])
beta.tab_CpGI.PlusMinus1500_TSS <- beta.tab_CpGI[rownames(beta.tab_CpGI) %in% tss.1500.probes,] #these are the CpGI probes within +/-1500 bp from TSS
beta.tab_CpGI.PlusMinus1500_TSS <- na.omit(beta.tab_CpGI.PlusMinus1500_TSS)
beta.tab_CpGI<-data.frame(beta.tab_CpGI)
beta.tab_CpGI.PlusMinus1500_TSS<-data.frame(beta.tab_CpGI.PlusMinus1500_TSS)
dim(beta.tab_CpGI.PlusMinus1500_TSS)

#Write TSS +/-1500 probes to file with gene names
if (length(rownames(beta.tab_CpGI.PlusMinus1500_TSS)) > 0) {
  beta.tab_CpGI.PlusMinus1500_TSS <- cbind.data.frame(beta.tab_CpGI.PlusMinus1500_TSS, "Gene"=rep(NA, times=nrow(beta.tab_CpGI.PlusMinus1500_TSS)))
  for(j in 1:nrow(beta.tab_CpGI.PlusMinus1500_TSS)) {
    probe.id <- rownames(beta.tab_CpGI.PlusMinus1500_TSS)[j]
    gene <- probe.to.gene$gene[probe.to.gene$probe==probe.id]
    if(length(gene) == 1) { beta.tab_CpGI.PlusMinus1500_TSS$Gene[j] <- as.character(gene) }
  }
  write.csv(beta.tab_CpGI.PlusMinus1500_TSS,paste(folder, "tss_probes.csv",sep="") , quote=F, row.names=T) #row.names set to T to get the probe ids.
}


# #Analyze methylation data and create genelist
# 
# analysis_methylation(beta.tab_CpGI.PlusMinus1500_TSS, sheet_genelist, folder, "TSS") 
# analysis_methylation(body, sheet_genelist, folder, "Body") 
# analysis_methylation(shore, sheet_genelist, folder, "Shore") 
# analysis_methylation(shelf, sheet_genelist, folder, "Shelf") 
# 
# #Make delta plots
# delta_plot(beta.tab_CpGI.PlusMinus1500_TSS, sheet_deltaplot, folder, "TSS") 
# delta_plot(body, sheet_deltaplot, folder, "Body") 
# delta_plot(shelf, sheet_deltaplot, folder, "Shelf") 
# delta_plot(shore, sheet_deltaplot, folder, "Shore") 

#save image
save(beta.tab_CpGI,beta.tab_CpGI.PlusMinus1500_TSS,file="michelle_new-annotation.RData")

# transition
Dir <- folder
subDir_file <- "files"
type <- "TSS"

head(beta.tab)
head(beta.tab_old)
dim(beta.tab)
dim(beta.tab_old)
head(beta.tab_CpGI.PlusMinus1500_TSS)
head(beta.tab_old_tss)

tss <- beta.tab_CpGI.PlusMinus1500_TSS
tss_reorder <- cbind.data.frame(tss$C.1M, tss$CSC.1M, tss$C.6M, tss$CSC.6M, tss$C.10M, tss$CSC.10M, tss$C.15M, tss$CSC.15M)
colnames(tss_reorder) <- c("c1", "csc1", "c6", "csc6", "c10", "csc10", "c15", "csc15")
tss <- tss_reorder
colnames(tss)
rownames(tss) <- rownames(beta.tab_CpGI.PlusMinus1500_TSS)
cbind(tss$c1, tss$csc1)
#Get CpG island probes within +/-1500 bp from TSS that get methylated/demethylated by more than 20% (age/proliferation dependent methylation), and then the probes that further get methylated/demethyalted upon CSC exposure.
#10M
#treatment specific changes

# #Get CpG island probes that get methylated specifically upon further exposure 
# #abs() <= 0.05 will pick probes that do not undergo much changes in methylation upon culturing/ageing.
# #mine
# stable.hypermethylation <- tss[which(abs(tss$c10-tss$c1)<=.05 & (tss$csc10-tss$c10)>=.2),] 
# #ash
# stable.hypermethylation <- tss[which(abs(tss$c10-tss$c6)<=.05 & (tss$csc10-tss$c10)>=.2),] 
# stable.hypomethylation <- tss[which(abs(tss$c10-tss$c1)<=.05 & (tss$csc10-tss$c10)<=-.2),] 
# 
# dim(stable.hypermethylation)
# dim(stable.hypomethylation)
# 
# #Get probes that range in methylation between +/-0.05 to +/-0.2 during aging/proliferation but are methylated/demethylated >= 0.2 
# #mine
# intermediate.hypermethylation <- tss[which(abs(tss$c10-tss$c1)>=.05 & abs(tss$c10-tss$c1)<=.2 & (tss$csc10-tss$c10)>=.2),] 
# #ash
# intermediate.hypermethylation <- tss[which(abs(tss$c10-tss$c6)>=.05 & abs(tss$c10-tss$c6)<=.2 & (tss$csc10-tss$c10)>=.2),] 
# intermediate.hypomethylation <- tss[which(abs(tss$c10-tss$c1)>=.05 & abs(tss$c10-tss$c1)<=.2 & (tss$csc10-tss$c10)<=-.2),] 

#Get CpG island probes that get methylated specifically upon further exposure 
#abs() <= 0.05 will pick probes that do not undergo much changes in methylation upon culturing/ageing.
#mine
stable.hypermethylation <- tss[which(abs(tss$c10-tss$c1)<=.05 & (tss$csc10-tss$c10)>=.2),] 
#ash
# stable.hypermethylation <- tss[which(abs(tss$c10-tss$c1)<=.05 & (tss$csc10-tss$c10)>=.2),] 
stable.hypomethylation <- tss[which(abs(tss$c10-tss$c1)<=.05 & (tss$csc10-tss$c10)<=-.2),] 

dim(stable.hypermethylation)
dim(stable.hypomethylation)

#Get probes that range in methylation between +/-0.05 to +/-0.2 during aging/proliferation but are methylated/demethylated >= 0.2 
#mine
intermediate.hypermethylation <- tss[which(abs(tss$c10-tss$c1)>=.05 & abs(tss$c10-tss$c1)<=.2 & (tss$csc10-tss$c10)>=.2),] 
#ash
# intermediate.hypermethylation <- tss[which(abs(tss$c10-tss$c6)>=.05 & abs(tss$c10-tss$c6)<=.2 & (tss$csc10-tss$c10)>=.2),] 
intermediate.hypomethylation <- tss[which(abs(tss$c10-tss$c1)>=.05 & abs(tss$c10-tss$c1)<=.2 & (tss$csc10-tss$c10)<=-.2),] 


dim(intermediate.hypermethylation)
dim(intermediate.hypomethylation)

ten.stable.original <- as.matrix(read.table(paste0( system.dir, "/users/shwang26/Michelle/MethylationData/methylatedGeneLists/new/age_stable_methylated_genes_at_10months_new_annotation.txt")))
ten.stable.original <- sort(unique(unlist(strsplit(ten.stable.original, ";"))))
ten.stable.original

ten.intermediate.original <- as.matrix(read.table(paste0( system.dir, "/users/shwang26/Michelle/MethylationData/methylatedGeneLists/new/age_intermediate_methylated_genes_at_10months_new_annotation.txt")))
ten.intermediate.original <- sort(unique(unlist(strsplit(ten.intermediate.original, ";"))))
ten.intermediate.original

stable.hypermethylation
head(probe.to.gene)
if (length(rownames(stable.hypermethylation)) > 0) {
  stable.hypermethylation <- cbind.data.frame(stable.hypermethylation, "Gene"=rep(NA, times=nrow(stable.hypermethylation)))
  for(j in 1:nrow(stable.hypermethylation)) {
    probe.id <- rownames(stable.hypermethylation)[j]
    gene <- probe.to.gene$gene[probe.to.gene$probe==probe.id]
    if(length(gene) == 1) { stable.hypermethylation$Gene[j] <- as.character(gene) }
  }
  write.csv(stable.hypermethylation, paste(Dir, subDir_file,"/", type, "_age.stable.hypermethylation_10M.csv",sep=""), quote=F, row.names=T) #row.names set to T to get the probe ids.
}

if (length(rownames(intermediate.hypermethylation)) > 0) {
  intermediate.hypermethylation <- cbind.data.frame(intermediate.hypermethylation, "Gene"=rep(NA, times=nrow(intermediate.hypermethylation)))
  for(j in 1:nrow(intermediate.hypermethylation)) {
    probe.id <- rownames(intermediate.hypermethylation)[j]
    gene <- probe.to.gene$gene[probe.to.gene$probe==probe.id]
    if(length(gene) == 1) { intermediate.hypermethylation$Gene[j] <- as.character(gene) }
  }
  write.csv(intermediate.hypermethylation, paste(Dir, subDir_file,"/", type, "_age.intermediate.hypermethylation_10M.csv",sep=""), quote=F, row.names=T) #row.names set to T to get the probe ids.
}


ten.stable.new <- stable.hypermethylation$Gene
ten.stable.new <- sort(unique(unlist(strsplit(ten.stable.new, ";"))))
stable.csc10csc6 <- intersect(ten.stable.original, ten.stable.new)
write.csv(stable.csc10csc6, paste(Dir, subDir_file,"/", type, "_common.stable.csv",sep=""), quote=F, row.names=T) #row.names set to T to get the probe ids.
length(ten.stable.original)
length(ten_new)
unique.stable.old <- setdiff(ten.stable.original, ten.stable.new)
write.csv(stable.csc10csc6, paste(Dir, subDir_file,"/", type, "_unique.stable.old.csv",sep=""), quote=F, row.names=T) #row.names set to T to get the probe ids.
unique.stable.new <- setdiff(ten.stable.new, ten.stable.original)
write.csv(stable.csc10csc6, paste(Dir, subDir_file,"/", type, "_unique.stable.new.csv",sep=""), quote=F, row.names=T) #row.names set to T to get the probe ids.


ten.intermediate.new <- intermediate.hypermethylation$Gene
ten.intermediate.new <- sort(unique(unlist(strsplit(ten.intermediate.new, ";"))))
intermediate.csc10csc6 <- intersect(ten.intermediate.original, ten.intermediate.new)
write.csv(stable.csc10csc6, paste(Dir, subDir_file,"/", type, "_common.intermediate.csv",sep=""), quote=F, row.names=T) #row.names set to T to get the probe ids.
length(ten.intermediate.original)
length(ten_new)
unique.intermediate.old <- setdiff(ten.intermediate.original, ten.intermediate.new)
write.csv(stable.csc10csc6, paste(Dir, subDir_file,"/", type, "_unique.intermediate.old.csv",sep=""), quote=F, row.names=T) #row.names set to T to get the probe ids.
unique.intermediate.new <- setdiff(ten.intermediate.new, ten.intermediate.original)
write.csv(stable.csc10csc6, paste(Dir, subDir_file,"/", type, "_unique.intermediate.new.csv",sep=""), quote=F, row.names=T) #row.names set to T to get the probe ids.