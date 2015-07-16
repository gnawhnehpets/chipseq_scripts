##########  
# Purpose:
##########
print("CORRECTED")
library(limma)
library(gplots)
library(heatmap.plus)
library(IlluminaHumanMethylation450k.db)
library(minfi)
library(IlluminaHumanMethylation450kmanifest)
system.dir = "/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/"
# system.dir="Z:/users/shwang26/"
#Load images created by Ash
a<-load(paste0( system.dir, "/Michelle/ash/jhu-methylation-scripts/R_image/all.probes.hg19.RData"))
b<-load(paste0( system.dir, "/Michelle/ash/jhu-methylation-scripts/R_image/cgi.xy.filter.RData"))
c<-load(paste0( system.dir, "/Michelle/ash/jhu-methylation-scripts/R_image/probesWithin5000TSS.RData"))
d<-load(paste0( system.dir, "/Michelle/ash/jhu-methylation-scripts/R_image/probe.to.gene.RData"))
e<-load(paste0( system.dir, "/Michelle/ash/jhu-methylation-scripts/R_image/annotation.RData"))
e
#Load object, new beta tables generated with new data
load(file=paste0(system.dir, "/Michelle/Robjects/michelle_new-annotation_old.RData")) #REP1

class(beta.tab_CpGI) #all cpgi
class(beta.tab_CpGI.PlusMinus1500_TSS) #cpgi promoters

#####################################
rep1_mds_beta.tab_CpGI <- as.matrix(beta.tab_CpGI)
head(rep1_mds_beta.tab_CpGI)
colnames(rep1_mds_beta.tab_CpGI) <- paste0("rep1_", colnames(rep1_mds_beta.tab_CpGI))
rep1_mds_beta.tab_CpGI <- as.matrix(rep1_mds_beta.tab_CpGI)
head(rep1_mds_beta.tab_CpGI)
#####################################
rep1_mds_beta.tab_CpGI.PlusMinus1500_TSS <- as.matrix(beta.tab_CpGI.PlusMinus1500_TSS[,-9])
head(rep1_mds_beta.tab_CpGI.PlusMinus1500_TSS)
colnames(rep1_mds_beta.tab_CpGI.PlusMinus1500_TSS) <- paste0("rep1_", colnames(rep1_mds_beta.tab_CpGI.PlusMinus1500_TSS))
rep1_mds_beta.tab_CpGI.PlusMinus1500_TSS <- as.matrix(rep1_mds_beta.tab_CpGI.PlusMinus1500_TSS)
head(rep1_mds_beta.tab_CpGI.PlusMinus1500_TSS)

#####################################

mdsPlot(rep1_mds_beta.tab_CpGI, pch=2, pal=c("red3", "green4"), sampGroups=factor(gsub("\\.\\d+", "", colnames(rep1_mds_beta.tab_CpGI))), legendPos="topleft", legendNCol=1, sampNames=colnames(rep1_mds_beta.tab_CpGI))
mdsPlot(rep1_mds_beta.tab_CpGI, pch=2, pal=c("red3", "green4"), sampGroups=factor(gsub("\\.\\d+", "", colnames(rep1_mds_beta.tab_CpGI))), legendPos="topleft", legendNCol=1)

mdsPlot(rep1_mds_beta.tab_CpGI.PlusMinus1500_TSS, pch=2, pal=c("red3", "green4"), sampGroups=factor(gsub("\\.\\d+", "", colnames(rep1_mds_beta.tab_CpGI))), legendPos="topleft", legendNCol=1, sampNames=colnames(rep1_mds_beta.tab_CpGI))
mdsPlot(rep1_mds_beta.tab_CpGI.PlusMinus1500_TSS, pch=2, pal=c("red3", "green4"), sampGroups=factor(gsub("\\.\\d+", "", colnames(rep1_mds_beta.tab_CpGI))), legendPos="topleft", legendNCol=1)


load(file=paste0(system.dir, "/Michelle/Robjects/michelle_new-annotation.RData")) #REP2

#####################################
rep2_mds_beta.tab_CpGI <- as.matrix(beta.tab_CpGI)
head(rep2_mds_beta.tab_CpGI)
colnames(rep2_mds_beta.tab_CpGI) <- paste0("rep2_", gsub("M", "", colnames(rep2_mds_beta.tab_CpGI)))
colnames(rep2_mds_beta.tab_CpGI) <- gsub("\\_C\\.", "\\_UNT\\.", colnames(rep2_mds_beta.tab_CpGI))

rep2_mds_beta.tab_CpGI <- as.matrix(rep2_mds_beta.tab_CpGI)
head(rep2_mds_beta.tab_CpGI)
#####################################
rep2_mds_beta.tab_CpGI.PlusMinus1500_TSS <- as.matrix(beta.tab_CpGI.PlusMinus1500_TSS[,-9])
head(rep2_mds_beta.tab_CpGI.PlusMinus1500_TSS)
colnames(rep2_mds_beta.tab_CpGI.PlusMinus1500_TSS) <- paste0("rep2_", gsub("M", "", colnames(rep2_mds_beta.tab_CpGI.PlusMinus1500_TSS)))
colnames(rep2_mds_beta.tab_CpGI.PlusMinus1500_TSS) <- gsub("\\_C\\.", "\\_UNT\\.", colnames(rep2_mds_beta.tab_CpGI.PlusMinus1500_TSS))
rep2_mds_beta.tab_CpGI.PlusMinus1500_TSS <- as.matrix(rep2_mds_beta.tab_CpGI.PlusMinus1500_TSS)
head(rep2_mds_beta.tab_CpGI.PlusMinus1500_TSS)

#####################################

mdsPlot(rep2_mds_beta.tab_CpGI, pch=2, pal=c("red3", "green4"), sampGroups=factor(gsub("\\.\\d+", "", colnames(rep2_mds_beta.tab_CpGI))), legendPos="topleft", legendNCol=1, sampNames=colnames(rep2_mds_beta.tab_CpGI))
mdsPlot(rep2_mds_beta.tab_CpGI, pch=2, pal=c("red3", "green4"), sampGroups=factor(gsub("\\.\\d+", "", colnames(rep2_mds_beta.tab_CpGI))), legendPos="topleft", legendNCol=1)

mdsPlot(rep2_mds_beta.tab_CpGI.PlusMinus1500_TSS, pch=2, pal=c("red3", "green4"), sampGroups=factor(gsub("\\.\\d+", "", colnames(rep2_mds_beta.tab_CpGI))), legendPos="topleft", legendNCol=1, sampNames=colnames(rep2_mds_beta.tab_CpGI))
mdsPlot(rep2_mds_beta.tab_CpGI.PlusMinus1500_TSS, pch=2, pal=c("red3", "green4"), sampGroups=factor(gsub("\\.\\d+", "", colnames(rep2_mds_beta.tab_CpGI))), legendPos="topleft", legendNCol=1)


#####################################
#####################################
#####################################
#ALL PROBES
length(rownames(rep1_mds_beta.tab_CpGI))
length(rownames(rep2_mds_beta.tab_CpGI))
length(intersect(rownames(rep1_mds_beta.tab_CpGI), rownames(rep2_mds_beta.tab_CpGI)))
common.probes <- intersect(rownames(rep1_mds_beta.tab_CpGI), rownames(rep2_mds_beta.tab_CpGI))
rep1_mds_beta <- rep1_mds_beta.tab_CpGI[common.probes, ]
rep2_mds_beta <- rep2_mds_beta.tab_CpGI[common.probes, ]
head(rep1_mds_beta)
head(rep2_mds_beta)
dim(rep1_mds_beta.tab_CpGI)
dim(rep2_mds_beta.tab_CpGI)
class(rep1_mds_beta.tab_CpGI)
class(rep2_mds_beta.tab_CpGI)

mds_combined_CpGI <- cbind(rep1_mds_beta, rep2_mds_beta)

head(mds_combined_CpGI)
jpeg(filename=paste0(system.dir, "Michelle/Analysis/pca/new/pca_rep1&rep2_sample.names_allprobes.jpeg"), height=600, width=900)
mdsPlot(mds_combined_CpGI, pch=1, pal=c("red", "gray4", "orange", "gray40"), sampGroups=factor(gsub("\\.\\d+", "", colnames(mds_combined_CpGI))), legendPos="topleft", legendNCol=1, sampNames=colnames(mds_combined_CpGI))
dev.off()
jpeg(filename=paste0(system.dir, "Michelle/Analysis/pca/new/pca_rep1&rep2_points_allprobes.jpeg"), height=600, width=900)
mdsPlot(mds_combined_CpGI, pch=17, pal=c("red", "gray4", "orange", "gray40"), sampGroups=factor(gsub("\\.\\d+", "", colnames(mds_combined_CpGI))), legendPos="topleft", legendNCol=1)
dev.off()

# BODY PROBES
length(body.probes)
names(body.probes)
body <- intersect(common.probes, names(body.probes))
rep1_mds_beta <- rep1_mds_beta.tab_CpGI[body, ]
rep2_mds_beta <- rep2_mds_beta.tab_CpGI[body, ]
head(rep1_mds_beta)
head(rep2_mds_beta)
dim(rep1_mds_beta)
dim(rep2_mds_beta)
class(rep1_mds_beta.tab_CpGI)
class(rep2_mds_beta.tab_CpGI)

mds_combined_CpGI1 <- cbind(rep1_mds_beta, rep2_mds_beta)


jpeg(filename=paste0(system.dir, "Michelle/Analysis/pca/new/pca_rep1&rep2_sample.names_bodyprobes.jpeg"), height=600, width=900)
mdsPlot(mds_combined_CpGI1, pch=1, pal=c("red", "gray4", "orange", "gray40"), sampGroups=factor(gsub("\\.\\d+", "", colnames(mds_combined_CpGI))), legendPos="topleft", legendNCol=1, sampNames=colnames(mds_combined_CpGI))
dev.off()
jpeg(filename=paste0(system.dir, "Michelle/Analysis/pca/new/pca_rep1&rep2_points_bodyprobes.jpeg"), height=600, width=900)
mdsPlot(mds_combined_CpGI1, pch=17, pal=c("red", "gray4", "orange", "gray40"), sampGroups=factor(gsub("\\.\\d+", "", colnames(mds_combined_CpGI))), legendPos="topleft", legendNCol=1)
dev.off()

# SHORE PROBES
length(shore.probes)
names(shore.probes)
shore <- intersect(common.probes, names(shore.probes))
rep1_mds_beta <- rep1_mds_beta.tab_CpGI[shore, ]
rep2_mds_beta <- rep2_mds_beta.tab_CpGI[shore, ]
head(rep1_mds_beta)
head(rep2_mds_beta)
dim(rep1_mds_beta)
dim(rep2_mds_beta)
class(rep1_mds_beta.tab_CpGI)
class(rep2_mds_beta.tab_CpGI)

mds_combined_CpGI2 <- cbind(rep1_mds_beta, rep2_mds_beta)


jpeg(filename=paste0(system.dir, "Michelle/Analysis/pca/new/pca_rep1&rep2_sample.names_shoreprobes.jpeg"), height=600, width=900)
mdsPlot(mds_combined_CpGI2, pch=1, pal=c("red", "gray4", "orange", "gray40"), sampGroups=factor(gsub("\\.\\d+", "", colnames(mds_combined_CpGI))), legendPos="topleft", legendNCol=1, sampNames=colnames(mds_combined_CpGI))
dev.off()
jpeg(filename=paste0(system.dir, "Michelle/Analysis/pca/new/pca_rep1&rep2_points_shoreprobes.jpeg"), height=600, width=900)
mdsPlot(mds_combined_CpGI2, pch=17, pal=c("red", "gray4", "orange", "gray40"), sampGroups=factor(gsub("\\.\\d+", "", colnames(mds_combined_CpGI))), legendPos="topleft", legendNCol=1)
dev.off()
e
# SHELF PROBES
length(shore.probes)
names(shore.probes)
shelf <- intersect(common.probes, names(shelf.probes))
rep1_mds_beta <- rep1_mds_beta.tab_CpGI[shelf, ]
rep2_mds_beta <- rep2_mds_beta.tab_CpGI[shelf, ]
head(rep1_mds_beta)
head(rep2_mds_beta)
dim(rep1_mds_beta)
dim(rep2_mds_beta)
class(rep1_mds_beta.tab_CpGI)
class(rep2_mds_beta.tab_CpGI)

mds_combined_CpGI2 <- cbind(rep1_mds_beta, rep2_mds_beta)


jpeg(filename=paste0(system.dir, "Michelle/Analysis/pca/new/pca_rep1&rep2_sample.names_shelfprobes.jpeg"), height=600, width=900)
mdsPlot(mds_combined_CpGI2, pch=1, pal=c("red", "gray4", "orange", "gray40"), sampGroups=factor(gsub("\\.\\d+", "", colnames(mds_combined_CpGI))), legendPos="topleft", legendNCol=1, sampNames=colnames(mds_combined_CpGI))
dev.off()
jpeg(filename=paste0(system.dir, "Michelle/Analysis/pca/new/pca_rep1&rep2_points_shelfprobes.jpeg"), height=600, width=900)
mdsPlot(mds_combined_CpGI2, pch=17, pal=c("red", "gray4", "orange", "gray40"), sampGroups=factor(gsub("\\.\\d+", "", colnames(mds_combined_CpGI))), legendPos="topleft", legendNCol=1)
dev.off()

e
# TSS PROBES
length(tss.probes)
names(tss.probes)
tss <- intersect(common.probes, names(tss.probes))
rep1_mds_beta <- rep1_mds_beta.tab_CpGI[tss, ]
rep2_mds_beta <- rep2_mds_beta.tab_CpGI[tss, ]
head(rep1_mds_beta)
head(rep2_mds_beta)
dim(rep1_mds_beta)
dim(rep2_mds_beta)
class(rep1_mds_beta.tab_CpGI)
class(rep2_mds_beta.tab_CpGI)

mds_combined_CpGI2 <- cbind(rep1_mds_beta, rep2_mds_beta)


jpeg(filename=paste0(system.dir, "Michelle/Analysis/pca/new/pca_rep1&rep2_sample.names_tssprobes.jpeg"), height=600, width=900)
mdsPlot(mds_combined_CpGI2, pch=1, pal=c("red", "gray4", "orange", "gray40"), sampGroups=factor(gsub("\\.\\d+", "", colnames(mds_combined_CpGI))), legendPos="topleft", legendNCol=1, sampNames=colnames(mds_combined_CpGI))
dev.off()
jpeg(filename=paste0(system.dir, "Michelle/Analysis/pca/new/pca_rep1&rep2_points_tssprobes.jpeg"), height=600, width=900)
mdsPlot(mds_combined_CpGI2, pch=17, pal=c("red", "gray4", "orange", "gray40"), sampGroups=factor(gsub("\\.\\d+", "", colnames(mds_combined_CpGI))), legendPos="topleft", legendNCol=1)
dev.off()

#TSS 1500
length(rownames(rep1_mds_beta.tab_CpGI))
length(rownames(rep2_mds_beta.tab_CpGI))
length(intersect(rownames(rep1_mds_beta.tab_CpGI), rownames(rep2_mds_beta.tab_CpGI)))
common.tss1500.probes <- intersect(rownames(rep1_mds_beta.tab_CpGI.PlusMinus1500_TSS), rownames(rep2_mds_beta.tab_CpGI.PlusMinus1500_TSS))
rep1_mds_beta <- rep1_mds_beta.tab_CpGI.PlusMinus1500_TSS[common.tss1500.probes, ]
rep2_mds_beta <- rep2_mds_beta.tab_CpGI.PlusMinus1500_TSS[common.tss1500.probes, ]
head(rep1_mds_beta)
head(rep2_mds_beta)
dim(rep1_mds_beta.tab_CpGI)
dim(rep2_mds_beta.tab_CpGI)
class(rep1_mds_beta.tab_CpGI)
class(rep2_mds_beta.tab_CpGI)

mds_combined_CpGI <- cbind(rep1_mds_beta, rep2_mds_beta)

head(mds_combined_CpGI)
jpeg(filename=paste0(system.dir, "Michelle/Analysis/pca/new/pca_rep1&rep2_sample.names_tss1500probes.jpeg"), height=600, width=900)
mdsPlot(mds_combined_CpGI, pch=1, pal=c("red", "gray4", "orange", "gray40"), sampGroups=factor(gsub("\\.\\d+", "", colnames(mds_combined_CpGI))), legendPos="topleft", legendNCol=1, sampNames=colnames(mds_combined_CpGI))
dev.off()
jpeg(filename=paste0(system.dir, "Michelle/Analysis/pca/new/pca_rep1&rep2_points_tss1500probes.jpeg"), height=600, width=900)
mdsPlot(mds_combined_CpGI, pch=17, pal=c("red", "gray4", "orange", "gray40"), sampGroups=factor(gsub("\\.\\d+", "", colnames(mds_combined_CpGI))), legendPos="topleft", legendNCol=1)
dev.off()