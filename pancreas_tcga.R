# Analysis of TCGA pancreatic data

########################  
# Load libraries
########################
#rm(list=ls()); ls()
library(limma)  
library(gplots)
library(heatmap.plus)
dir = "/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/tcga/"
dir = "/amber2/scratch/baylin/shwang/tcga/"
setwd(dir)

# Illumina 450K probe annotation
library(IlluminaHumanMethylation450k.db)
probe.loc.in.gene <- as.list(IlluminaHumanMethylation450kPROBELOCATION) #Map a probe ID with its location in the gene
CpGloc <- as.list(IlluminaHumanMethylation450kCPGILOCATION) #IlluminaHumanMethylation450kCPGILOCATION maps between Illumina probe IDs and the UCSC CpG island features, if any, with which they are associated
island.probes <- unlist(sapply(c(1:length(CpGloc)), FUN=function(i) {return(grep("island" , CpGloc[i], ignore.case=T, value=T))})) #gets island probes
tss.probes <- unlist(sapply(c(1:length(probe.loc.in.gene)), FUN=function(i) {return(grep("TSS" , probe.loc.in.gene[i], ignore.case=T, value=T))})) #gets probes that are near TSS (200 or 1500 bp near TSS)
chromosome <- as.list(IlluminaHumanMethylation450kCHR)
x.probes <- unlist(sapply(c(1:length(probe.loc.in.gene)), FUN=function(i) {return(grep("X" , chromosome[i], ignore.case=T, value=T))})) #gets X chromosome probes
y.probes <- unlist(sapply(c(1:length(probe.loc.in.gene)), FUN=function(i) {return(grep("Y" , chromosome[i], ignore.case=T, value=T))})) #gets Y chromosome probes
island.tss.probes <- intersect(names(island.probes), names(tss.probes))
length(which(chromosome %in% y.probes))
length(which(chromosome %in% x.probes))
# remove sex chrom
autosome.probes <- names(chromosome[-c(which(chromosome %in% y.probes), which(chromosome %in% x.probes))])
# load("./Robjects/probesWithin5000TSS.RData")
# head(probesWithin5000TSS)
# save(probesWithin5000TSS, autosome.probes, chromosome, probe.loc.in.gene, CpGloc, island.probes, tss.probes, island.tss.probes, x.probes, y.probes, file="./Robjects/annotationobjects.rda")
load(file="./Robjects/annotationobjects.rda")
load(file="./Robjects/merged.all.rda")
# indices of tumors/normals
merged.tumors <- grep("tumors",colnames(merged.all[,-c(1:2)]))
merged.normals <- grep("normals",colnames(merged.all[,-c(1:2)]))
merged.tumors.names <- colnames(merged.all)[merged.tumors]
merged.normals.names <- colnames(merged.all)[merged.normals]
length(merged.tumors.names)
length(merged.normals.names)
colnames(merged.all)
colnames(merged.all)[grep("breast", colnames(merged.all))]
colnames(merged.all)[grep("colon", colnames(merged.all))]
colnames(merged.all)[grep("lung", colnames(merged.all))]
colnames(merged.all)[grep("pancreas", colnames(merged.all))]



########################
#Functions
########################
fun.loadTCGA.Level3_Data <- function(files){#This function reads the level3 Infinium data, orders the probes in the same order across teh samples, and outputs a table of beta values. files is the vector of file names of level3 data downloaded from TCGA.
  dat <- lapply(files, function(filename){
    a <- scan(filename, what=list("Sample"="", "Probes"="", "beta"=0, "Gene"="", NULL, NULL), sep="\t", skip=1)
    a <- a[-1*c(5,6)]
    #    names(a) <- c("Probes", "M.signal", "U.signal", "Detection_P_Value")
    return(a)
  })
  sample.names <- sapply(c(1:length(dat)), FUN=function(i){unique(dat[[i]]$Sample)})
  common.probes <- lapply(c(1:length(dat)), FUN=function(i){dat[[i]]$Probes})
  common.probes <- Reduce(intersect, common.probes)
  #Create a new list dat1 containing beta values from all samples aranged in the same order and corresponding to the order of common.probes
  dat.beta <- lapply(c(1:length(dat)), FUN=function(i){
    x <- match(common.probes, dat[[i]]$Probes) #gets the positions in the 2nd argument that matches those in the 1st argument and in te order of teh 1st argument- hence it gets probes for each sample in the order of the common.probes
    beta.x <- dat[[i]]$beta[x]
    #Probes.x <- dat[[i]]$Probes[x]
    #Gene.x <- dat[[i]]$Gene[x]
    #xx <- list("Probes"=Probes.x, "Gene"=Gene.x, "beta"=beta.x)
    return(beta.x)
  })
  names(dat.beta) <- sample.names
  beta.tab <- do.call(cbind, dat.beta)
  #Get the Probes and Genes in the right order corresponding to common.probes from one of the samples (doesn't matter which one because the beta values are aligned to common.probes above)
  x <- match(common.probes, dat[[1]]$Probes) #as above
  Probes.x <- dat[[1]]$Probes[x]
  Gene.x <- dat[[1]]$Gene[x]
  beta.tab <- cbind.data.frame("Probes"=Probes.x, "Gene"=Gene.x, beta.tab)
  return(beta.tab)
}
# colon.tcga.files <- paste("/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/tcga/Pancreatic_grant/TCGA_data/Colon/DNA_Methylation/JHU_USC__HumanMethylation450/Level_3/", dir("/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/tcga/Pancreatic_grant/TCGA_data/Colon/DNA_Methylation/JHU_USC__HumanMethylation450/Level_3"), sep="")
# colon.tcga.files[1]
# Colon.Data_Meth <- funlvl3(colon.tcga.files[1])
# files <- colon.tcga.files[1]
# dir.to.files <- colon.formatted.dir
################################################################################################################################
################################################################################################################################
################################################################################################################################
fun_convert_tcga_to_newformat <- function(files, dir.to.files){
  path.dir <- dir.to.files
  print(path.dir)
  for(i in 1:length(files)){
    print(paste("File #",i," formatting start...", sep=""))
    barcode <- lapply(files[i], function(filename){
      header <- scan(filename, what=list("H"="", "barcode"="", NULL, NULL, NULL), nlines=1, sep="\t")
      header <- header[[2]]
      return(header)
    })
  #   barcode
    #import_dat <- lapply(paste("/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/tcga/tcgatest/",dir("/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/tcga/tcgatest/"), sep=""), read.delim)
    # read-in tcga data
    import_dat <- lapply(files[i], read.delim)
  #   length(import_dat)
    # get it into __methylation_analysis.txt format
    tmpdat <- import_dat[[1]][-1,]  
    tmpdat <- cbind(barcode, import_dat[[1]][-1,])
    colnames(tmpdat) <- c("barcode", "probe name", "beta value", "gene symbol", "chromosome", "position")
  #   head(tmpdat, n=5)
    # designate new filename
    tmpfile=paste("jhu-usc.edu__HumanMethylation450__", barcode, "__methylation_analysis.txt",sep="")
    filename <- paste(path.dir, tmpfile, sep="")
    print(filename)
  #   filename
    write.table(tmpdat, file=filename, append=FALSE, quote=FALSE, sep="\t", na="NA", row.names=FALSE, col.names=TRUE)
  #   filename
  #   getwd()
    print(paste("File #",i," complete.", sep=""))
  }
}

pancreas.tcga.files <- paste(paste(dir, "paad/DNA_Methylation/JHU_USC__HumanMethylation450/Level_3/", sep=""), dir(paste(dir, "paad/DNA_Methylation/JHU_USC__HumanMethylation450/Level_3", sep="")), sep="")
pancreas.formatted.dir <- paste(dir, "/paad/DNA_Methylation/JHU_USC__HumanMethylation450/Level_3formatted/",sep="")
fun_convert_tcga_to_newformat(pancreas.tcga.files, pancreas.formatted.dir)
pancreas.tcga.files.formatted <- paste(paste(dir, "/paad/DNA_Methylation/JHU_USC__HumanMethylation450/Level_3formatted/", sep=""), dir(paste(dir, "/paad/DNA_Methylation/JHU_USC__HumanMethylation450/Level_3formatted/", sep="")), sep="")
pancreas.data.meth <- fun.loadTCGA.Level3_Data(pancreas.tcga.files.formatted)
pancreas.data.meth <- fun.loadTCGA.Level3_Data(pancreas.tcga.files.formatted[100:195])
save(pancreas.data.meth, file="./Robjects/pancreas.data.meth") #filtered
save(pancreas.data.meth, file="./Robjects/pancreas.data.meth2") #all 197 samples
pancreas.data.meth <- as.matrix(pancreas.data.meth)
rm(pancreas.data.meth)
gc(); gc(); gc()

panc.tumors <- grep("01A", colnames(pancreas.data.meth), value=T)
panc.normals <- grep("11A", colnames(pancreas.data.meth), value=T)

#Annotate the normal and tumor samples
gc();gc();gc();gc();
# load(file="./Robjects/breast.data.meth")
# load(file="./Robjects/colon.data.meth")
# load(file="./Robjects/leukemia.data.meth")
# load(file="./Robjects/lung.data.meth")
# load(file="./Robjects/pancreas.data.meth")
load(file="./Robjects/pancreas.data.meth2")
##########################################################################################################

********************************************************************************
#PCA
********************************************************************************
x <- pancreas.data.meth[3:ncol(pancreas.data.meth)]
PCA <- prcomp(t(x))
dim(PCA$x)

# save(PCA, file="./Robjects/pancreas.PCA.rda")
load(file="./Robjects/pancreas.PCA.rda")
# save(PCA, file="./Robjects/merged.all.PCA.rda")
load(file="./Robjects/merged.all.PCA.rda")
# save(PCA, file="./Robjects/merged.filtered.PCA.rda")
load(file="./Robjects/merged.filtered.PCA.rda")
gc();gc();gc();
# save(merged.tumors, file="./Robjects/merged.tumors.rda")
# save(merged.normals, file="./Robjects/merged.normals.rda")
load(file="./Robjects/merged.tumors.rda")
load(file="./Robjects/merged.normals.rda")

pdf("./Robjects/PCA_merged.all.pdf",width=900,height=600)
title = "PCA of all cancers (filtered) colored by tumor/normal"
pdf("./Robjects/PCA_merged.filtered.pdf",width=900,height=600)
title = "PCA of all cancers (filtered) colored by tumor/normal"
plot(PCA$x[, c("PC1", "PC2")], main=title , xlab=paste("PC1: ",signif((((PCA$sdev)^2)[1]/(sum((PCA$sdev)^2))),2),sep=""), ylab=paste("PC2: ",signif((((PCA$sdev)^2)[2]/(sum((PCA$sdev)^2))),2),sep=""), pch=19, cex=0.5)
# points(matrix(PCA$x[panc.tumors, c("PC1", "PC2")], ncol=2), pch=19, col=rgb(1, 0, 0, max = 1), cex=0.5)
# points(matrix(PCA$x[panc.normals, c("PC1", "PC2")], ncol=2), pch=19, col=rgb(0, 1, 0, max = 1), cex=0.5)
points(matrix(PCA$x[merged.tumors, c("PC1", "PC2")], ncol=2), pch=19, col=rgb(1, 0, 0, max = 1), cex=0.5)
points(matrix(PCA$x[merged.normals, c("PC1", "PC2")], ncol=2), pch=19, col=rgb(0, 1, 0, max = 1), cex=0.5)
legend("topleft", legend=c("tumor", "normal"), col=c("red", "green"), cex=1, pch=19)
dev.off()



# ORIG
********************************************************************************
#Heatmaps
********************************************************************************
  #Heatmaps1
x <- pancreas.data.meth[3:ncol(pancreas.data.meth)]
rownames(x) <- pancreas.data.meth[,1]
rownames(x)
dim(x)
x <- na.omit(x)
x.tumor <- x[,panc.tumors]
dim(x); dim(x.tumor)
x.sd <- apply(x, 1, sd)
x.tumor.sd <- apply(x.tumor, 1, sd)
# plot(density(x.sd))
# plot(density(x.tumor.sd))
dim(x); dim(x.tumor)
x.all <- x[which(x.sd > 0.2), ]
x.all.tumVar <- x[which(x.tumor.sd > 0.2), ]
dim(x.all); dim(x.all.tumVar)
colsidecol <- vector()
for (i in 1:ncol(x)) {colsidecol[i] <- ifelse(colnames(x)[i] %in% panc.tumors, "red", "green")}
pdf(file="./Robjects/x.all.pdf")
heatmap.2(as.matrix(x.all), col=greenred(50), trace="none", density.info="none", hclustfun = function(x)hclust(x, method = "ward"), margins=c(10,10), cexRow=0.8, cexCol=0.7, key=T, ColSideColors=colsidecol)
dev.off()
pdf(file="./Robjects/x.all.tumVar.pdf")
heatmap.2(as.matrix(x.all.tumVar), col=greenred(50), trace="none", density.info="none", hclustfun = function(x)hclust(x, method = "ward"), margins=c(10,10), cexRow=0.8, cexCol=0.7, key=T, ColSideColors=colsidecol)
dev.off()

#Heatmaps2 - filtered
index1 <- grep("IGF2", pancreas.data.meth[,2]); index2 <- grep("SOCS1", pancreas.data.meth[,2]); index3 <- grep("RUNX3", pancreas.data.meth[,2]); index4 <- grep("CACNA1G", pancreas.data.meth[,2]); index5 <- grep("NEUROG1", pancreas.data.meth[,2])
index.all <- sort(c(index1, index2, index3, index4, index5))
# STEP 2
x <- pancreas.data.meth[index.all,]
# #remove sex chromosome
# x <- [-x.probes,]
# x <- [-y.probes,]
x <- x[3:ncol(x)]
# x <- pancreas.data.meth[3:ncol(pancreas.data.meth)]
dim(x)
x <- na.omit(x)
dim(c)
x.tumor <- x[,panc.tumors]
dim(x); dim(x.tumor)
x.sd <- apply(x, 1, sd)
x.tumor.sd <- apply(x.tumor, 1, sd)
# plot(density(x.sd))
# plot(density(x.tumor.sd))
dim(x); dim(x.tumor)
x.all <- x[which(x.sd > 0.2), ]
x.all.tumVar <- x[which(x.tumor.sd > 0.2), ]
dim(x.all); dim(x.all.tumVar)
colsidecol <- vector()
for (i in 1:ncol(x)) {colsidecol[i] <- ifelse(colnames(x)[i] %in% panc.tumors, "red", "green")}
pdf(file="./Robjects/x.filtered.pdf")
heatmap.2(as.matrix(x.all), col=greenred(50), trace="none", density.info="none", hclustfun = function(x)hclust(x, method = "ward"), margins=c(10,10), cexRow=0.8, cexCol=0.7, key=T, ColSideColors=colsidecol)
dev.off()
pdf(file="./Robjects/x.filtered.tumVar.pdf")
heatmap.2(as.matrix(x.all.tumVar), col=greenred(50), trace="none", density.info="none", hclustfun = function(x)hclust(x, method = "ward"), margins=c(10,10), cexRow=0.8, cexCol=0.7, key=T, ColSideColors=colsidecol)
dev.off()

#####################
#Heatmaps3
#####################
x <- pancreas.data.meth
rownames(x) <- x$Probes
x <- x[island.tss.probes, 3:ncol(x)]
dim(x)
x <- na.omit(x)
x.tumor <- x[,panc.tumors]
dim(x); dim(x.tumor)
x.sd <- apply(x, 1, sd)
x.tumor.sd <- apply(x.tumor, 1, sd)
dim(x); dim(x.tumor)
x.all <- x[which(x.sd > 0.2), ]
x.all.tumVar <- x[which(x.tumor.sd > 0.2), ]
dim(x.all); dim(x.all.tumVar)
colsidecol <- vector()
for (i in 1:ncol(x)) {colsidecol[i] <- ifelse(colnames(x)[i] %in% panc.tumors, "red", "green")}
pdf(file="./Robjects/x.all.cpgineartss.pdf")
heatmap.2(as.matrix(x.all), col=greenred(50), trace="none", density.info="none", hclustfun = function(x)hclust(x, method = "ward"), margins=c(10,10), cexRow=0.8, cexCol=0.7, key=T, ColSideColors=colsidecol)
dev.off()
pdf(file="./Robjects/x.all.tumVar.cpgineartss.pdf")
heatmap.2(as.matrix(x.all.tumVar), col=greenred(50), trace="none", density.info="none", hclustfun = function(x)hclust(x, method = "ward"), margins=c(10,10), cexRow=0.8, cexCol=0.7, key=T, ColSideColors=colsidecol)
dev.off()
#####################
#Heatmaps4 - filtered
#####################
index1 <- grep("IGF2", pancreas.data.meth[,2]); index2 <- grep("SOCS1", pancreas.data.meth[,2]); index3 <- grep("RUNX3", pancreas.data.meth[,2]); index4 <- grep("CACNA1G", pancreas.data.meth[,2]); index5 <- grep("NEUROG1", pancreas.data.meth[,2])
index.all <- sort(c(index1, index2, index3, index4, index5))
x <- pancreas.data.meth[index.all,]
rownames(x) <- x$Probes
x <- x[island.tss.probes, 3:ncol(x)]
dim(x)
x <- na.omit(x)
x.tumor <- x[,panc.tumors]
dim(x); dim(x.tumor)
x.sd <- apply(x, 1, sd)
x.tumor.sd <- apply(x.tumor, 1, sd)
dim(x); dim(x.tumor)
x.all <- x[which(x.sd > 0.2), ]
x.all.tumVar <- x[which(x.tumor.sd > 0.2), ]
dim(x.all); dim(x.all.tumVar)
colsidecol <- vector()
for (i in 1:ncol(x)) {colsidecol[i] <- ifelse(colnames(x)[i] %in% panc.tumors, "red", "green")}
pdf(file="./Robjects/x.filtered.cpgineartss.pdf")
heatmap.2(as.matrix(x.all), col=greenred(50), trace="none", density.info="none", hclustfun = function(x)hclust(x, method = "ward"), margins=c(10,10), cexRow=0.8, cexCol=0.7, key=T, ColSideColors=colsidecol)
dev.off()
pdf(file="./Robjects/x.filtered.tumVar.cpgineartss.pdf")
heatmap.2(as.matrix(x.all.tumVar), col=greenred(50), trace="none", density.info="none", hclustfun = function(x)hclust(x, method = "ward"), margins=c(10,10), cexRow=0.8, cexCol=0.7, key=T, ColSideColors=colsidecol)
dev.off()

#
#
#
#
#
#
#
#
#
#

#PCA with island-TSS 
rm(x)
x <- pancreas.data.meth[3:ncol(pancreas.data.meth)]
rownames(x) <- pancreas.data.meth$Probes
dim(x)
x <- x[intersect(autosome.probes, island.tss.probes), ]
dim(x)
x <- na.omit(x)
dim(x)
PCA <- prcomp(t(x))
pdf("./Robjects/pca-island_TSS-probes.pdf", height = 10, width = 10)
plot(PCA$x[, c("PC1", "PC2")], main="PCA of Panc Cancers colored by tumor/normal\nisland-TSS probes", xlab=paste("PC1: ",signif((((PCA$sdev)^2)[1]/(sum((PCA$sdev)^2))),2),sep=""), ylab=paste("PC2: ",signif((((PCA$sdev)^2)[2]/(sum((PCA$sdev)^2))),2),sep=""), pch=19, cex=0.5)
points(matrix(PCA$x[panc.tumors, c("PC1", "PC2")], ncol=2), pch=19, col=rgb(1, 0, 0, max = 1), cex=0.5)
points(matrix(PCA$x[panc.normals, c("PC1", "PC2")], ncol=2), pch=19, col=rgb(0, 1, 0, max = 1), cex=0.5)
legend("topright", legend=c("tumor", "normal"), col=c("red", "green"), cex=1, pch=19)
dev.off()


#Heatmaps (all variable probes)
x <- pancreas.data.meth[3:ncol(pancreas.data.meth)]
dim(x)
x <- na.omit(x)
x.tumor <- x[,panc.tumors]
dim(x); dim(x.tumor)
x.sd <- apply(x, 1, sd)
x.tumor.sd <- apply(x.tumor, 1, sd)
plot(density(x.sd))
plot(density(x.tumor.sd))
dim(x); dim(x.tumor)
x.all <- x[which(x.sd > 0.2), ]
x.all.tumVar <- x[which(x.tumor.sd > 0.2), ]
dim(x.all); dim(x.all.tumVar)
colsidecol <- vector()
for (i in 1:ncol(x)) {colsidecol[i] <- ifelse(colnames(x)[i] %in% panc.tumors, "red", "green")}
pdf(file="./Robjects/x.all.new.pdf")
heatmap.2(as.matrix(x.all), col=greenred(50), trace="none", density.info="none", hclustfun = function(x)hclust(x, method = "ward"), margins=c(10,10), cexRow=0.8, cexCol=0.7, key=T, ColSideColors=colsidecol)
dev.off()
pdf(file="./Robjects/x.all.tumVar.new.pdf")
heatmap.2(as.matrix(x.all.tumVar), col=greenred(50), trace="none", density.info="none", hclustfun = function(x)hclust(x, method = "ward"), margins=c(10,10), cexRow=0.8, cexCol=0.7, key=T, ColSideColors=colsidecol)
dev.off()

#Heatmap of only CpG-island probes near the TSS
rm(x)
x <- pancreas.data.meth[3:ncol(pancreas.data.meth)]
rownames(x) <- pancreas.data.meth$Probes
x <- x[intersect(autosome.probes, island.tss.probes), ]
x <- na.omit(x)
x.tumor <- x[,panc.tumors]
dim(x); dim(x.tumor)
x.sd <- apply(x, 1, sd)
x.tumor.sd <- apply(x.tumor, 1, sd)
# plot(density(x.sd))
# plot(density(x.tumor.sd))
dim(x); dim(x.tumor)
x.all <- x[which(x.sd > 0.15), ]
x.all.tumVar <- x[which(x.tumor.sd > 0.15), ]
dim(x.all); dim(x.all.tumVar)
colsidecol <- vector()
for (i in 1:ncol(x)) {colsidecol[i] <- ifelse(colnames(x)[i] %in% panc.tumors, "red", "green")}
## Heatmap of probes variable among tunors and normals
pdf("./Robjects/x.all.tss.new.pdf")
heatmap.2(as.matrix(x.all), col=greenred(50), trace="none", density.info="none", hclustfun = function(x)hclust(x, method = "ward.D"), margins=c(10,10), cexRow=0.8, cexCol=0.7, key=T, ColSideColors=colsidecol, main="Hierarchical clustering by probes variable among pancreatic tumors and normals\nCpG-island-TSS probes", cex.main=0.8)
dev.off()
## Heatmap of probes variable among tunors and normals
pdf("./Robjects/x.all.tumVar.tss.new.pdf")
heatmap.2(as.matrix(x.all.tumVar), col=greenred(50), trace="none", density.info="none", hclustfun = function(x)hclust(x, method = "ward.D"), margins=c(10,10), cexRow=0.8, cexCol=0.7, key=T, ColSideColors=colsidecol, main="Hierarchical clustering by probes variable among pancreatic tumors only\nCpG-island-TSS probes")
dev.off()

#Heatmap of only CpG-island probes near the TSS
#FILTERED
index1 <- grep("IGF2", pancreas.data.meth[,2]); index2 <- grep("SOCS1", pancreas.data.meth[,2]); index3 <- grep("RUNX3", pancreas.data.meth[,2]); index4 <- grep("CACNA1G", pancreas.data.meth[,2]); index5 <- grep("NEUROG1", pancreas.data.meth[,2])
index.all <- sort(c(index1, index2, index3, index4, index5))
cimp.probes <- pancreas.data.meth[index.all,1]
rm(x)
x <- pancreas.data.meth[3:ncol(pancreas.data.meth)]
rownames(x) <- pancreas.data.meth$Probes
x <- x[intersect(intersect(autosome.probes, island.tss.probes), cimp.probes), ]
x <- na.omit(x)
x.tumor <- x[,panc.tumors]
dim(x); dim(x.tumor)
x.sd <- apply(x, 1, sd)
x.tumor.sd <- apply(x.tumor, 1, sd)
# plot(density(x.sd))
# plot(density(x.tumor.sd))
dim(x); dim(x.tumor)
x.all <- x[which(x.sd > 0.15), ]
x.all.tumVar <- x[which(x.tumor.sd > 0.15), ]
dim(x.all); dim(x.all.tumVar)
colsidecol <- vector()
for (i in 1:ncol(x)) {colsidecol[i] <- ifelse(colnames(x)[i] %in% panc.tumors, "red", "green")}
## Heatmap of probes variable among tunors and normals
pdf("./Robjects/x.filtered.tss.new.pdf")
heatmap.2(as.matrix(x.all), col=greenred(50), trace="none", density.info="none", hclustfun = function(x)hclust(x, method = "ward.D"), margins=c(10,10), cexRow=0.8, cexCol=0.7, key=T, ColSideColors=colsidecol, main="Hierarchical clustering by probes variable among pancreatic tumors and normals\nCpG-island-TSS probes", cex.main=0.8)
dev.off()
## Heatmap of probes variable among tunors and normals
pdf("./Robjects/x.filtered.tumVar.tss.new.pdf")
heatmap.2(as.matrix(x.all.tumVar), col=greenred(50), trace="none", density.info="none", hclustfun = function(x)hclust(x, method = "ward.D"), margins=c(10,10), cexRow=0.8, cexCol=0.7, key=T, ColSideColors=colsidecol, main="Hierarchical clustering by probes variable among pancreatic tumors only\nCpG-island-TSS probes")
dev.off()