##########
# Purpose:
##########
# To create side annotation color bars

########
# Usage:
########
# Rscript create_ChIPseq_rowsidecolor_annotation_commandprompt.R "stable.10M" 200

Sys.time(); print(""); print("")

args <- commandArgs(trailingOnly=TRUE)
genelist.name <- as.character(args[1])
bin.size <- as.numeric(args[2]) #10bp or 200bp
print(paste0("genelist: ", genelist.name))
print(paste0("bin: ", bin.size))


dir<-"/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/"
# dir <- "Z:/users/shwang26/"
# dir<-"/amber2/scratch/baylin/shwang/"
source(paste0(dir,"Michelle/Rscripts/ChIP-SeqLibraryOfFunctions_original_newFunctions.R"))
# source(paste0(dir,"Michelle/Rscripts/ChIP-SeqLibraryOfFunctions_original.R"))
options(bitmapType='cairo') 
# coverage_files_dir = dir,"Michelle/normalizedBED"

coverage_files_dir = paste0(dir,"Michelle/BED_files/Coverage_TSS_", bin.size, "bp_bin/normalizedBED_", bin.size, "bp_bin/")
# plot_results_dir = dir,"Michelle/BED_Files/plots/normalized/15M_stable"
plot_results_dir = paste0(dir,"Michelle/BED_files/Coverage_TSS_", bin.size, "bp_bin/normalizedBED_", bin.size, "bp_bin/outputdir/")
dir.create(plot_results_dir)
setwd(plot_results_dir)
getwd()

# genelist.name="all"
if(genelist.name!="all"){
     version="hg19"
     if(genelist.name=="stable.10M"){
          genes <- read.table(paste0(system.dir,"Michelle/MethylationData/methylatedGeneLists/new/age_stable_methylated_genes_at_10months_new_annotation.txt"), header=FALSE)
     }
     #10M intermediate
     if(genelist.name=="intermediate.10M"){
          genes <- read.table(paste0(system.dir,"Michelle/MethylationData/methylatedGeneLists/new/age_intermediate_methylated_genes_at_10months_new_annotation.txt"), header=FALSE)
     }
     # stricter .4 cutoff
     #10M stable
     if(genelist.name=="stable.10M.new"){
          genes <- read.table(paste0(system.dir,"Michelle/MethylationData/methylatedGeneLists/new/age_stable_methylated_genes_at_10months_new_cutoff.txt"), header=FALSE)
     }
     #10M intermediate
     if(genelist.name=="intermediate.10M.new"){
          genes <- read.table(paste0(system.dir,"Michelle/MethylationData/methylatedGeneLists/new/age_intermediate_methylated_genes_at_10months_new_cutoff.txt"), header=FALSE)
     }
     #highly.expressed
     if(genelist.name=="highly.expressed"){
          genes <- read.table(paste0(system.dir,"Michelle/MethylationData/otherGenelists/list of high expression genes.txt"))
     }
     
     if(genelist.name=="viral.defense"){
          genes <- read.table(paste0(system.dir,"Michelle/MethylationData/otherGenelists/viraldefense.txt"))
     }
     if(genelist.name=="age.dependent"){
          genes <- read.table(paste0(system.dir,"Michelle/MethylationData/methylatedGeneLists/new/aging.dependent.methylated.genes_1M10M.txt"))
     }
     if(genelist.name=="stable.1M.newdata"){
          print("stable.1M.newdata genelist selected")
          genes <- read.table(paste0(system.dir,"Michelle/MethylationData/methylatedGeneLists/new/TSS_treatment.specific.hypermethylation_stable.1M.newdata_genelist.txt"))
     }
     if(genelist.name=="stable.6M.newdata"){
          print("stable.6M.newdata genelist selected")
          genes <- read.table(paste0(system.dir,"Michelle/MethylationData/methylatedGeneLists/new/TSS_treatment.specific.hypermethylation_stable.6M.newdata_genelist.txt"))
     }
     if(genelist.name=="stable.10M.newdata"){
          print("stable.10M.newdata genelist selected")
          genes <- read.table(paste0(system.dir,"Michelle/MethylationData/methylatedGeneLists/new/TSS_treatment.specific.hypermethylation_stable.10M.newdata_genelist.txt"))
     }
     if(genelist.name=="stable.15M.newdata"){
          print("stable.15M.newdata genelist selected")
          genes <- read.table(paste0(system.dir,"Michelle/MethylationData/methylatedGeneLists/new/TSS_treatment.specific.hypermethylation_stable.15M.newdata_genelist.txt"))
     }
     genes <- as.list(sort(unique(unlist(strsplit(unlist(as.matrix(genes)), ";")))))
     genes <- list(unlist(genes))
     names(genes) <- genelist.name
     goi <- fun.genelist.info(goi.list=genes, genelist.names=names(genes), goi.list.type="hgnc_symbol", version=version)
     goi <- lapply(c(1:length(goi)), FUN=function(i){
          x <- goi[[i]]
          return(x[!(duplicated(x[,"Gene"])), ])
     })
     names(goi) <-names(genes)
     goi
     goi[[1]]
}

# version = "hg19-UCSC"
if(genelist.name=="all"){
     version=="hg19-UCSC"
     hg19ucsc <- fun.genelist.info_allGenes_biomaRt(version="hg19-UCSC", 
                                                    hg19UCSCGeneAnnotations="hg19UCSCGeneAnnotaions.txt", 
                                                    hg19UCSCGeneAnnotationsPath=paste0(system.dir, "/Michelle/Required_Files/Annotations/hg19Data/"))
     dat <- cbind(hg19ucsc$Gene, hg19ucsc$Chromosome, hg19ucsc$TSS, hg19ucsc$TES, hg19ucsc$Strand, hg19ucsc$Gene)
     colnames(dat) <- c("Gene", "Chromosome", "TSS", "TES", "Strand", "ExternalGeneID")
     datlist <- list()
     datlist[[1]] <- as.data.frame(dat)
     goi <- datlist
     # goi[[1]][c(1:6),]
#      sample(1:nrow(goi[[1]]), 500)
     # nrow(goi[[1]])
     goi[[1]] <- goi[[1]][sample(1:nrow(goi[[1]]), 500),]
}

dim(goi[[1]])
# C10D
coverage_files = dir(coverage_files_dir)[grep("C10D", dir(coverage_files_dir))[1:6]]
input_coverage_file = coverage_files[6]
input_coverage_file
coverage_files

#metric
genelist.info=goi
coverage_files=coverage_files
input_coverage_file=input_coverage_file
coverage_files.dir=coverage_files_dir
RegAroundTSS=10000
bin=bin.size
chr.prefix.chromosome=T
plot.Directory=plot_results_dir
#plot.Directory=dir,"Michelle/BED_Files/plots"


library(gplots, lib.loc=lib.path)

if((RegAroundTSS/bin) %% 2 == 0){
  num.bins <- (RegAroundTSS/bin)+1
  int <- seq(from=-(RegAroundTSS/(bin*2)), to=(RegAroundTSS/(bin*2)), by=1)*bin
} else {
  num.bins <- (RegAroundTSS/bin)
  int <-  seq(from=-((num.bins-1)/2), to= (num.bins-1)/2, by=1)*bin
}
# Read the input coverage file
if(exists("input_coverage_file"))
{
  input.coverage <- read.table(file.path(coverage_files.dir, input_coverage_file), header=F, sep="\t")
  colnames(input.coverage) <- c("Chromosome", "start", "end", "Seq_tags")
} else 
{
  print("INPUT COVERAGE DATA ABSENT")
}

i <- dir(coverage_files.dir)[grep("^C10DH3K4", dir(coverage_files.dir))]
i

chip.coverage <- read.table(file.path(coverage_files.dir, i), header=F, sep="\t")
colnames(chip.coverage) <- c("Chromosome", "start", "end", "Seq_tags")
chipCoverage.TSSAverageList <- list() # holds the ChIP average values object for genes in genelist.info
inputCoverage.TSSAverageList <- list() # holds the Input average values object for genes in genelist.info
g=1 #for troubleshooting
bin
#********************CHECK chipseq-functions_original_newfunctions
chipCoverage.TSS <- fun.GetGeneCoverage(genelist.info=genelist.info[[g]], coverage.data=chip.coverage, RegAroundTSS=RegAroundTSS, bin=bin, chr.prefix.chromosome=chr.prefix.chromosome, version=version)
inputCoverage.TSS <- fun.GetGeneCoverage(genelist.info=genelist.info[[g]], coverage.data=input.coverage, RegAroundTSS=RegAroundTSS, bin=bin, chr.prefix.chromosome=chr.prefix.chromosome, version=version)

## Calculate average ChIP coverage
chipCoverage.TSSAverage <- apply(chipCoverage.TSS, 1, mean, na.rm=T)
chipCoverage.TSSAverageList[[g]] <- chipCoverage.TSSAverage
names(chipCoverage.TSSAverageList)[g] <- paste(i, names(genelist.info)[g], sep="-")
#chipCoverage.TSSAverageList
## Calculate average input coverage
inputCoverage.TSSAverage <- apply(inputCoverage.TSS, 1, mean, na.rm=T)
## If any of the inputCoverage.TSSAverage vales is 0, replace it with the minimum value
inputCoverage.TSSAverage <- replace(inputCoverage.TSSAverage, which(inputCoverage.TSSAverage == 0), min(inputCoverage.TSSAverage[which(inputCoverage.TSSAverage != 0)]))
inputCoverage.TSSAverageList[[g]] <- inputCoverage.TSSAverage
names(inputCoverage.TSSAverageList)[g] <- paste(i, names(genelist.info)[g], sep="-")

# Plot coverage (raw seq reads and ratio to input of seq reads) as heatmap plot
# Generate matices for heatmaps
x.trnposed <- t(chipCoverage.TSS)
ratioToInp.x.trnposed <- t(chipCoverage.TSS/inputCoverage.TSSAverage) #ratio of seq reads to average of input for this set of genes

colnames(x.trnposed) <- colnames(ratioToInp.x.trnposed) <- int
#write.table(int, file=paste0(system.dir, "Michelle/Required_Files/Annotations/chipseq.heatmap/int_",bin,"bp-2.txt"), row.names=FALSE, col.names=FALSE, sep="\t")
#rownames(x.trnposed) <- rownames(ratioToInp.x.trnposed) <- genelist.info[[g]]$hgnc_symbol
rownames(x.trnposed) <- rownames(ratioToInp.x.trnposed) <- genelist.info[[g]]$Gene
numberOfColors.x.trnposed <- colorRampPalette(c("white", "black"))(round(range(x.trnposed, na.rm=T)[2]))
numberOfColors.ratioToInp.x.trnposed <- colorRampPalette(c("white", "black"))(round(range(ratioToInp.x.trnposed, na.rm=T)[2]))

# get annotation_c10d_h3k4
dim(x.trnposed)
# trnposed200 <- x.trnposed
dim(x.trnposed[complete.cases(x.trnposed),])
dim(x.trnposed)
jpeg(paste0(genelist.name,"_heatmap.",bin.size,"bp_h3k4.test.jpeg"), height=600, width=900, quality=100)
# heatmap.x <- heatmap.2(x.trnposed, Rowv=T, Colv=F, scale="none", col=numberOfColors.x.trnposed, trace="none", dendrogram="row", cexRow=0.2, main="Heatmap of raw seq reads", na.rm=TRUE)
heatmap.x <- heatmap.3(x.trnposed[complete.cases(x.trnposed),], 
                       Rowv=T, 
                       Colv=F, 
                       scale="none", 
                       col=numberOfColors.x.trnposed, 
                       trace="none", 
                       dendrogram="row", 
                       cexRow=1, 
                       main="Heatmap of raw seq reads", na.rm=TRUE)
dev.off()


# heatmap.x <- heatmap.3(x.trnposed)
names(heatmap.x)
heatmap.x$rowDendrogram
heatmap.x$rowInd
hm <- as.hclust(heatmap.x$rowDendrogram)
#names(hm)
#hm$order
cut <- cutree(hm, k=4)
#dendrogram.order <- cut[hm$order]
dendrogram.order <- cut
length(dendrogram.order)
dendrogram.order[1]
row.annotation.color <- vector()
# for (i in 1:length(dendrogram.order)) {
# row.annotation.color[i] <- ifelse(dendrogram.order[i]==2, "green4", "red4")
# }
for (i in 1:length(dendrogram.order)) {
  if(dendrogram.order[i]==1){
    row.annotation.color[i] <- "gray10"
  }
  if(dendrogram.order[i]==2){
    row.annotation.color[i] <- "green3"
  }
  if(dendrogram.order[i]==3){
    row.annotation.color[i] <- "dodgerblue"
  }
  if(dendrogram.order[i]==4){
    row.annotation.color[i] <- "red"
  }
}

row.annotation.color <- t(as.matrix(row.annotation.color))
dim(row.annotation.color)
colnames(row.annotation.color) <- names(dendrogram.order)
rownames(row.annotation.color) <- "C10D_H3K4"
# dim(row.annotation.color)
# head(row.annotation.color)
h3k4.annotation <- row.annotation.color
x.trnposed_c10d.h3k4 <- x.trnposed[complete.cases(x.trnposed),]

#CHECK

heatmap.3(x.trnposed_c10d.h3k4[complete.cases(x.trnposed_c10d.h3k4),], Rowv=T, Colv=F, scale="none", col=numberOfColors.x.trnposed, RowSideColors=h3k4.annotation, trace="none", dendrogram="row", cexRow=0.2, main="Heatmap of raw seq reads", na.rm=TRUE)

#####################################################################################################
#####################################################################################################

# H3K27
# C10D
i <- dir(coverage_files.dir)[grep("^C10DH3K27", dir(coverage_files.dir))]
i

chip.coverage <- read.table(file.path(coverage_files.dir, i), header=F, sep="\t")
colnames(chip.coverage) <- c("Chromosome", "start", "end", "Seq_tags")
chipCoverage.TSSAverageList <- list() # holds the ChIP average values object for genes in genelist.info
inputCoverage.TSSAverageList <- list() # holds the Input average values object for genes in genelist.info
g=1 #for troubleshooting
chipCoverage.TSS <- fun.GetGeneCoverage(genelist.info=genelist.info[[g]], coverage.data=chip.coverage, RegAroundTSS=RegAroundTSS, bin=bin, chr.prefix.chromosome=chr.prefix.chromosome, version=version)
inputCoverage.TSS <- fun.GetGeneCoverage(genelist.info=genelist.info[[g]], coverage.data=input.coverage, RegAroundTSS=RegAroundTSS, bin=bin, chr.prefix.chromosome=chr.prefix.chromosome, version=version)
## Calculate average ChIP coverage
chipCoverage.TSSAverage <- apply(chipCoverage.TSS, 1, mean, na.rm=T)
chipCoverage.TSSAverageList[[g]] <- chipCoverage.TSSAverage
names(chipCoverage.TSSAverageList)[g] <- paste(i, names(genelist.info)[g], sep="-")
#chipCoverage.TSSAverageList
## Calculate average input coverage
inputCoverage.TSSAverage <- apply(inputCoverage.TSS, 1, mean, na.rm=T)
## If any of the inputCoverage.TSSAverage vales is 0, replace it with the minimum value
inputCoverage.TSSAverage <- replace(inputCoverage.TSSAverage, which(inputCoverage.TSSAverage == 0), min(inputCoverage.TSSAverage[which(inputCoverage.TSSAverage != 0)]))
inputCoverage.TSSAverageList[[g]] <- inputCoverage.TSSAverage
names(inputCoverage.TSSAverageList)[g] <- paste(i, names(genelist.info)[g], sep="-")

# Plot coverage (raw seq reads and ratio to input of seq reads) as heatmap plot
## Create plot name for raw heatmaps
# Generate matices for heatmaps
x.trnposed <- t(chipCoverage.TSS)
ratioToInp.x.trnposed <- t(chipCoverage.TSS/inputCoverage.TSSAverage) #ratio of seq reads to average of input for this set of genes

colnames(x.trnposed) <- colnames(ratioToInp.x.trnposed) <- int
rownames(x.trnposed) <- rownames(ratioToInp.x.trnposed) <- genelist.info[[g]]$hgnc_symbol
numberOfColors.x.trnposed <- colorRampPalette(c("white", "black"))(round(range(x.trnposed, na.rm=T)[2]))
numberOfColors.ratioToInp.x.trnposed <- colorRampPalette(c("white", "black"))(round(range(ratioToInp.x.trnposed, na.rm=T)[2]))

# get annotation_c10d_h3k27
jpeg(paste0(genelist.name,"_heatmap.",bin.size,"bp_h3k27.test.jpeg"), height=600, width=900, quality=100)
heatmap.x <- heatmap.2(x.trnposed[complete.cases(x.trnposed),], Rowv=T, Colv=F, scale="none", col=numberOfColors.x.trnposed, trace="none", dendrogram="row", cexRow=0.2, main="Heatmap of raw seq reads", na.rm=TRUE)
dev.off()
names(heatmap.x)
heatmap.x$rowDendrogram
heatmap.x$rowInd
hm <- as.hclust(heatmap.x$rowDendrogram)
#names(hm)
#hm$order
cut <- cutree(hm, k=4)
#dendrogram.order <- cut[hm$order]
dendrogram.order <- cut
length(dendrogram.order)
dendrogram.order[1]
row.annotation.color <- vector()
# for (i in 1:length(dendrogram.order)) {
#   row.annotation.color[i] <- ifelse(dendrogram.order[i]==2, "green4", "red4")
# }

for (i in 1:length(dendrogram.order)) {
  if(dendrogram.order[i]==1){
    row.annotation.color[i] <- "gray10"
  }
  if(dendrogram.order[i]==2){
    row.annotation.color[i] <- "green3"
  }
  if(dendrogram.order[i]==3){
    row.annotation.color[i] <- "dodgerblue"
  }
  if(dendrogram.order[i]==4){
    row.annotation.color[i] <- "red"
  }
}
row.annotation.color

row.annotation.color <- t(as.matrix(row.annotation.color))
dim(row.annotation.color)
colnames(row.annotation.color) <- names(dendrogram.order)
rownames(row.annotation.color) <- "C10D_H3K27"
# dim(row.annotation.color)
# head(row.annotation.color)
h3k27.annotation <- row.annotation.color
x.trnposed_c10d.h3k27 <- x.trnposed[complete.cases(x.trnposed),]

#CHECK
heatmap.3(x.trnposed_c10d.h3k27, Rowv=T, Colv=F, scale="none", col=numberOfColors.x.trnposed, RowSideColors=h3k27.annotation, trace="none", dendrogram="row", cexRow=.2, main="Heatmap of raw seq reads", na.rm=TRUE)

c10d.10M.ann.bar <- rbind(h3k4.annotation, h3k27.annotation)
head(c10d.10M.ann.bar)

# NEW DOUBLE BAR
getwd()
jpeg(paste0("C10D.h3k4_",names(genes),"_",bin,"bp_annotation.calibration2.jpeg"), height=900, width=1600, quality=100)
par(mar=c(10,5,5,5))
heatmap.3(x.trnposed_c10d.h3k4, Rowv=T, Colv=F, scale="none", col=numberOfColors.x.trnposed, RowSideColors=c10d.10M.ann.bar , trace="none", dendrogram="row", cexRow=3, main="Heatmap of raw seq reads", na.rm=TRUE, cex=2)
dev.off()

jpeg(paste0("C10D.h3k27_",names(genes),"_",bin,"bp_annotation.calibration2.jpeg"), height=900, width=1600, quality=100)
par(mar=c(10,5,5,5))
heatmap.3(x.trnposed_c10d.h3k27, Rowv=T, Colv=F, scale="none", col=numberOfColors.x.trnposed, RowSideColors=c10d.10M.ann.bar , trace="none", dendrogram="row", cexRow=3, main="Heatmap of raw seq reads", na.rm=TRUE, cex=2)
dev.off()

# save(c10d.10M.ann.bar, file=paste0("c10d_",names(genes),"_",bin,"bp_ann.bar.Rdata"))
paste0("c10d_",genelist.name,"_",bin,"bp_ann.bar.Rdata")
save(c10d.10M.ann.bar, file=paste0("c10d_",genelist.name,"_",bin,"bp_ann.bar.Rdata"))

heatmap.3(x.trnposed_c10d.h3k27, Rowv=T, Colv=F, scale="none", col=numberOfColors.x.trnposed, RowSideColors=c10d.10M.ann.bar, trace="none", dendrogram="row", cexRow=.2, main="Heatmap of raw seq reads", na.rm=TRUE)

