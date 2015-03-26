#/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/
#/amber2/scratch/baylin/shwang/
#/Volumes/onc-analysis$/users/stephenhwang/

bin.size=10
source("/amber2/scratch/baylin/shwang/Michelle/Rscripts/ChIP-SeqLibraryOfFunctions_original_newFunctions.R")
options(bitmapType='cairo') 
# coverage_files_dir = "/amber2/scratch/baylin/shwang/Michelle/normalizedBED"

coverage_files_dir = paste0("/amber2/scratch/baylin/shwang/Michelle/BED_files/Coverage_TSS_", bin.size, "bp_bin/normalizedBED_", bin.size, "bp_bin/")
# plot_results_dir = "/amber2/scratch/baylin/shwang/Michelle/BED_Files/plots/normalized/15M_stable"
plot_results_dir = paste0("/amber2/scratch/baylin/shwang/Michelle/BED_files/Coverage_TSS_", bin.size, "bp_bin/normalizedBED_", bin.size, "bp_bin/outputdir/")
setwd(plot_results_dir)

######################################################################################################

#10M stable
genes.name <- "stable.10M"
genes <- read.table("/amber2/scratch/baylin/shwang/Michelle/MethylationData/methylatedGeneLists/new/age_stable_methylated_genes_at_10months_new_annotation.txt", header=FALSE)

#10M intermediate
genes.name <- "intermediate.10M"
genes <- read.table("/amber2/scratch/baylin/shwang/Michelle/MethylationData/methylatedGeneLists/new/age_intermediate_methylated_genes_at_10months_new_annotation.txt", header=FALSE)

# stricter .4 cutoff
#10M stable
genes.name <- "new_stable.10M"
genes <- read.table("/amber2/scratch/baylin/shwang/Michelle/MethylationData/methylatedGeneLists/new/age_stable_methylated_genes_at_10months_new_cutoff.txt", header=FALSE)

#10M intermediate
genes.name <- "new_intermediate.10M"
genes <- read.table("/amber2/scratch/baylin/shwang/Michelle/MethylationData/methylatedGeneLists/new/age_intermediate_methylated_genes_at_10months_new_cutoff.txt", header=FALSE)

######################################################################################################

genes <- as.list(sort(unique(unlist(strsplit(unlist(as.matrix(genes)), ";")))))
genes <- list(unlist(genes))
names(genes) <- genes.name

version = "hg19"
goi <- fun.genelist.info(goi.list=genes, genelist.names=names(genes), goi.list.type="hgnc_symbol", version=version)
goi <- lapply(c(1:length(goi)), FUN=function(i){
  x <- goi[[i]]
  return(x[!(duplicated(x[,"Gene"])), ])
})
names(goi) <-names(genes)
goi

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
#plot.Directory="/amber2/scratch/baylin/shwang/Michelle/BED_Files/plots"


library(gplots, lib.loc=lib.path)
dir.create(plot.Directory)
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
## Create plot name for raw heatmaps
heat_plot_raw <- file.path(plot.Directory, paste(i, names(genelist.info)[g], "heat_plot_raw.jpeg", sep="-"))
## Create plot name for RATIO heatmaps
heat_plot_RATIO <- file.path(plot.Directory, paste(i, names(genelist.info)[g], "heat_plot_RATIO.jpeg", sep="-"))

# Generate matices for heatmaps
x.trnposed <- t(chipCoverage.TSS)
ratioToInp.x.trnposed <- t(chipCoverage.TSS/inputCoverage.TSSAverage) #ratio of seq reads to average of input for this set of genes

colnames(x.trnposed) <- colnames(ratioToInp.x.trnposed) <- int
#rownames(x.trnposed) <- rownames(ratioToInp.x.trnposed) <- genelist.info[[g]]$hgnc_symbol
rownames(x.trnposed) <- rownames(ratioToInp.x.trnposed) <- genelist.info[[g]]$Gene
numberOfColors.x.trnposed <- colorRampPalette(c("white", "black"))(round(range(x.trnposed, na.rm=T)[2]))
numberOfColors.ratioToInp.x.trnposed <- colorRampPalette(c("white", "black"))(round(range(ratioToInp.x.trnposed, na.rm=T)[2]))

# get annotation_c10d_h3k4
jpeg("test_heatmap.x10bp.jpeg", height=600, width=900, quality=100)
# heatmap.x <- heatmap.2(x.trnposed, Rowv=T, Colv=F, scale="none", col=numberOfColors.x.trnposed, trace="none", dendrogram="row", cexRow=0.2, main="Heatmap of raw seq reads", na.rm=TRUE)
heatmap.x <- heatmap.3(x.trnposed, Rowv=T, Colv=F, scale="none", col=numberOfColors.x.trnposed, trace="none", dendrogram="row", cexRow=0.2, main="Heatmap of raw seq reads", na.rm=TRUE)
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
dim(row.annotation.color)
head(row.annotation.color)
h3k4.annotation <- row.annotation.color
x.trnposed_c10d.h3k4 <- x.trnposed

#CHECK
heatmap.3(x.trnposed_c10d.h3k4, Rowv=T, Colv=F, scale="none", col=numberOfColors.x.trnposed, RowSideColors=h3k4.annotation, trace="none", dendrogram="row", cexRow=0.2, main="Heatmap of raw seq reads", na.rm=TRUE)

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
heat_plot_raw <- file.path(plot.Directory, paste(i, names(genelist.info)[g], "heat_plot_raw.jpeg", sep="-"))
#heat_plot_raw <- file.path(plot.Directory, paste(i, names(genelist.info)[g], "x_raw.jpeg", sep="-"))
## Create plot name for RATIO heatmaps
heat_plot_RATIO <- file.path(plot.Directory, paste(i, names(genelist.info)[g], "heat_plot_RATIO.jpeg", sep="-"))
#heat_plot_RATIO <- file.path(plot.Directory, paste(i, names(genelist.info)[g], "x_RATIO.jpeg", sep="-"))

# Generate matices for heatmaps
x.trnposed <- t(chipCoverage.TSS)
ratioToInp.x.trnposed <- t(chipCoverage.TSS/inputCoverage.TSSAverage) #ratio of seq reads to average of input for this set of genes

colnames(x.trnposed) <- colnames(ratioToInp.x.trnposed) <- int
rownames(x.trnposed) <- rownames(ratioToInp.x.trnposed) <- genelist.info[[g]]$hgnc_symbol
numberOfColors.x.trnposed <- colorRampPalette(c("white", "black"))(round(range(x.trnposed, na.rm=T)[2]))
numberOfColors.ratioToInp.x.trnposed <- colorRampPalette(c("white", "black"))(round(range(ratioToInp.x.trnposed, na.rm=T)[2]))

# get annotation_c10d_h3k27
jpeg("test_heatmap.x27.jpeg", height=600, width=900, quality=100)
heatmap.x <- heatmap.2(x.trnposed, Rowv=T, Colv=F, scale="none", col=numberOfColors.x.trnposed, trace="none", dendrogram="row", cexRow=0.2, main="Heatmap of raw seq reads", na.rm=TRUE)
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
dim(row.annotation.color)
head(row.annotation.color)
h3k27.annotation <- row.annotation.color
x.trnposed_c10d.h3k27 <- x.trnposed

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
save(c10d.10M.ann.bar, file=paste0("c10d_",genes.name,"_",bin,"bp_ann.bar.Rdata"))