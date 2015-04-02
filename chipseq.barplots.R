################################################################################
################################################################################
# this script will generate the chipseq barplots + txt files containing genes
# that are part of each barplot group
################################################################################
################################################################################
system.dir="/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/"
# system.dir="/Volumes/onc-analysis$/stephenhwang/"
# system.dir="/amber2/scratch/baylin/shwang/"
# system.dir="Z:/users/shwang26/"
# source(paste0(system.dir, "Michelle/Rscripts/ChIP-SeqLibraryOfFunctions_original.R"))
source(paste0(system.dir, "Michelle/Rscripts/ChIP-SeqLibraryOfFunctions_original_newFunctions.R"))
date <- gsub(" \\d+\\:\\d+\\:\\d+.*", "", Sys.time())
newdir <- paste0(system.dir, "Michelle/Analysis/chipseq_plot/", date)
newdir
# Create new working directory
dir.create(newdir)
setwd(newdir)
getwd()

################################################################################
################################################################################
# track which directories files are being saved to
################################################################################
################################################################################
log.line <- paste0(Sys.time(), "\t", getwd())
log.line
write(log.line, file=paste0(system.dir, "Michelle/save_log.txt"), append=TRUE)

################################################################################
################################################################################
# pull in genelist
################################################################################
################################################################################

# Where is hg19UCSCGeneAnnotaions.txt generated from?
genelist.info <- fun.genelist.info_allGenes_biomaRt(version="hg19-UCSC", hg19UCSCGeneAnnotations="hg19UCSCGeneAnnotaions.txt", hg19UCSCGeneAnnotationsPath=paste0(system.dir, "Michelle/Required_Files/Annotations/hg19Data"))

# load in genelists
load(file=paste0(system.dir, "Michelle/Robjects/genelists.rda"))

################################################################################
################################################################################
# generate plots
################################################################################
################################################################################
methylatedGeneListsDir=paste0(system.dir, "Michelle/MethylationData/methylatedGeneLists")


fun.find.unique.genes.between.group(C10D_EZH2_genes, CSC10D_EZH2_genes, "C10D", "CSC10D_EZH2")
fun.find.unique.genes.between.group(C3M_EZH2_genes, CSC3M_EZH2_genes, "C3M", "CSC3M_EZH2")
fun.find.unique.genes.between.group(C10M_EZH2_genes, CSC10M_EZH2_genes, "C10M", "CSC10M_EZH2")
fun.find.unique.genes.between.group(C10D_K4.K27_genes, CSC10D_K4.K27_genes, "C10D", "CSC10D_K4.K27")
fun.find.unique.genes.between.group(C3M_K4.K27_genes, CSC3M_K4.K27_genes, "C3M", "CSC3M_K4.K27")
fun.find.unique.genes.between.group(C10M_K4.K27_genes, CSC10M_K4.K27_genes, "C10M", "CSC10M_K4.K27")

# create barplots
for(i in dir(methylatedGeneListsDir)){ # for each file
     tab=read.table(file.path(methylatedGeneListsDir, i), header=F, sep="\t", stringsAsFactors=F) # read in
     # Remove duplication of gene names from the list Michelle provided
     for(j in grep(";", tab[,1])){
          tab[j,1] = strsplit(tab[j,1], split=";")[[1]][1]
     }  
     
     dir.create(paste0(getwd(), "/all3_10M_Comparisons"))
     
     fun.BarPlot(set1=C10D_K4.K27.genelist.info$Gene, 
                 set2=CSC10D_K4.K27.genelist.info$Gene, 
                 set3=C10M_K4.K27.genelist.info$Gene, 
                 set4=CSC10M_K4.K27.genelist.info$Gene, 
                 set5=tab[,1], 
                 set6=C3M_K4.K27.genelist.info$Gene, 
                 set7=CSC3M_K4.K27.genelist.info$Gene, 
                 set1Name="C10D_K4.K27", 
                 set2Name="CSC10D_K4.K27", 
                 set3Name="C10M_K4.K27", 
                 set4Name="CSC10M_K4.K27", 
                 set5Name=i, 
                 set6Name="C3M_K4.K27", 
                 set7Name="CSC3M_K4.K27", 
                 barPlotName=paste(i, "barPlotK4K27", sep="_"), 
                 plotDir=paste0(system.dir, "Michelle/Analysis/Plots/all3_Comparisons"))
     
     fun.BarPlot(set1=C10D_GenesWithPeaks_K4$Gene, 
                 set2=CSC10D_GenesWithPeaks_K4$Gene, 
                 set3=C10M_GenesWithPeaks_K4$Gene, 
                 set4=CSC10M_GenesWithPeaks_K4$Gene, 
                 set5=tab[,1], 
                 set6=C3M_GenesWithPeaks_K4$Gene, 
                 set7=CSC3M_GenesWithPeaks_K4$Gene, 
                 set1Name="C10D_K4", 
                 set2Name="CSC10D_K4", 
                 set3Name="C10M_K4", 
                 set4Name="CSC10M_K4", 
                 set5Name=i, 
                 set6Name="C3M_K4", 
                 set7Name="CSC3M_K4", 
                 barPlotName=paste(i, "barPlotK4", sep="_"), 
                 plotDir=paste0(system.dir, "Michelle/Analysis/Plots/all3_Comparisons"))
     
     fun.BarPlot(set1=C10D_GenesWithPeaks_K27$Gene, 
                 set2=CSC10D_GenesWithPeaks_K27$Gene, 
                 set3=C10M_GenesWithPeaks_K27$Gene, 
                 set4=CSC10M_GenesWithPeaks_K27$Gene, 
                 set5=tab[,1], 
                 set6=C3M_GenesWithPeaks_K27$Gene, 
                 set7=CSC3M_GenesWithPeaks_K27$Gene, 
                 set1Name="C10D_K27", 
                 set2Name="CSC10D_K27", 
                 set3Name="C10M_K27", 
                 set4Name="CSC10M_K27", 
                 set5Name=i, 
                 set6Name="C3M_K27", 
                 set7Name="CSC3M_K27", 
                 barPlotName=paste(i, "barPlotK27", sep="_"), 
                 plotDir=paste0(dir,"Michelle/Analysis/Plots/all3_Comparisons"))
     
     fun.BarPlot(set1=C10D_GenesWithPeaks_EZH2$Gene, 
                 set2=CSC10D_GenesWithPeaks_EZH2$Gene, 
                 set3=C10M_GenesWithPeaks_EZH2$Gene, 
                 set4=CSC10M_GenesWithPeaks_EZH2$Gene, 
                 set5=tab[,1], 
                 set6=C3M_GenesWithPeaks_EZH2$Gene, 
                 set7=CSC3M_GenesWithPeaks_EZH2$Gene, 
                 set1Name="C10D_EZH2", 
                 set2Name="CSC10D_EZH2", 
                 set3Name="C10M_EZH2", 
                 set4Name="CSC10M_EZH2", 
                 set5Name=i, 
                 set6Name="C3M_EZH2", 
                 set7Name="CSC3M_EZH2", 
                 barPlotName=paste(i, "barPlotEZH2", sep="_"), 
                 plotDir=paste0(dir,"Michelle/Analysis/Plots/all3_Comparisons")
}

# generate barplots based off all comparisons
for(i in dir(methylatedGeneListsDir)){ # for each file
     tab=read.table(file.path(methylatedGeneListsDir, i), header=F, sep="\t", stringsAsFactors=F) # read in
     # Remove duplication of gene names from the list Michelle provided
     for(j in grep(";", tab[,1])){
          tab[j,1] = strsplit(tab[j,1], split=";")[[1]][1]
     }  
     dir.create(paste0(getwd(), "/all3_10M_Comparisons"))
     fun.allComparisons(set1=C10D_GenesWithPeaks_K4$Gene, 
                        set2=C3M_GenesWithPeaks_K4$Gene, 
                        set3=C10M_GenesWithPeaks_K4$Gene, 
                        set4=CSC10D_GenesWithPeaks_K4$Gene, 
                        set5=tab[,1], 
                        set6=C3M_GenesWithPeaks_K4$Gene, 
                        set7=CSC3M_GenesWithPeaks_K4$Gene, 
                        set8=C10D_GenesWithPeaks_K27$Gene, 
                        set9=C3M_GenesWithPeaks_K27$Gene, 
                        set10=C10M_GenesWithPeaks_K27$Gene, 
                        set11=CSC10D_GenesWithPeaks_K27$Gene, 
                        set12=C3M_GenesWithPeaks_K27$Gene, 
                        set13=CSC3M_GenesWithPeaks_K27$Gene, 
                        set14=C10D_GenesWithPeaks_EZH2$Gene, 
                        set15=C3M_GenesWithPeaks_EZH2$Gene, 
                        set16=C10M_GenesWithPeaks_EZH2$Gene, 
                        set17=CSC10D_GenesWithPeaks_EZH2$Gene, 
                        set18=C3M_GenesWithPeaks_EZH2$Gene, 
                        set19=CSC3M_GenesWithPeaks_EZH2$Gene, 
                        set1Name="C10D_K4", 
                        set2Name="CSC10D_K4", 
                        set3Name="C10M_K4", 
                        set4Name="CSC10M_K4", 
                        set5Name=i, 
                        set6Name="C3M_K4", 
                        set7Name="CSC3M_K4", 
                        set8Name="C10D_K27", 
                        set9Name="CSC10D_K27", 
                        set10Name="C10M_K27", 
                        set11Name="CSC10M_K27", 
                        set12Name="C3M_K27", 
                        set13Name="CSC3M_K27", 
                        set14Name="C10D_EZH2", 
                        set15Name="CSC10D_EZH2", 
                        set16Name="C10M_EZH2", 
                        set17Name="CSC10M_EZH2", 
                        set18Name="C3M_EZH2", 
                        set19Name="CSC3M_EZH2", 
                        barPlotName=paste(i, "barPlot_allMarkers", sep="_"), 
                        plotDir=paste0(dir,"Michelle/Analysis/Plots/all3_Comparisons")
}


################################################################################
################################################################################
# custom barplots
################################################################################
################################################################################

methylatedGeneListsDir=paste0(system.dir,"Michelle/MethylationData/methylatedGeneLists/new/")
ls()
# Create barplots for 
dir(methylatedGeneListsDir)
genelists <- c("intermediate.10M", "intermediate.10M.new", "stable.10M", "stable.10M.new")
filename <- dir(methylatedGeneListsDir)[2]
genelist.name <- genelists[2]
filename
genelist.name
meth.genes=read.table(file.path(methylatedGeneListsDir, filename), header=F, sep="\t", stringsAsFactors=F) # read in
# Remove duplication of gene names from the list Michelle provided
for(j in grep(";", meth.genes[,1])){
     meth.genes[j,1] = strsplit(meth.genes[j,1], split=";")[[1]][1]
}  
# file.name <- gsub("months.*", "months", gsub("_methylated_genes_at", "", gsub("age_", "", filename)))
# file.name
meth.genes <- as.matrix(meth.genes)
# colnames(meth.genes) <- file.name
colnames(meth.genes) <- genelist.name
head(meth.genes)
# all genes
mark = cbind("EZH2", "H3K4", "H3K27", "K4.K27")
# mark <- mark[2]
# mark

# pause <- function(x)
# {
#      p1 <- proc.time()
#      Sys.sleep(x)
#      proc.time() - p1 # The cpu usage should be negligible
# }
# pause(3.7)

for(i in 1:length(mark)){
     if(mark[i]=="EZH2")     {
          print("EZH2 hit")
          s1 <- unique(C10D_EZH2_genes)
          s2 <- unique(CSC10D_EZH2_genes)
          s3 <- unique(C3M_EZH2_genes)
          s4 <- unique(CSC3M_EZH2_genes)
          s5 <- unique(C10M_EZH2_genes)
          s6 <- unique(CSC10M_EZH2_genes)
     }
     if(mark[i]=="H3K4"){
          print("H3K4 hit")
          s1 <- unique(C10D_K4_genes)
          s2 <- unique(CSC10D_K4_genes)
          s3 <- unique(C3M_K4_genes)
          s4 <- unique(CSC3M_K4_genes)
          s5 <- unique(C10M_K4_genes)
          s6 <- unique(CSC10M_K4_genes)          
     }
     if(mark[i]=="H3K27"){
          print("H3K27 hit")
          s1 <- unique(C10D_K27_genes)
          s2 <- unique(CSC10D_K27_genes)
          s3 <- unique(C3M_K27_genes)
          s4 <- unique(CSC3M_K27_genes)
          s5 <- unique(C10M_K27_genes)
          s6 <- unique(CSC10M_K27_genes)          
     }
     if(mark[i]=="K4.K27"){
          print("K4.K27 hit")
          s1 <- unique(C10D_K4.K27_genes)
          s2 <- unique(CSC10D_K4.K27_genes)
          s3 <- unique(C3M_K4.K27_genes)
          s4 <- unique(CSC3M_K4.K27_genes)
          s5 <- unique(C10M_K4.K27_genes)
          s6 <- unique(CSC10M_K4.K27_genes)
     }
     name1 <- "C10D"
     name2 <- "CSC10D"
     name3 <- "C3M"
     name4 <- "CSC3M"
     name5 <- "C10M"
     name6 <- "CSC10M"
     ################################
     # all.genes
     ################################
     # num of genes per sample
     a1 <-length(s1)
     a2 <-length(s2)
     a3 <-length(s3)
     a4 <-length(s4)
     a5 <-length(s5)
     a6 <-length(s6)
     # genes common to control, csc genes
     all.con_genes <- intersect(s5, intersect(s1,s3))
     all.csc_genes <- intersect(s6, intersect(s2,s4))
     all.con <- length(all.con_genes)
     all.csc <- length(all.csc_genes)
     # genes that are unique to that particular CON timepoint
     unique.to.s1 <- setdiff(s1, all.con_genes) #c10d
     unique.to.s2 <- setdiff(s2, all.csc_genes) #csc10d
     unique.to.s3 <- setdiff(s3, all.con_genes) #c3m
     unique.to.s4 <- setdiff(s4, all.csc_genes) #csc3m
     unique.to.s5 <- setdiff(s5, all.con_genes) #c10m
     unique.to.s6 <- setdiff(s6, all.csc_genes) #csc10m
     
     # num of genes unique to control/csc genes
     unique1 <- a1-all.con
     unique2 <- a2-all.csc
     unique3 <- a3-all.con
     unique4 <- a4-all.csc
     unique5 <- a5-all.con
     unique6 <- a6-all.csc
     percentage.1 <- paste0(round((unique1/a1)*100,1), "%")
     percentage.2 <- paste0(round((unique2/a2)*100,1), "%")
     percentage.3 <- paste0(round((unique3/a3)*100,1), "%")
     percentage.4 <- paste0(round((unique4/a4)*100,1), "%")
     percentage.5 <- paste0(round((unique5/a5)*100,1), "%")
     percentage.6 <- paste0(round((unique6/a6)*100,1), "%")
     percentage.1
     percentage.2
     percentage.3
     percentage.4
     percentage.5
     percentage.6
     ################################
     # meth.genes
     ################################
     
     # genes common to genelist and meth.genelist
     m1 <- intersect(meth.genes, s1)
     m2 <- intersect(meth.genes, s2)
     m3 <- intersect(meth.genes, s3)
     m4 <- intersect(meth.genes, s4)
     m5 <- intersect(meth.genes, s5)
     m6 <- intersect(meth.genes, s6)
     # num of genes common to genelist and meth.genelist
     n1=length(m1)
     n2=length(m2)
     n3=length(m3)
     n4=length(m4)
     n5=length(m5)
     n6=length(m6)
     # num of methylated.genes common to control/csc timepoints
     m.con_genes <- intersect(m5, intersect(m1,m3))
     m.csc_genes <- intersect(m6, intersect(m2,m4))
     m.con <- length(m.con_genes)
     m.csc <- length(m.csc_genes)
     
     # num of genes unique to 
     unique.m1 <- n1-m.con
     unique.m2 <- n2-m.csc
     unique.m3 <- n3-m.con
     unique.m4 <- n4-m.csc
     unique.m5 <- n5-m.con
     unique.m6 <- n6-m.csc
     # find the methylated.genes that are unique to that control sample timepoint
     unique.meth.genes.to.1 <- setdiff(m1, m.con_genes) #c10d
     unique.meth.genes.to.2 <- setdiff(m2, m.csc_genes) #csc10d
     unique.meth.genes.to.3 <- setdiff(m3, m.csc_genes) #c3m
     unique.meth.genes.to.4 <- setdiff(m4, m.csc_genes) #csc3m
     unique.meth.genes.to.5 <- setdiff(m5, m.csc_genes) #c10m
     unique.meth.genes.to.6 <- setdiff(m6, m.csc_genes) #csc10m
     percentage.m1 <- paste0(round((unique.m1/n1)*100,1), "%")
     percentage.m2 <- paste0(round((unique.m2/n2)*100,1), "%")
     percentage.m3 <- paste0(round((unique.m3/n3)*100,1), "%")
     percentage.m4 <- paste0(round((unique.m4/n4)*100,1), "%")
     percentage.m5 <- paste0(round((unique.m5/n5)*100,1), "%")
     percentage.m6 <- paste0(round((unique.m6/n6)*100,1), "%")
     percentage.m1
     percentage.m2
     percentage.m3
     percentage.m4
     percentage.m5
     percentage.m6
     
     # writeout common all.genes across all timepoints
     write.table(c(paste0("# common.all.genes_CONTROL_", mark[i], "_" ,genelist.name), sort(all.con_genes)), file=paste0("common.all.genes_CONTROL_", mark[i], "_" , genelist.name, ".txt"), row.names=FALSE, col.names=FALSE, quote=FALSE)
     write.table(c(paste0("# common.all.genes_CSC_", mark[i], "_" ,genelist.name), sort(all.csc_genes)), file=paste0("common.all.genes_CSC_", mark[i], "_" , genelist.name, ".txt"), row.names=FALSE, col.names=FALSE, quote=FALSE)
     
     # writeout unique all.genes at each timepoint
     write.table(c(paste0("# unique.all.genes_", mark[i], "_" ,genelist.name, "_", name1), sort(unique.to.s1)), file=paste0("unique.all.genes_", mark[i], "_" , genelist.name, "_", name1,".txt"), row.names=FALSE, col.names=FALSE, quote=FALSE)
     write.table(c(paste0("# unique.all.genes_", mark[i], "_" ,genelist.name, "_", name2), sort(unique.to.s2)), file=paste0("unique.all.genes_", mark[i], "_" , genelist.name, "_", name2,".txt"), row.names=FALSE, col.names=FALSE, quote=FALSE)
     write.table(c(paste0("# unique.all.genes_", mark[i], "_" ,genelist.name, "_", name3), sort(unique.to.s3)), file=paste0("unique.all.genes_", mark[i], "_" , genelist.name, "_", name3,".txt"), row.names=FALSE, col.names=FALSE, quote=FALSE)
     write.table(c(paste0("# unique.all.genes_", mark[i], "_" ,genelist.name, "_", name4), sort(unique.to.s4)), file=paste0("unique.all.genes_", mark[i], "_" , genelist.name, "_", name4,".txt"), row.names=FALSE, col.names=FALSE, quote=FALSE)
     write.table(c(paste0("# unique.all.genes_", mark[i], "_" ,genelist.name, "_", name5), sort(unique.to.s5)), file=paste0("unique.all.genes_", mark[i], "_" , genelist.name, "_", name5,".txt"), row.names=FALSE, col.names=FALSE, quote=FALSE)
     write.table(c(paste0("# unique.all.genes_", mark[i], "_" ,genelist.name, "_", name6), sort(unique.to.s6)), file=paste0("unique.all.genes_", mark[i], "_" , genelist.name, "_", name6,".txt"), row.names=FALSE, col.names=FALSE, quote=FALSE)
     
     # writeout common meth.genes across all timepoints
     write.table(c(paste0("# common.meth.genes_CONTROL_", mark[i], "_" ,genelist.name), sort(m.con_genes)), file=paste0("common.meth.genes_CONTROL_", mark[i], "_" , genelist.name, ".txt"), row.names=FALSE, col.names=FALSE, quote=FALSE)
     write.table(c(paste0("# common.meth.genes_CSC_", mark[i], "_" ,genelist.name), sort(m.csc_genes)), file=paste0("common.meth.genes_CSC_", mark[i], "_" , genelist.name, ".txt"), row.names=FALSE, col.names=FALSE, quote=FALSE)
     
     # writeout unique meth.genes at each timepoint
     write.table(c(paste0("# unique.meth.genes_", mark[i], "_" ,genelist.name, "_", name1), sort(unique.meth.genes.to.1)), file=paste0("unique.meth.genes_", mark[i], "_" , genelist.name, "_", name1,".txt"), row.names=FALSE, col.names=FALSE, quote=FALSE)
     write.table(c(paste0("# unique.meth.genes_", mark[i], "_" ,genelist.name, "_", name2), sort(unique.meth.genes.to.2)), file=paste0("unique.meth.genes_", mark[i], "_" , genelist.name, "_", name2,".txt"), row.names=FALSE, col.names=FALSE, quote=FALSE)
     write.table(c(paste0("# unique.meth.genes_", mark[i], "_" ,genelist.name, "_", name3), sort(unique.meth.genes.to.3)), file=paste0("unique.meth.genes_", mark[i], "_" , genelist.name, "_", name3,".txt"), row.names=FALSE, col.names=FALSE, quote=FALSE)
     write.table(c(paste0("# unique.meth.genes_", mark[i], "_" ,genelist.name, "_", name4), sort(unique.meth.genes.to.4)), file=paste0("unique.meth.genes_", mark[i], "_" , genelist.name, "_", name4,".txt"), row.names=FALSE, col.names=FALSE, quote=FALSE)
     write.table(c(paste0("# unique.meth.genes_", mark[i], "_" ,genelist.name, "_", name5), sort(unique.meth.genes.to.5)), file=paste0("unique.meth.genes_", mark[i], "_" , genelist.name, "_", name5,".txt"), row.names=FALSE, col.names=FALSE, quote=FALSE)
     write.table(c(paste0("# unique.meth.genes_", mark[i], "_" ,genelist.name, "_", name6), sort(unique.meth.genes.to.6)), file=paste0("unique.meth.genes_", mark[i], "_" , genelist.name, "_", name6,".txt"), row.names=FALSE, col.names=FALSE, quote=FALSE)
     
     ###########################################
     # PLOT 8
     # BOTH CONTROL & CSC - unique&common methylated.genes
     ###########################################
     ylim_custom1 <- c(0, max(cbind(n1, n2, n3, n4, n5, n6)))
     jpeg(paste0("unique&common.meth.genes_",mark[i],"_",genelist.name,".jpeg"), height = 900, width = 1600, quality=100)
     tmpmat <- rbind(c(0,1,2,3,0,4,5,0,6,7,8,0))
     layout(tmpmat, widths=c(1,3,3.5,3.5,3,3.5,3.5,3,3.5,3.5,11,2))
     par(mar=c(3,2,5,2))
     # stack1 - C10D, CSC10D
     #1
     barplot(t(c(1,2)), main="", col=NA, border="NA", axes=FALSE, names.arg=rep('',2), xpd=TRUE)
     mtext(text='# of genes',side=2,line=1, cex=2)
     #2
     # stack3 - C10D, CSC10D
     b <- barplot(t(cbind(c(unique.m1), c(m.con))), col=c("lightskyblue", "darkgoldenrod1") , names.arg=c(name1), cex.names=2.5, cex.axis=2.5, cex.lab=2.5, ylim=ylim_custom1, beside=FALSE) #ylab="# of genes", 
     text(b, 1, percentage.m1, cex=2.5, pos=3)
     b <- barplot(t(cbind(c(unique.m2), c(m.csc))), col=c("royalblue", "darkgoldenrod4") , names.arg=c(name2), cex.names=2.5, ylim=ylim_custom1, yaxt='n', beside=FALSE)
     text(b, 1, percentage.m2, cex=2.5, pos=3)
     # stack2 - C3M, CSC3M
     b <- barplot(t(cbind(c(unique.m3), c(m.con))), col=c("lightskyblue", "darkgoldenrod1") , names.arg=c(name3), cex.names=2.5, ylim=ylim_custom1, yaxt='n', beside=FALSE)
     text(b, 1, percentage.m3, cex=2.5, pos=3)
     b <- barplot(t(cbind(c(unique.m4), c(m.csc))), col=c("royalblue", "darkgoldenrod4") , names.arg=c(name4), cex.names=2.5, ylim=ylim_custom1, yaxt='n', beside=FALSE)
     text(b, 1, percentage.m4, cex=2.5, pos=3)
     # stack3 - C10M, CSC10M
     b <- barplot(t(cbind(c(unique.m5), c(m.con))), col=c("lightskyblue", "darkgoldenrod1") , names.arg=c(name5), cex.names=2.5, ylim=ylim_custom1, yaxt='n', beside=FALSE)
     text(b, 1, percentage.m5, cex=2.5, pos=3)
     b <- barplot(t(cbind(c(unique.m6), c(m.csc))), col=c("royalblue", "darkgoldenrod4") , names.arg=c(name6), cex.names=2.5, ylim=ylim_custom1, yaxt='n', beside=FALSE)
     text(b, 1, percentage.m6, cex=2.5, pos=3)
     # dummy plot
     barplot(t(c(1,2)), main="", col=NA, border="NA", axes=FALSE, names.arg=rep('',2), xpd=TRUE)
     legend("center", fill=c("lightskyblue", "royalblue", "darkgoldenrod1", "white", "darkgoldenrod4", "white"), legend=c("unique CON meth.genes", "unique CSC meth.genes", "common meth.genes", "across all CON timepoints","common meth.genes", "across all CSC timepoints"), inset=c(-1.6,0), border="white", box.col="white", cex=2.7)
     title(paste0(mark[i], " - ",genelist.name,"\nunique&common meth.genes"), cex.main=3, outer=FALSE, line=-2)
     dev.off()
     
     ###########################################
     # PLOT 9
     # BOTH CONTROL & CSC - unique&common all.genes
     ###########################################
     ylim_custom1 <- c(0, max(cbind(a1, a2, a3, a4, a5, a6)))
     ylim_custom1
     jpeg(paste0("unique&common.all.genes_",mark[i],"_",genelist.name,".jpeg"), height = 900, width = 1600, quality=100)
     tmpmat <- rbind(c(0,1,2,3,0,4,5,0,6,7,8,0))
     layout(tmpmat, widths=c(1,3,3.5,3.5,3,3.5,3.5,3,3.5,3.5,11,2))
     par(mar=c(3,2,5,2))
     # stack1 - C10D, CSC10D
     #1
     barplot(t(c(1,2)), main="", col=NA, border="NA", axes=FALSE, names.arg=rep('',2), xpd=TRUE)
     mtext(text='# of genes',side=2,line=1, cex=2)
     #2
     # stack3 - C10D, CSC10D
     b <- barplot(t(cbind(c(unique1), c(all.con))), col=c("green2", "orange") , names.arg=c(name1), cex.names=2.5, cex.axis=2.5, cex.lab=2.5, ylim=ylim_custom1, beside=FALSE) #ylab="# of genes", 
     text(b, 1, percentage.1, cex=2.5, pos=3)
     b <- barplot(t(cbind(c(unique2), c(all.csc))), col=c("green4", "red") , names.arg=c(name2), cex.names=2.5, ylim=ylim_custom1, yaxt='n', beside=FALSE)
     text(b, 1, percentage.2, cex=2.5, pos=3)
     # stack2 - C3M, CSC3M
     b <- barplot(t(cbind(c(unique3), c(all.con))), col=c("green2", "orange") , names.arg=c(name3), cex.names=2.5, ylim=ylim_custom1, yaxt='n', beside=FALSE)
     text(b, 1, percentage.3, cex=2.5, pos=3)
     b <- barplot(t(cbind(c(unique4), c(all.csc))), col=c("green4", "red") , names.arg=c(name4), cex.names=2.5, ylim=ylim_custom1, yaxt='n', beside=FALSE)
     text(b, 1, percentage.4, cex=2.5, pos=3)
     # stack3 - C10M, CSC10M
     b <- barplot(t(cbind(c(unique5), c(all.con))), col=c("green2", "orange") , names.arg=c(name5), cex.names=2.5, ylim=ylim_custom1, yaxt='n', beside=FALSE)
     text(b, 1, percentage.5, cex=2.5, pos=3)
     b <- barplot(t(cbind(c(unique6), c(all.csc))), col=c("green4", "red") , names.arg=c(name6), cex.names=2.5, ylim=ylim_custom1, yaxt='n', beside=FALSE)
     text(b, 1, percentage.6, cex=2.5, pos=3)
     # dummy plot
     barplot(t(c(1,2)), main="", col=NA, border="NA", axes=FALSE, names.arg=rep('',2), xpd=TRUE)
     legend("center", fill=c("green2", "green4", "orange", "white", "red", "white"), legend=c("unique CON genes", "unique CSC genes", "common genes across", "all CON timepoints","common genes across", "all CSC timepoints"), inset=c(-1.6,0), border="white", box.col="white", cex=2.7)
     title(paste0(mark[i], " - ",genelist.name,"\nunique&common meth.genes"), cex.main=3, outer=FALSE, line=-2)
     dev.off()
     
     Sys.sleep(3)
}
