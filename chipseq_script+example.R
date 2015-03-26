#/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/
#/amber2/scratch/baylin/shwang/
#/Volumes/onc-analysis$/users/stephenhwang/

# Notes
# Untreated vs smoke-treated cells
# Summary: create a model in lab for chronical smoking in lung cell lining
#         tumor from this lining have genetic and epigenetic changes but do not understand well
# initial analysis - genome wide
# small bar - list of genes that are methylated; how many of those genes that are methylated are
#         k4 marked at the time point
#         basically
# in culture, over 3months, across the genome, the genes take an increase in k27
# from TSS, +/- 2500 bp
#         that wide for k27
#         k4 - wouldn't have that big of a window
#              nucleosome starts around 2kbp from TSS
# what changes here is bivalency; how did it change?
#         genes with both marks is increasing
#         

# Read in chipseq functions and genelist annotations
# source("/amber2/scratch/baylin/shwang/Michelle/ChIP-SeqLibraryOfFunctions/ChIP-SeqLibraryOfFunctions_originalcompletecases2.R")
# source("/amber2/scratch/baylin/shwang/Michelle/ChIP-SeqLibraryOfFunctions/ChIP-SeqLibraryOfFunctions_original.R")
source("/amber2/scratch/baylin/shwang/Michelle/Rscripts/ChIP-SeqLibraryOfFunctions_original_newFunctions.R")
options(bitmapType='cairo') 
bin.size=10 #10bp or 200bp
# coverage_files_dir = "/amber2/scratch/baylin/Hari/Michelle/BED_Files/new/Coverage_200bp_bin"
#1
genelist.name <- "stable.10M"
#plot_results_dir = "/amber2/scratch/baylin/shwang/Michelle/BED_Files/plots/coverage_200bp_bin/10M_intermediate"
#genes <- read.table("/amber2/scratch/baylin/Hari/Michelle/MethylationData/methylatedGeneLists/age_intermediate_methylated_genes_at_10months_new_annotation.txt", header=FALSE)
            
plot_results_dir = paste0("/amber2/scratch/baylin/shwang/Michelle/BED_files/Coverage_TSS_",bin.size,"bp_bin/normalizedBED_",bin.size,"bp_bin/outputdir/stable.10M/")
genes <- read.table("/amber2/scratch/baylin/shwang/Michelle/MethylationData/methylatedGeneLists/new/age_stable_methylated_genes_at_10months_new_annotation.txt")

#2
genelist.name <- "intermediate.10M"
#plot_results_dir = "/amber2/scratch/baylin/shwang/Michelle/BED_Files/plots/coverage_200bp_bin/10M_stable"
#genes <- read.table("/amber2/scratch/baylin/Hari/Michelle/MethylationData/methylatedGeneLists/age_stable_methylated_genes_at_10months_new_annotation.txt", header=FALSE)# 
plot_results_dir = "/amber2/scratch/baylin/shwang/michelle/BED_files/Coverage_TSS_200bp_bin/normalizedBED_200bp_bin/outputdir/intermediate.10M/"
genes <- read.table("/amber2/scratch/baylin/shwang/michelle/MethylationData/methylatedGeneLists/new/age_intermediate_methylated_genes_at_10months_new_annotation.txt")

#3
genelist.name <- "stable.10M.new"
#plot_results_dir = "/amber2/scratch/baylin/shwang/Michelle/BED_Files/plots/coverage_200bp_bin/15M_intermediate"
#genes <- read.table("/amber2/scratch/baylin/Hari/Michelle/MethylationData/methylatedGeneLists/age_intermediate_methylated_genes_at_15months_new_annotation.txt", header=FALSE)
plot_results_dir = "/amber2/scratch/baylin/shwang/michelle/BED_files/Coverage_TSS_200bp_bin/normalizedBED_200bp_bin/outputdir/stable.10M.new/"
genes <- read.table("/amber2/scratch/baylin/shwang/michelle/MethylationData/methylatedGeneLists/new/age_stable_methylated_genes_at_10months_new_cutoff.txt")

#4
genelist.name <- "intermediate.10M.new"
# plot_results_dir = "/amber2/scratch/baylin/shwang/Michelle/BED_Files/plots/coverage_200bp_bin/15M_stable"
#genes <- read.table("/amber2/scratch/baylin/Hari/Michelle/MethylationData/methylatedGeneLists/age_stable_methylated_genes_at_15months_new_annotation.txt", header=FALSE)
plot_results_dir = "/amber2/scratch/baylin/shwang/michelle/BED_files/Coverage_TSS_200bp_bin/normalizedBED_200bp_bin/outputdir/intermediate.10M.new/"
genes <- read.table("/amber2/scratch/baylin/shwang/michelle/MethylationData/methylatedGeneLists/new/age_intermediate_methylated_genes_at_10months_new_cutoff.txt")


#####################################################################################################################
# Location of normalized BED files
#coverage_files_dir = "/amber2/scratch/baylin/shwang/Michelle/normalizedBED"

coverage_files_dir <- paste0("/amber2/scratch/baylin/shwang/Michelle/BED_files/Coverage_TSS_",bin.size,"bp_bin/normalizedBED_",bin.size,"bp_bin/")
coverage_files_dir
# Read in genelist
genes <- as.list(sort(unique(unlist(strsplit(unlist(as.matrix(genes)), ";")))))
genes <- list(unlist(genes))
names(genes) <- genelist.name

version = "hg19"
goi <- fun.genelist.info(goi.list=genes, genelist.names=genelist.name, goi.list.type="hgnc_symbol", version=version)

############################################################
#Get info about genes of interest from biomaRt
############################################################
#Retain one entry per gene (hgnc_symbol)
goi <- lapply(c(1:length(goi)), FUN=function(i){
     x <- goi[[i]]
     return(x[!(duplicated(x[,"Gene"])), ])
})
names(goi) <- names(genes)
goi

# EXAMPLE
# options(bitmapType='cairo') 
# coverage_files.dir = "/amber2/scratch/baylin/Hari/JHUSB01002/JHUSB01002_000_analysis/BED_Files/new/Coverage_TSS_10bp_bin"
# coverage_files = c("JHUSB01002_002_K4_SS02.sorted.bam_FILTERED.bed_Coverage_TSS_10bp.bed", "JHUSB01002_003_K27_SS03.sorted.bam_FILTERED.bed_Coverage_TSS_10bp.bed", "JHUSB01002_004_2o_SS04.sorted.bam_FILTERED.bed_Coverage_TSS_10bp.bed", "JHUSB01002_005_H2AZ_SS05.sorted.bam_FILTERED.bed_Coverage_TSS_10bp.bed", "JHUSB01002_001_INP_SS01.sorted.bam_FILTERED.bed_Coverage_TSS_10bp.bed")
# input_coverage_file = "JHUSB01002_001_INP_SS01.sorted.bam_FILTERED.bed_Coverage_TSS_10bp.bed"
# fun.average_heat.plots(genelist.info=goi, 
#                        coverage_files=coverage_files, 
#                        input_coverage_file=input_coverage_file, 
#                        coverage_files.dir=coverage_files.dir, 
#                        RegAroundTSS=10000, 
#                        bin=10, 
#                        chr.prefix.chromosome=T, 
#                        plot.Directory="/amber2/scratch/baylin/shwang/Michelle/BED_Files/example/outputdir")
# 

############################################################
# Make plots
## Make plots using coverage data generated form filtered bam files
############################################################
#The coverage data used here is 10bp binned data across the TSS
## The argument RegAroundTSS is set to 10000, which means that coverage will be plotted for ±5000bp from the TSS
# coverage_files.dir = "/amber2/scratch/baylin/Hari/JHUSB01002/JHUSB01002_000_analysis/BED_Files/new/Coverage_TSS_10bp_bin"

# # For troubleshooting
# genelist.info=goi
# coverage_files=coverage_files
# input_coverage_file=input_coverage_file
# coverage_files.dir=coverage_files_dir
# RegAroundTSS=10000
# bin=bin.size
# chr.prefix.chromosome=T
# whichgenelist=genelist.name
# # plot.Directory="/amber2/scratch/baylin/shwang/michelle/Analysis/chipseq_plot/updated_code_3.19.2015/"
# plot.Directory=paste0("/amber2/scratch/baylin/shwang/michelle/BED_files/Coverage_TSS_200bp_bin/normalizedBED_200bp_bin/outputdir/", genelist.name)
# # plot.Directory="/amber2/scratch/baylin/shwang/Michelle/BED_Files/plots"


# C10D
coverage_files <- dir(coverage_files_dir)[grep("C10D", dir(coverage_files_dir))[1:6]]
coverage_files
input_coverage_file = coverage_files[6]
input_coverage_file
fun.average_heat.plots(genelist.info=goi, coverage_files=coverage_files, input_coverage_file=input_coverage_file, coverage_files.dir=coverage_files_dir, RegAroundTSS=10000, bin=bin.size, chr.prefix.chromosome=T, plot.Directory=paste0(plot_results_dir, "C10D"), version=version, whichgenelist=genelist.name)

# CSC10D
coverage_files = dir(coverage_files_dir)[grep("C10D", dir(coverage_files_dir))[7:12]]
input_coverage_file = coverage_files[6]
# coverage_files <- coverage_files[-6]
input_coverage_file
coverage_files
fun.average_heat.plots(genelist.info=goi, coverage_files=coverage_files, input_coverage_file=input_coverage_file, coverage_files.dir=coverage_files_dir, RegAroundTSS=10000, bin=bin.size, chr.prefix.chromosome=T, plot.Directory=paste0(plot_results_dir, "CSC10D"), version=version, whichgenelist=genelist.name)

# C3M
coverage_files = dir(coverage_files_dir)[grep("C3M", dir(coverage_files_dir))[1:6]]
input_coverage_file = coverage_files[6]
# coverage_files <- coverage_files[-6]
input_coverage_file
coverage_files
fun.average_heat.plots(genelist.info=goi, coverage_files=coverage_files, input_coverage_file=input_coverage_file, coverage_files.dir=coverage_files_dir, RegAroundTSS=10000, bin=bin.size, chr.prefix.chromosome=T, plot.Directory=paste0(plot_results_dir, "C3M"), version=version, whichgenelist=genelist.name)

# CSC3M
coverage_files = dir(coverage_files_dir)[grep("C3M", dir(coverage_files_dir))[7:12]] #CSC3M
input_coverage_file = coverage_files[6]
# coverage_files <- coverage_files[-6]
input_coverage_file
coverage_files
fun.average_heat.plots(genelist.info=goi, coverage_files=coverage_files, input_coverage_file=input_coverage_file, coverage_files.dir=coverage_files_dir, RegAroundTSS=10000, bin=bin.size, chr.prefix.chromosome=T, plot.Directory=paste0(plot_results_dir, "CSC3M"), version=version, whichgenelist=genelist.name)

# C10M
coverage_files = dir(coverage_files_dir)[grep("C10M", dir(coverage_files_dir))[1:6]] 
input_coverage_file = coverage_files[6]
# coverage_files <- coverage_files[-6]
input_coverage_file
coverage_files
fun.average_heat.plots(genelist.info=goi, coverage_files=coverage_files, input_coverage_file=input_coverage_file, coverage_files.dir=coverage_files_dir, RegAroundTSS=10000, bin=bin.size, chr.prefix.chromosome=T, plot.Directory=paste0(plot_results_dir, "C10M"), version=version, whichgenelist=genelist.name)

# CSC10M
coverage_files = dir(coverage_files_dir)[grep("C10M", dir(coverage_files_dir))[7:12]]
input_coverage_file = coverage_files[6]
# coverage_files <- coverage_files[-6]
input_coverage_file
coverage_files
fun.average_heat.plots(genelist.info=goi, coverage_files=coverage_files, input_coverage_file=input_coverage_file, coverage_files.dir=coverage_files_dir, RegAroundTSS=10000, bin=bin.size, chr.prefix.chromosome=T, plot.Directory=paste0(plot_results_dir, "CSC10M"), version=version, whichgenelist=genelist.name)

