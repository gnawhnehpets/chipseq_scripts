# Read in chipseq functions and genelist annotations
system.dir='/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/'
# system.dir='Z:/users/shwang26/'
# system.dir='/amber2/scratch/baylin/shwang/'

source(paste0(system.dir, "Michelle/Rscripts/ChIP-SeqLibraryOfFunctions_original_newFunctions.R"))
options(bitmapType='cairo') 
bin.size=200 #10bp or 200bp

# genelists
# genelist.name <- "stable.10M"
# genelist.name <- "stable.10M.new"
# genelist.name <- "intermediate.10M"
# genelist.name <- "intermediate.10M.new"
genelist.name <- "highly.expressed"
# genelist.name <- "all"
# genelist.name <- ""
whichdata <- "normalized"

if(genelist.name=="all"){
     print("all genes requested")
     plot_results_dir = paste0(system.dir, "Michelle/BED_files/Coverage_TSS_",bin.size,"bp_bin/",whichdata,"BED_",bin.size,"bp_bin/outputdir/",genelist.name,"/")
     plot_results_dir
     version="hg19-UCSC"
}else{
     version="hg19"
     if(genelist.name=="stable.10M"){
          print("stable.10M genelist selected")
          plot_results_dir = paste0(system.dir, "Michelle/BED_files/Coverage_TSS_",bin.size,"bp_bin/",whichdata,"BED_",bin.size,"bp_bin/outputdir/stable.10M/")
          genes <- read.table(paste0(system.dir, "Michelle/MethylationData/methylatedGeneLists/new/age_stable_methylated_genes_at_10months_new_annotation.txt"))
     }
     if(genelist.name=="intermediate.10M"){
          print("intermediate.10M genelist selected")
          plot_results_dir = paste0(system.dir, "Michelle/BED_files/Coverage_TSS_",bin.size,"bp_bin/",whichdata,"BED_",bin.size,"bp_bin/outputdir/intermediate.10M/")
          genes <- read.table(paste0(system.dir, "Michelle/MethylationData/methylatedGeneLists/new/age_intermediate_methylated_genes_at_10months_new_annotation.txt"))
     }     
     if(genelist.name=="stable.10M.new"){
          print("stable.10M.new genelist selected")
          plot_results_dir = paste0(system.dir, "Michelle/BED_files/Coverage_TSS_",bin.size,"bp_bin/",whichdata,"BED_",bin.size,"bp_bin/outputdir/stable.10M.new/")
          genes <- read.table(paste0(system.dir, "Michelle/MethylationData/methylatedGeneLists/new/age_stable_methylated_genes_at_10months_new_cutoff.txt"))
     }
     if(genelist.name=="intermediate.10M.new"){
          print("intermediate.10M.new genelist selected")
          plot_results_dir = paste0(system.dir, "Michelle/BED_files/Coverage_TSS_",bin.size,"bp_bin/",whichdata,"BED_",bin.size,"bp_bin/outputdir/intermediate.10M.new/")
          genes <- read.table(paste0(system.dir, "Michelle/MethylationData/methylatedGeneLists/new/age_intermediate_methylated_genes_at_10months_new_cutoff.txt"))
     }
     if(genelist.name=="highly.expressed"){
          print("highly.expressed genelist selected")
          plot_results_dir = paste0(system.dir, "Michelle/BED_files/Coverage_TSS_",bin.size,"bp_bin/",whichdata,"BED_",bin.size,"bp_bin/outputdir/",genelist.name,"/")
          plot_results_dir
          genes <- read.table(paste0(system.dir, "Michelle/MethylationData/otherGenelists/list of high expression genes.txt"))
     }
     h3k4.maxbreak=1500
     h3k27.maxbreak=250
     dnmt.maxbreak=50
     ezh2.maxbreak=120
     inp.maxbreak=50
     h3.maxbreak=50
     # Read in genelist
     genes <- as.list(sort(unique(unlist(strsplit(unlist(as.matrix(genes)), ";")))))
     genes <- list(unlist(genes))
     names(genes) <- genelist.name
}

#####################################################################################################################
# Location of normalized BED files
# Create directory if !exist
dir.create(plot_results_dir)
coverage_files_dir <- paste0(system.dir, "Michelle/BED_files/Coverage_TSS_",bin.size,"bp_bin/",whichdata,"BED_",bin.size,"bp_bin/")
coverage_files_dir

# # Read in genelist
# genes <- as.list(sort(unique(unlist(strsplit(unlist(as.matrix(genes)), ";")))))
# genes <- list(unlist(genes))
# names(genes) <- genelist.name

############################################################
#Get info about genes of interest from biomaRt
############################################################
# Add this into the 'if' statements above ??
# annotation for just genes of interest
if(version=="hg19"){
     goi <- fun.genelist.info(goi.list=genes, genelist.names=genelist.name, goi.list.type="hgnc_symbol", version=version)
}
# annotation for all genes
if(version=="hg19-UCSC"){
     hg19ucsc <- fun.genelist.info_allGenes_biomaRt(version="hg19-UCSC", 
                                                    hg19UCSCGeneAnnotations="hg19UCSCGeneAnnotaions.txt", 
                                                    hg19UCSCGeneAnnotationsPath=paste0(system.dir, "/Michelle/Required_Files/Annotations/hg19Data/"))
     dat <- cbind(hg19ucsc$Gene, hg19ucsc$Chromosome, hg19ucsc$TSS, hg19ucsc$TES, hg19ucsc$Strand, hg19ucsc$Gene)
     colnames(dat) <- c("Gene", "Chromosome", "TSS", "TES", "Strand", "ExternalGeneID")
     datlist <- list()
     datlist[[1]] <- as.data.frame(dat)
     goi <- datlist
}
# if(version==NULL){
#      version="hg19"
#      print ("warning: need to designate a gene name to the 'genelist.name' object")
#      datlist <- list()
#      datlist[[1]] <- as.data.frame(genelist.name)
#      goi <- datlist
# }

#Retain one entry per gene (hgnc_symbol)
goi <- lapply(c(1:length(goi)), FUN=function(i){
     x <- goi[[i]]
     return(x[!(duplicated(x[,"Gene"])), ])
})
names(goi) <- genelist.name
head(goi[[1]])
dim(goi[[1]]) #

############################################################
# Make plots
## Make plots using coverage data generated form filtered bam files
############################################################
#The coverage data used here is 10bp binned data across the TSS
## The argument RegAroundTSS is set to 10000, which means that coverage will be plotted for ±5000bp from the TSS
# # For troubleshooting
# genelist.info=goi
# coverage_files=coverage_files
# input_coverage_file=input_coverage_file
# coverage_files.dir=coverage_files_dir
# RegAroundTSS=10000
# bin=bin.size
# chr.prefix.chromosome=T
# whichgenelist=genelist.name
# plot.Directory=paste0(system.dir, "Michelle/BED_files/Coverage_TSS_200bp_bin/normalizedBED_200bp_bin/outputdir/", genelist.name)


# C10D
dir(coverage_files_dir)
coverage_files <- dir(coverage_files_dir)[grep("C10D", dir(coverage_files_dir))[1:6]]
coverage_files
input_coverage_file = coverage_files[6]
input_coverage_file
fun.average_heat.plots(genelist.info=goi, coverage_files=coverage_files, input_coverage_file=input_coverage_file, coverage_files.dir=coverage_files_dir, RegAroundTSS=10000, bin=bin.size, chr.prefix.chromosome=T, plot.Directory=paste0(plot_results_dir, "C10D"), version=version, whichgenelist=genelist.name, h3k4.max=h3k4.maxbreak, h3k27.max=h3k27.maxbreak, dnmt.max=dnmt.maxbreak, ezh2.max=ezh2.maxbreak, inp.max=inp.maxbreak, h3.max=h3.maxbreak)

# CSC10D
coverage_files = dir(coverage_files_dir)[grep("C10D", dir(coverage_files_dir))[7:12]]
input_coverage_file = coverage_files[6]
# coverage_files <- coverage_files[-6]
input_coverage_file
coverage_files
fun.average_heat.plots(genelist.info=goi, coverage_files=coverage_files, input_coverage_file=input_coverage_file, coverage_files.dir=coverage_files_dir, RegAroundTSS=10000, bin=bin.size, chr.prefix.chromosome=T, plot.Directory=paste0(plot_results_dir, "CSC10D"), version=version, whichgenelist=genelist.name, h3k4.max=h3k4.maxbreak, h3k27.max=h3k27.maxbreak, dnmt.max=dnmt.maxbreak, ezh2.max=ezh2.maxbreak, inp.max=inp.maxbreak, h3.max=h3.maxbreak)

# C3M
coverage_files = dir(coverage_files_dir)[grep("C3M", dir(coverage_files_dir))[1:6]]
input_coverage_file = coverage_files[6]
# coverage_files <- coverage_files[-6]
input_coverage_file
coverage_files
fun.average_heat.plots(genelist.info=goi, coverage_files=coverage_files, input_coverage_file=input_coverage_file, coverage_files.dir=coverage_files_dir, RegAroundTSS=10000, bin=bin.size, chr.prefix.chromosome=T, plot.Directory=paste0(plot_results_dir, "C3M"), version=version, whichgenelist=genelist.name, h3k4.max=h3k4.maxbreak, h3k27.max=h3k27.maxbreak, dnmt.max=dnmt.maxbreak, ezh2.max=ezh2.maxbreak, inp.max=inp.maxbreak, h3.max=h3.maxbreak)

# CSC3M
coverage_files = dir(coverage_files_dir)[grep("C3M", dir(coverage_files_dir))[7:12]] #CSC3M
input_coverage_file = coverage_files[6]
# coverage_files <- coverage_files[-6]
input_coverage_file
coverage_files
fun.average_heat.plots(genelist.info=goi, coverage_files=coverage_files, input_coverage_file=input_coverage_file, coverage_files.dir=coverage_files_dir, RegAroundTSS=10000, bin=bin.size, chr.prefix.chromosome=T, plot.Directory=paste0(plot_results_dir, "CSC3M"), version=version, whichgenelist=genelist.name, h3k4.max=h3k4.maxbreak, h3k27.max=h3k27.maxbreak, dnmt.max=dnmt.maxbreak, ezh2.max=ezh2.maxbreak, inp.max=inp.maxbreak, h3.max=h3.maxbreak)

# C10M
coverage_files = dir(coverage_files_dir)[grep("C10M", dir(coverage_files_dir))[1:6]] 
input_coverage_file = coverage_files[6]
# coverage_files <- coverage_files[-6]
input_coverage_file
coverage_files
fun.average_heat.plots(genelist.info=goi, coverage_files=coverage_files, input_coverage_file=input_coverage_file, coverage_files.dir=coverage_files_dir, RegAroundTSS=10000, bin=bin.size, chr.prefix.chromosome=T, plot.Directory=paste0(plot_results_dir, "C10M"), version=version, whichgenelist=genelist.name, h3k4.max=h3k4.maxbreak, h3k27.max=h3k27.maxbreak, dnmt.max=dnmt.maxbreak, ezh2.max=ezh2.maxbreak, inp.max=inp.maxbreak, h3.max=h3.maxbreak)

# CSC10M
coverage_files = dir(coverage_files_dir)[grep("C10M", dir(coverage_files_dir))[7:12]]
input_coverage_file = coverage_files[6]
# coverage_files <- coverage_files[-6]
input_coverage_file
coverage_files
fun.average_heat.plots(genelist.info=goi, coverage_files=coverage_files, input_coverage_file=input_coverage_file, coverage_files.dir=coverage_files_dir, RegAroundTSS=10000, bin=bin.size, chr.prefix.chromosome=T, plot.Directory=paste0(plot_results_dir, "CSC10M"), version=version, whichgenelist=genelist.name, h3k4.max=h3k4.maxbreak, h3k27.max=h3k27.maxbreak, dnmt.max=dnmt.maxbreak, ezh2.max=ezh2.maxbreak, inp.max=inp.maxbreak, h3.max=h3.maxbreak)
