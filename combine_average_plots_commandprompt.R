##########
# Purpose:
##########
# To plot all of the average TSS genecoverage values for all samples on a single plot

################
# Prerequisites:
################
# 1.) run create_ChIPseq_rowsidecolor_annotation.R
# 2.) run chipseq_script.R
# 3.) run consolidate_files.sh # this will consolidate all average plots + values.txt into separate folders
# system.dir="Z:/users/shwang26/"

args <- commandArgs(trailingOnly=TRUE)
whichgenelist <- as.character(args[1])
bin.size <- as.numeric(args[2]) #10bp or 200bp
# whichgenelist<-"rep1.agerelated.10M.new.noncpg"
# whichgenelist <- "random.cpg"
bin.size=200
print(paste0("genelist: ", whichgenelist))
print(paste0("bin: ", bin.size))

system.dir="/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/"
source(paste0(system.dir, "Michelle/Rscripts/ChIP-SeqLibraryOfFunctions_newFunctions.sort.R"))

marks <- c("H3K4", "H3K27", "EZH2", "DNMT1", "Input")
timepoints <- c("C10D", "CSC10D", "C3M", "CSC3M", "C10M", "CSC10M")

# Generate Input-normalized values for each group at every timepoint
for(i in marks){
     print("################################")
     print(i)
     print("################################")
     input.dir <- paste0(system.dir, "/Michelle/BED_files/Coverage_TSS_200bp_bin/normalizedBED_200bp_bin/outputdir/", whichgenelist, "/Input/")
     marks.dir <- paste0(system.dir, "/Michelle/BED_files/Coverage_TSS_200bp_bin/normalizedBED_200bp_bin/outputdir/", whichgenelist, "/", i, "/")
     for(j in timepoints){
          # read in Input file for that particular timepoint
          input.file <- dir(input.dir)[grep(paste0("^",j ,".*average_values.txt"), dir(input.dir))]
          path.to.input <- paste0(input.dir, input.file)
          input.file <- as.matrix(read.table(path.to.input, header=TRUE))
          # read in sample file for that particular timepoint
          marks.file <- dir(marks.dir)[grep(paste0("^",j,".*average_values.txt"), dir(marks.dir))]
          path.to.mark <- paste0(marks.dir, marks.file)
          path.to.mark
          marks.file <- as.matrix(read.table(path.to.mark, header=TRUE))
          # normalize to input
          ratio.file <- marks.file/input.file
          print(head(ratio.file))
          write.table(ratio.file, file=paste0(system.dir, "/Michelle/BED_files/Coverage_TSS_200bp_bin/normalizedBED_200bp_bin/outputdir/", whichgenelist, "/", i, "/", j, "_",i, "_ratio_values.txt"), row.names=FALSE, col.names=TRUE, quote=FALSE)
     }
}

# Combine average values into a single plot
# Combine Input-normalized values into a single plot
# sample.type=c("DNMT1", "EZH2", "H3", "H3K4", "H3K27", "Input")
for(i in marks){
     output.dir <- paste0(system.dir, "Michelle/BED_files/Coverage_TSS_",bin.size,"bp_bin/normalizedBED_",bin.size,"bp_bin/outputdir/", whichgenelist, "/", i, "/")
     print(output.dir)
#      generate_combined_avg_plots(sample.name=i, name.of.genelist=whichgenelist, bin=bin.size, outputDirectory=output.dir)
     generate_combined_avg_plots_both(sample.name=i, name.of.genelist=whichgenelist, bin=bin.size, outputDirectory=output.dir)
     print(i)
}
