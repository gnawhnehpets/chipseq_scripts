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
genelist.name <- as.character(args[1])
bin.size <- as.numeric(args[2]) #10bp or 200bp
print(paste0("genelist: ", genelist.name))
print(paste0("bin: ", bin.size))

system.dir="/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/"
source(paste0(system.dir, "Michelle/Rscripts/ChIP-SeqLibraryOfFunctions_original_newFunctions.R"))

bin.size=200
whichgenelist="stable.10M.newdata" # CHANGE DEPENDING ON GENELIST USED

sample.type=c("DNMT1", "EZH2", "H3", "H3K4", "H3K27", "INPUT")
for(i in sample.type){
     output.dir <- paste0(system.dir, "Michelle/BED_files/Coverage_TSS_",bin.size,"bp_bin/normalizedBED_",bin.size,"bp_bin/outputdir/", whichgenelist, "/", i, "/")
     generate_combined_avg_plots(sample.name=i, name.of.genelist=whichgenelist, bin=bin.size, outputDirectory=output.dir)
     print(i)
}
dir(directory)
