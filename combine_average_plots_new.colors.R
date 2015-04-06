#/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/
#/amber2/scratch/baylin/shwang/
#/Volumes/onc-analysis$/users/stephenhwang/
# system.dir="Z:/users/shwang26/"

system.dir="/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/"
source(paste0(system.dir, "Michelle/Rscripts/ChIP-SeqLibraryOfFunctions_original_newFunctions.R"))

bin.size=10
whichgenelist="stable.10M.new" # CHANGE DEPENDING ON GENELIST USED

sample.type=c("DNMT1", "EZH2", "H3", "H3K4", "H3K27", "INPUT")
for(i in sample.type){
     output.dir <- paste0(system.dir, "Michelle/BED_files/Coverage_TSS_",bin.size,"bp_bin/normalizedBED_",bin.size,"bp_bin/outputdir/", whichgenelist, "/", i, "/")
     generate_combined_avg_plots(sample.name=i, name.of.genelist=whichgenelist, bin=bin.size, outputDirectory=output.dir)
     print(i)
}
dir(directory)
