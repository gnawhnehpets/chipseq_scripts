#/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/
#/amber2/scratch/baylin/shwang/
#/Volumes/onc-analysis$/users/stephenhwang/
# dir="Z:/users/shwang26/"
dir="/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/"
source(paste0(dir, "Michelle/Rscripts/ChIP-SeqLibraryOfFunctions_original_newFunctions.R"))
bin.size=200
sample.type=c("DNMT1", "EZH2", "H3", "H3K4", "H3K27", "INPUT")
whichgenelist="stable.10M.new" # CHANGE DEPENDING ON GENELIST USED
for(i in sample.type){
     output.dir <- paste0(dir, "Michelle/BED_files/Coverage_TSS_",bin.size,"bp_bin/unnormalizedBED_",bin.size,"bp_bin/outputdir/", whichgenelist, "/", i, "/")
     generate_combined_avg_plots(sample.name=i, name.of.genelist=whichgenelist, bin=bin.size, outputDirectory=output.dir)
     print(i)
}
dir(directory)
# generate_combined_avg_plots <- function(sample.name, name.of.genelist){
#      directory = paste0(dir, "michelle/BED_files/Coverage_TSS_200bp_bin/normalizedBED_200bp_bin/outputdir/", name.of.genelist, "/", sample.name, "/")
#      directory
#      files <- list.files(directory, pattern=".*R1.*.txt")
#      files <- files[c(1,3,2,4,6,5)]
#      files
#      # read in average values into list
#      datalist <- list()
#      d<- data.matrix(read.delim(paste0(directory,files[1]), header=TRUE))
#      class(d)
#      for(i in 1:length(files)){
#           tmp <- data.matrix(read.table(paste0(directory,files[i]), header=TRUE))
#           datalist[[i]]  <- tmp
#           colnames(datalist[[i]]) <- gsub("_R1.*", "", colnames(datalist[[i]]))
#      }
#      class(datalist[[5]])
#      datalist
#      int <- as.numeric(unlist(read.delim(paste0(directory, "int.txt"), sep="\t", header=FALSE)))
#      length(int)
#      int
#      
#      # Average plot
#      filename=paste0(directory, "combined_avg_plots_", name.of.genelist, "_", sample.name, ".jpeg")
#      filename
#      jpeg(filename, height=600, width=900, quality=100)
#      tmpmat <- rbind(c(0,1,1,1,2))
#      layout(tmpmat, widths=c(.1,.1))
#      #PLOT1
#      par(mar=c(5,7,5,0))
#      plot(0~0, type="n", xlim=range(int), ylim=c(0,range(datalist, na.rm=T)[2]), xlab="Dist. from TSS", ylab="Seq Counts", cex.lab=1.8, cex.axis=1.8, cex.main=2, main=paste0("Combined Average Plot for ", sample.name, " timepoints (", name.of.genelist, ")"))
#      plot.col <- vector()
#      plot.lty <- vector()
#      sample.name <- vector()
#      for(u in 1:length(datalist)){
#           sample.name <- c(sample.name, colnames(datalist[[u]]))
#           # Start loop average plot
#           plot.col[u] <- u
#           if(!any(datalist[[u]] %in% NA)){
#                lines(int, datalist[[u]], lty=1, lwd=4, col=u)
#                plot.lty[u] <- 1
#           } else{lines(int, rep(0, num.bins), lty=3, lwd=4, col=u)
#                  text(x=int[1], y=0, labels=c("coverage data for atleast two genes not present"), adj=0, cex=0.8)
#                  plot.lty[u] <- 3}
#      } # end plotting average plot
#      #PLOT2 (dummy plot)
#      barplot(t(c(1,2)), main="", col=NA, border="NA", axes=FALSE, names.arg=rep('',2), xpd=TRUE)
#      legend("left", fill=c(1:length(datalist)), legend=sample.name, border="white", box.col="white", cex=2)
#      dev.off()
# }
