#/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/
#/amber2/scratch/baylin/shwang/
#/Volumes/onc-analysis$/users/stephenhwang/

# RAW READS: BARPLOT 
getwd()
# setwd("/Volumes/onc-analysis$/users/shwang26/michelle/BED_files/Coverage_TSS_200bp_bin/unnormalizedBED_200bp_bin/")
bin.size=10
setwd(paste0("/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/michelle/BED_files/Coverage_TSS_",bin.size,"bp_bin/normalizedBED_",bin.size,"bp_bin/"))
txt <- read.delim("bed_stat.txt", header=FALSE)
txt
# txt <- txt[-c(37,38),]
txt <- txt[c(1:6, 13:18, 7:12, 19:24, 31:36, 25:30),]
txt
rownames(txt) <- seq(1:36)
txt
num.of.lines=unique(txt[,3])
dat <- data.matrix(txt[,2])
rownames(dat) <- txt[,1]
class(dat)
rownames(dat) <- gsub("_R1.*", "", rownames(dat))
dat
mod_dat <- dat/1000000
colnames(mod_dat) <- "reads"

jpeg(paste0("readcount_",bin.size,"bp_barplot.jpeg"), height=600, width=900)
# custom.ylim=c(range(mod_dat)[2]/1.01, range(mod_dat)[2])
custom.ylim=c(0, range(mod_dat)[2])
par(mar=c(10,5,5,2))
barplot(t(mod_dat), beside=F, las=2, main=paste0("read counts (normalized_",bin.size,"bp_bin)\nnumber of ",bin.size,"bp bin=", num.of.lines), ylab= "# of reads (in million bp)", ylim=custom.ylim, xpd=FALSE)
dev.off()

