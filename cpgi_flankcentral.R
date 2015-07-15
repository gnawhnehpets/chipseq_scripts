dir <- "~/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/Michelle/MethylationData/methylatedGeneLists/compiled/rep1/rep1_cscmethylated_new_db.2/beta.values/"
setwd(dir)
beta <- read.table("TSS_cscmethylated_db.2_hypermethylation_10M_rep1_sb_beta.values.txt")
head(beta)
probes <- rownames(beta)
probes

cpgi <- read.table("~/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/Michelle/Required_Files/Annotations/hg19Data/hg19_cpgIslandExt.txt", sep="\t")[,-1]
colnames(cpgi) <- c("chrom", "chromStart", "chromEnd", "name", "length", "cpgNum", "gcNum", "perCpg", "perGc", "obsExp")
head(cpgi)
library(IlluminaHumanMethylation450k.db)
d <- load("~/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/Michelle/ash/jhu-methylation-scripts/R_image/annotation.RData")
d
head(CpGloc)
names(CpGloc)
length(probes)
length(CpGloc[probes])

filtered_probes <- CpGloc[probes]
filtered_probes <- (unlist(filtered_probes))
class(filtered_probes)
head(filtered_probes)
filtered_probes <- filtered_probes[grep("Island", filtered_probes)]
head(filtered_probes)
length(filtered_probes)
class(filtered_probes)
strsplit(filtered_probes, ":")
unlist(strsplit(filtered_probes, ":"))

mat <- read.table(text = filtered_probes, sep=":", colClasses = "character")
head(mat)
mat2 <- read.table( text = mat[,2], sep="-", colClasses = "character")
mat2
dat <- cbind(mat[,1], mat2, mat[,3])
colnames(dat) <- c("chr", "start", "end", "annotation")
head(dat)

#rownames(dat) = colnames(filtered_probes)
head(dat)

head(dat, n=30)

head(cpgi)
nrow(dat)
for(i in 1:nrow(dat)){
     chrom <- as.character(dat$chr[i])
     start <- as.numeric(dat$start[i])
     end <- as.numeric(dat$end[i])
#      print(paste0("chrom: ", chrom))
#      print(paste0("start: ", start))
#      print(paste0("end: ", end))     
#      tail(cpgi[which(cpgi$chrom=="chr1"),][,1:3])
#      head(cpgi[which(cpgi$chrom=="chr2"),][,1:3])
     cpgi[1:3]
     index <- which(cpgi$chrom==chrom)
# print(paste0(class(start), " ", start))
     #print(head(index))
     for(j in 1:length(index)){
         if(start > as.numeric(cpgi$chromStart[index[j]])){
              if(end < as.numeric(cpgi$chromStart[index[j]])){
#                    print(paste0("chrom: ", chrom))
#                    print(paste0("start: ", start))
#                    print(paste0("end: ", end))
                   print("hit!")
              }
         }
     }
}


head(cpgi)
chrom
which(cpgi$chrom==chrom)
