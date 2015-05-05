
system.dir="/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/"
# system.dir="Z:/users/shwang26/"
# a1 <- as.matrix(read.table(paste0(system.dir, "Michelle/MethylationData/methylatedGeneLists/2015-05-01/rep1_ash/TSS_age.stable.hypermethylation_10M_rep1_genelist.txt"))
# a2 <- as.matrix(read.table(paste0(system.dir, "Michelle/MethylationData/methylatedGeneLists/new/age_intermediate_methylated_genes_at_10months_new_cutoff.txt"))
# b1 <- as.matrix(read.table(paste0(system.dir, "Michelle/MethylationData/methylatedGeneLists/2015-05-01/rep2_ash/TSS_age.stable.hypermethylation_10M_rep2_genelist.txt"))
# 
# c1 <- as.matrix(read.table(paste0(system.dir, "Michelle/MethylationData/methylatedGeneLists/2015-05-01/rep2_chipseq/TSS_age.stable.hypermethylation_10M_rep2_genelist.txt"))
# 
# length(a1)
# length(a2)
# length(b1)
# intersect(a1, a2)


# AGE-STABLE 
#.2
#rep1
a1 <- as.matrix(read.table(paste0(system.dir, "Michelle/MethylationData/methylatedGeneLists/2015-05-04/rep1_agestable/TSS_age.stable.2.hypermethylation_10M_rep1_genelist.txt")))
#rep2
a2 <- as.matrix(read.table(paste0(system.dir, "Michelle/MethylationData/methylatedGeneLists/2015-05-04/rep2_agestable/TSS_age.stable.2.hypermethylation_10M_rep2_genelist.txt")))
dir.create(paste0(system.dir, "Michelle/MethylationData/methylatedGeneLists/2015-05-04/common_agestable"))
intersect_a1a2 <- intersect(a1, a2)
length(a1); length(a2); length(intersect_a1a2)
a1a <- setdiff(a1, intersect_a1a2)
a2a <- setdiff(a2, intersect_a1a2)
length(a1a); length(a2a); length(intersect_a1a2)
write.table(a1a, paste0(system.dir, "Michelle/MethylationData/methylatedGeneLists/2015-05-04/common_agestable/TSS_rep1.age.stable.2.hypermethylation_10M_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.
write.table(a2a, paste0(system.dir, "Michelle/MethylationData/methylatedGeneLists/2015-05-04/common_agestable/TSS_rep2.age.stable.2.hypermethylation_10M_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.
write.table(intersect_a1a2, paste0(system.dir, "Michelle/MethylationData/methylatedGeneLists/2015-05-04/common_agestable/TSS_rep1&2.age.stable.2.hypermethylation_10M_common_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) #row.names set to T to get the probe ids.


#.4
#rep1 
a3 <- as.matrix(read.table(paste0(system.dir, "Michelle/MethylationData/methylatedGeneLists/2015-05-04/rep1_agestable/TSS_age.stable.4.hypermethylation_10M_rep1_genelist.txt")))
#rep2
a4 <- as.matrix(read.table(paste0(system.dir, "Michelle/MethylationData/methylatedGeneLists/2015-05-04/rep2_agestable/TSS_age.stable.4.hypermethylation_10M_rep2_genelist.txt")))
dir.create(paste0(system.dir, "Michelle/MethylationData/methylatedGeneLists/2015-05-04/common_agestable"))
intersect_a3a4 <- intersect(a3, a4)
length(a3); length(a4); length(intersect_a3a4)
a3 <- setdiff(a3, intersect_a3a4)
a4 <- setdiff(a4, intersect_a3a4)
length(a3); length(a4); length(intersect_a3a4)
write.table(a3, paste0(system.dir, "Michelle/MethylationData/methylatedGeneLists/2015-05-04/common_agestable/TSS_rep1.age.stable.4.hypermethylation_10M_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) 
write.table(a4, paste0(system.dir, "Michelle/MethylationData/methylatedGeneLists/2015-05-04/common_agestable/TSS_rep2.age.stable.4.hypermethylation_10M_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) 
write.table(intersect_a3a4, paste0(system.dir, "Michelle/MethylationData/methylatedGeneLists/2015-05-04/common_agestable/TSS_rep1&2.age.stable.4.hypermethylation_10M_common_genelist.txt",sep=""), quote=F, row.names=F, col.names=F) 


# stable.10M
a_orig <- as.matrix(read.table(paste0(system.dir, "Michelle/MethylationData/methylatedGeneLists/new/age_stable_methylated_genes_at_10months_new_annotation.txt")))
length(a_orig)
# stable.10M.new
a_orig.new <- as.matrix(read.table(paste0(system.dir, "Michelle/MethylationData/methylatedGeneLists/new/age_stable_methylated_genes_at_10months_new_cutoff.txt")))
length(a_orig.new)

length(a1); length(a_orig); length(intersect(a_orig, a1))
a1b <- setdiff(a1, a_orig)
a_origa <- setdiff(a_orig, a1)
length(a1b); length(a_origa); length(intersect(a_orig, a1))

a3b <- setdiff(a3, a_orig.new)
a_origb <- setdiff(a_orig.new, a3)
length(a3b); length(a_origb); length(intersect(a_orig.new, a3))
