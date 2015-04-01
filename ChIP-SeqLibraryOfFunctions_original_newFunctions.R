#/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/
#/amber2/scratch/baylin/shwang/
#/Volumes/onc-analysis$/users/stephenhwang/

options(bitmapType='cairo') 

##################################################################
##################################################################
################# Functions for ChIP-seq analyses ################
##################################################################
##################################################################
# dir="Z:/users/shwang26/"
# dir="/amber2/scratch/baylin/shwang/"
dir="/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/"
##################################################################
# For coverageBed, a BED file called windows.bed is required to estimate the coverage data in defined windows (eg. 200bp across the genome, or defined intervals around TSS, etc) across the genome. For example:
## chr1 0 200
## chr1 200 400

# Note that the coverage data cannot be estimated for every nucleotide position across the genome due to memory limitations even in the cluster.
# Note that the coverage.data should have start-end columns of the form 0-200, 200-400, 400-600 (and not 0-199, 200-399, etc) for the >= and < operators to pick exactly the required number of bins (in functions fun.GetGeneCoverage and fun.GetGECoverage )

#Create such a file in R using hg18_chromInfo.txt, which creates the bins/windows in which the sequencing coverage has to be summarized
##################################################################

##################################################################
# Set path to the library taht will be used for R packages
## This needs to be explicitly stated after migrating cluster to JHPCE. Until this is resolved, use this path.
##################################################################
.libPaths()
lib.path <- "/home/steve/R/i686-pc-linux-gnu-library/3.1"
# lib.path <- "/jhpce/shared/community/compiler/gcc/4.4.7/R/3.0.x/lib64/R/site-library"
# lib.path <- "/Library/Frameworks/R.framework/Versions/3.1/Resources/library"
# lib.path <- "C:/Users/steve/Documents/R/win-library/3.1"
# lib.path <- "C:/Program Files/R/R-3.1.2/library" #remote desktop

#.libPaths()
##################################################################
#Important functions
##################################################################
#Function to read table given the file and path to file
fun.read.table.path <- function(file, path, ...){
     read.table(file.path(path, file), ...)
} #this is a function to read tables without changing the working directory


#Function that creates a bed format dataframe of the bins/windows in which the sequencing coverage has to be summarized.
fun.windows <- function(chromInfo, bin, chromInfoFilePath, windowsBedFileOut, ...){
     #this function creates a dataframe containing regular bins/windows in which the sequencing coverage has to be summarized; chromInfo is the character vector describing the filename that contains the genome information (chromosomes and their sizes) downloaded from USCSC; chromInfoFilePath is the path of the dataframe; bin is the size of the bins/windows (eg. 200 for 200 bp); windowsBedFileOut is character vector specifying output file name along with path where the output has to be saved.
     options("scipen"=100, "digits"=4) #scipen set to arbitrary 100 so that scientific notations (like 1.8 e+06) are not printed in the output file below. If this happens further functions in BedTools do not work.
     
     chromInfo <- fun.read.table.path(file=chromInfo, path=chromInfoFilePath, ...)
     chromInfo <- chromInfo[-grep("_", chromInfo[,1]), ] #remove rows that contain characters of the form "chr1_random".
     colnames(chromInfo) <- c("Chromosome", "Chromosome_Length", "source")
     c1 <- list()
     for(i in 1:nrow(chromInfo)){
          a <- seq(from=0, to=chromInfo[i, 2],  by=bin)
          b <- seq(from=bin, to=chromInfo[i, 2],  by=bin)
          u=min(length(a), length(b))
          chr.no <- strsplit(as.character(chromInfo[i, 1]), "r")[[1]][2]
          c1[[i]] <- data.frame(rep(paste("chr", chr.no, sep=""), u), a[1:u], b[1:u])
     }
     windows.bed <- do.call(rbind, c1)
     write.table(windows.bed, file=windowsBedFileOut, sep="\t", quote=F, row.names=F, col.names=F)
}


#Function that creates a bed format dataframe of the bins/windows in which the sequencing coverage has to be summarized in a defined region around the TSS
fun.TSSwindows <- function(goi.list, RegAroundTSS, bin, allGenesTSSWindowsBedFilePath, allGenesTSSWindowsFileOut){
     #this function creates a dataframe containing the intervals around TSS in which the sequencing coverage has to be summarized; goi.list is the list of dataframes containing the required information from biomart about the genes of interest (should be a object of class list); allGenesTSSWindowsBedFilePath is the path where the output file (dataframe) should be stored; bin is the size of the bins/windows (eg. 200 for 200 bp); allGenesTSSWindowsFileOut is character vector specifying output file name
     # NOTE: The RegAroundTSS used here is the length of the region in bp on either side of the TSS (this 10000 means 5000 up- and 5000 downstream. In other fucntions, the same variable name is used for the total region around the TSS. Unfortunately confusing. Keep in mind.
     
     x <- do.call(rbind, goi.list)
     
     # Get the nearst upstream TSS that is a multiple of 10, called TSSNearestUp10
     ## The idea is that the boundaries of the windows generated should be multiples of 10 so that if there is a TSS near another TSS, both will have the same windows, hence the coverage data in same windows. This prevents occurrence of overlapping windows, which would have been teh case if the TSS were used directly. This maneuver is required for fun.average_heat.plots to work properly (when a region within +/- RegAroundTSS/2 is used to get coverage data using fun.GetGeneCoverage) such that the number of intervals is same for each TSS region in question. Having overlapping windows gives variable number of intervals and the data for such genes are lost during evaluation of the condition "if (length(gene.coverage) == num.bins)".
     xWithTSSNearestUp10 <- x
     for(i in 1:nrow(xWithTSSNearestUp10)){
          tss = xWithTSSNearestUp10$TSS[i]
          if((tss %% 10) != 0){
               tssNearestUp10 = tss - (tss %% 10)
          } else {
               tssNearestUp10 = tss
          }
          xWithTSSNearestUp10[i, "TSSNearestUp10"] <- tssNearestUp10
     }
     
     # Create bed file
     x.bed <- data.frame(xWithTSSNearestUp10$Chromosome, 
                         xWithTSSNearestUp10$TSSNearestUp10-(RegAroundTSS+(2*bin)),
                         xWithTSSNearestUp10$TSSNearestUp10+(RegAroundTSS+(2*bin)),
                         xWithTSSNearestUp10$Gene) #dataframe organized in bed format; TSSNearestUp10 +/- RegAroundTSS+(2*bin) computed to avoid any issues later on  while getting the coverage data in the central bin, considered as 0th bin, and upstream and downstream bins.
     #although colnames is not written out to the output file below, it is labeled here for convenience of understanding the script
     colnames(x.bed) <- c("Chromosome", "TSSNearestUp10-RegAroundTSS", "TSSNearestUp10+RegAroundTSS", "Gene")
     
     # Define windows
     # Create windows bed files for each chromosome and rbind them
     chr = unique(x.bed$Chromosome)
     windowsBedByChr <- lapply(chr, FUN=function(u){
          # Split x.bed by chromosomes 
          x.bed.u <- x.bed[x.bed[,"Chromosome"] %in% u, ]
          windowsBed.u <- list()
          for(i in 1:nrow(x.bed.u)){
               binStart <- seq(from=x.bed.u[i,"TSSNearestUp10-RegAroundTSS"], 
                               to=x.bed.u[i, "TSSNearestUp10+RegAroundTSS"],  by=bin)
               binEnd <- seq(from=x.bed.u[i, "TSSNearestUp10-RegAroundTSS"] + bin, 
                             to=x.bed.u[i, "TSSNearestUp10+RegAroundTSS"] + bin,  by=bin)
               windowsBed.u[[i]] <- data.frame(rep(as.character(u), length(binStart)), binStart, binEnd)
          }
          windowsBed.u <- do.call(rbind, windowsBed.u)
          colnames(windowsBed.u) <- c("Chromosome", "BinStart", "BinEnd")
          # Only retain unique window intervals (the interals will be repeated because if two TSS are nearby, each will contribute to the windows)
          windowsBed.u <- windowsBed.u[!duplicated(windowsBed.u$BinStart), ]
          return(windowsBed.u)
     })
     windowsBED <- do.call(rbind, windowsBedByChr)
     
     write.table(windowsBED, paste(allGenesTSSWindowsBedFilePath, allGenesTSSWindowsFileOut, sep="/"), sep="\t", quote=F, row.names=F, col.names=F)
}


#############################################################
#############################################################
############ Functions for plotting coverage data ###########
#############################################################
#############################################################


############################################################
#Required functions
############################################################
#Function to retrieve coverage data for a histone modification within a defined region up and downstream from the TSS for each gene in a gene list
fun.GetGeneCoverage <- function(genelist.info, coverage.data, bin=200, RegAroundTSS=10000, chr.prefix.chromosome=T, version){
     # This function extracts the coverage data in a defined region around the TSS (RegAroundTSS); genelist.info is the dataframe consisting the required information ("Gene", "Chromosome", "TSS", "Strand") about the genes (typically obtained from biomaRt); coverage.data is the coverage data in bed format; bin is the window size in which the coverage data is binned (default is 200 bp); RegAroundTSS is the +/-region around the TSS for which the coverage in the defined bin size has to be extracted (default is 10000bp, i.e. +/-5000 bp around TSS); chr.prefix.chromosome takes values T or F for whether or not a "chr" prefix needs to be added to the chromosome values in genelist.info. This needs to be done if the genelist.info input is directly derived from biomaRt in which case the chromosomes are of the form 1,2,3..., which needs to be converted to chr1, chr2, chr3... . If the chromosome values are already of the latter value, then assign chr.prefix.chromosome=F
     
     # Define num.bins and window.bp:
     ## num.bins is the number of coverage bins around the TSS that that has to be obtained given the RegAroundTSS and bin.
     ## window.bp is the length in basepairs used to calculate the upstream and downstream coordinates from the TSS within which the coverage data will be obtained.
     ## If RegAroundTSS/bin is even, add 1 to num.bins to make it odd. This is done because the coverage is plotted for equal number of bins from a central bin. So the total number of bins has to be odd. (%% is the modulo function)  
     if((RegAroundTSS/bin) %% 2 == 0){
          num.bins <- (RegAroundTSS/bin)+1
     } else {
          num.bins <- (RegAroundTSS/bin)
     }
     window.bp <- ((num.bins-1)/2)*bin
     
     # split genelist.info by chromosomes
     # split coverage data by chromosomes
     # then get coverage around each chromosome for each gene
     chr <- paste("chr", c(1:22, "X", "Y"), sep="")
     if(version=="hg19"){
          if(chr.prefix.chromosome==F){
               genelist.info$Chromosome <- paste("chr", genelist.info$Chromosome, sep="") #this needs to be done because the dataframe fed as genelist.info to this function, normally obtained using biomaRt, contains only numeric value as chromsome numbers (eg. 22) while the ChIP-seq coverage data contains chromsosomes of the form chr22.
          }
     }
     if(version=="hg18"){
          if(chr.prefix.chromosome==T){
               genelist.info$Chromosome <- paste("chr", genelist.info$Chromosome, sep="") #this needs to be done because the dataframe fed as genelist.info to this function, normally obtained using biomaRt, contains only numeric value as chromsome numbers (eg. 22) while the ChIP-seq coverage data contains chromsosomes of the form chr22.
          }
     }
     
     #   # Get the nearest upstream coordinate that is a multiple of 10 as in fun.TSSwindows
     #   ## THIS IS NOT REQUIORED FOR THIS FUNCTION
     #   for(i in 1:nrow(genelist.info)){
     #     tss = genelist.info$TSS[i]
     #     if((tss %% 10) != 0){
     #       tssNearestUp10 = tss - (tss %% 10)
     #     } else {
     #       tssNearestUp10 = tss
     #     }
     #     genelist.info[i, "TSSNearestUp10"] <- tssNearestUp10
     #   }
     #   
     
     #split object into list based off 'chrX'
     genelist.info <- lapply(chr, FUN=function(u){genelist.info[genelist.info[,"Chromosome"] %in% u, ]}) # splits genelist.info by chromosomes and creates a list (each dataframe object in list holds information about genes in one chromosome)
     names(genelist.info) <- chr
     coverage.data <- lapply(chr, FUN=function(u){coverage.data[coverage.data[,"Chromosome"] %in% u, ]})
     names(coverage.data) <- chr
     
     #for each chromosome...
     newlist <- lapply(chr, FUN=function(u){
          genelist.info.u <- genelist.info[[u]]
          
          #Get coverage data only if there is a gene in a particular iteration of a chromosome. 'if' condition used here so that the script 
          #does not stop here  with an error if there is no gene in the query belonging to a particular chromosome.
          if(dim(genelist.info.u)[1] > 0){
               coverage.data.u <- coverage.data[[u]]
               cov.chr <- sapply(1:nrow(genelist.info.u), FUN=function(x){
                    #           print(x)
                    if(genelist.info.u[x, "Strand"] == 1){
                         up <- genelist.info.u[x, "TSS"] - window.bp
                         down <- genelist.info.u[x, "TSS"] + window.bp
                         
                         #pick the bins up and downstream of the TSS window (set to the bin size, i.e. the bin contaiing the TSS coordinate is considered
                         #the central bin, and the coverage in the required number of bins up- and down-stream from the central bin will be extracted)
                         gene.coverage <- coverage.data.u[(coverage.data.u[,"start"] >= up - bin) & (coverage.data.u[,"end"] < down + bin), "Seq_tags"]
                         # The number of bins extracted for gene.coverage should be equal to num.bins. For eg. for 10000bp and 200 bp bin, it will be 51. 
                         # If for some reason the number of coverage bins do not correspond to this number, remove them from the analysis (thus retain 
                         # only genes having number of coverage bins= num.bins). Note that the coverage.data should have start-end columns of the form 
                         # 0-200, 200-400, 400-600 (and not 0-199, 200-399, etc) for the >= and < operators to pick exactly the required number of bins.
                         # if the TSS is at an end of a chromosome (or the data for that region is lacking), then report these as NA's so that the 
                         # plotting functions later on work.
                         #                print(paste0("strand==1: ", class(gene.coverage)))
                         if (length(gene.coverage) == num.bins){
                              gene.coverage=gene.coverage
                         }else{
                              gene.coverage=rep(NA,num.bins)
                         }
                    }else{
                         
                         up <- genelist.info.u[x, "TSS"] + window.bp
                         down <- genelist.info.u[x, "TSS"] - window.bp
                         #pick the bins up and downstream of the TSS window as above (note that the up and down definitions are changed here because of the transcript orientation).
                         #Note that the coverage.data should have start-end columns of the form 0-200, 200-400, 400-600 (and not 0-199, 200-399, etc) for the >= and < operators to pick exactly the required number of bins.
                         #                 gene.coverage <- vector()
                         gene.coverage <- coverage.data.u[(coverage.data.u[,"end"] < up + bin) & (coverage.data.u[,"start"] >= down - bin), c("start", "end", "Seq_tags")]
                         gene.coverage <- gene.coverage[order(gene.coverage[,"start"], decreasing=T), ] #this is done to have transcripts in the minus strand to be oriented from upstream to downstream
                         gene.coverage <- gene.coverage[, "Seq_tags"]
                         #               print(paste0("else:      ", class(gene.coverage)))
                         if (length(gene.coverage) == num.bins){
                              gene.coverage=gene.coverage
                         }else{
                              gene.coverage=rep(NA,num.bins)
                         }
                    }
                    return(gene.coverage)
               })
               return(cov.chr)
          }
     })
     coverage <- do.call(cbind, newlist)
}


#Function for getting genesets grouped by different percentiles of expression
fun.GetGenes.quantiles <- function(probs, dat){
     #Gets genesets grouped by different percentiles of expression; prob is the numeric vector of probabilities with values in [0,1]; dat is the gene expression table with 1st column as gene names and 2nd column as gene expression values; the output is an object of class list containing the list of genes of different expression percentiles
     quantiles <- quantile(dat[,2], probs=probs)
     x <- lapply(c(1:(length(quantiles)-1)), FUN=function(u){
          return(as.character(dat[dat[,2] >= quantiles[u] & dat[,2] < quantiles[u+1], 1]))
     })
     return(x)
}


# Function for getting information about gene list of interest (goi) from biomaRt
fun.genelist.info <- function(goi.list, genelist.names, goi.list.type, version){
     #Gets required information from biomaRt about genes in a list of genesets (goi.list) and outputs it as list of dataframes; genelist.names is a vector names that the output dataframes should have; goi.list.type is whether the genes in the goi.list are HGNC symbols ("hgnc_symbol") or RefSeq ids ("RefSeq_ID"); version denotes which version of genome should be used, eg. "hg18".
     
     library(biomaRt, lib.loc=lib.path)
     library(RCurl, lib.loc=lib.path)
     if(version=="hg18"){
          print("version hg18")
          #listMarts(host="may2009.archive.ensembl.org", path="/biomart/martservice",archive=FALSE)
          ensembl54 = useMart(host="may2009.archive.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl") #Get the hg18 version in biomart
          #listAttributes(ensembl54); #listFilters(ensembl54)
          x <- list()
          for(i in 1:length(goi.list)){
               genelist <- unique(goi.list[[i]])
               if(goi.list.type=="hgnc_symbol"){
                    m <- getBM(attributes=c("hgnc_symbol", "hgnc_curated_gene_name", "external_gene_id", "transcript_start", "transcript_end", "chromosome_name", "strand", "refseq_dna", "gene_biotype"), filters = c("hgnc_symbol"), values = genelist, mart = ensembl54)
               }
               if(goi.list.type=="RefSeq_ID"){
                    m <- getBM(attributes=c("hgnc_symbol", "hgnc_curated_gene_name", "external_gene_id", "transcript_start", "transcript_end", "chromosome_name", "strand", "refseq_dna", "gene_biotype"), filters = c("refseq_dna"), values = genelist, mart = ensembl54)
               }
               m <- m[m$gene_biotype %in% "protein_coding", ] #Retain only protein coding genes
               m <- m[!(m[,"refseq_dna"] %in% ""), ] #Retain only genes with a RefSeq ID.
               TSS <- sapply(c(1:length(m[,1])), FUN=function(u) {ifelse(m$strand[u]==1, m[u,4], m[u,5])})
               TES <- sapply(c(1:length(m[,1])), FUN=function(u) {ifelse(m$strand[u]==1, m[u,5], m[u,4])})
               genelist.info <- cbind.data.frame(m$hgnc_symbol, m$refseq_dna, m$chromosome_name, TSS, TES, m$strand, m$external_gene_id, stringsAsFactors=F) #table containing the accession nos. and TSS
               colnames(genelist.info) <- c("Gene", "Accession", "Chromosome", "TSS", "TES", "Strand", "ExternalGeneID")
               #Note: Some TSS and TES (transcription end site) are duplicated. Keep in mind.
               ## Some TSS and TES (transcription end site) are duplicated. Even some of the RefSeq ID (Accession) are duplicated. So retain only genes that have unique TSS (as we are mostly interested in the histone modifications around TSS)
               chr <- unique(genelist.info[,"Chromosome"])
               genelist.info.uniqueTSS <- lapply(chr, FUN=function(i){
                    genelist.info.chr <- genelist.info[genelist.info$Chromosome %in% i, ]
                    genelist.uniqueTSS <-  genelist.info.chr[!(duplicated(genelist.info.chr$TSS)), ]
                    return(genelist.uniqueTSS)
               })
               genelist.info <- do.call(rbind, genelist.info.uniqueTSS)
               # => This is the final genelist.info dataframe that contains necessary information about the given gene list
               
               assign(paste(genelist.names[i]), genelist.info)
               x[[i]] <- get(genelist.names[[i]])
               #assign and get commands need not be used for assigning genelist.info to the list object x, and instead genelist.info can be directly added to x.  These commands are used here because they ae cool commands to keep in mind! ;-)
          }
          names(x) <- genelist.names
          return(x)
     }
     if(version=="hg19"){
          print("version: hg19")
          chrom <- c(1:22, "X", "Y")
          ensembl54 = useMart(host='feb2014.archive.ensembl.org', biomart='ENSEMBL_MART_ENSEMBL', dataset="hsapiens_gene_ensembl") #Get the hg19 version in biomart
          x <- list()
          for(i in 1:length(goi.list)){
               genelist <- unique(goi.list[[i]])
               if(goi.list.type=="hgnc_symbol"){
                    m <- getBM(attributes=c("hgnc_symbol", "external_gene_id", "transcript_start", "transcript_end", "chromosome_name", "strand", "gene_biotype"), filters = c("chromosome_name"), values = chrom, mart = ensembl54)
               }
               if(goi.list.type=="RefSeq_ID"){
                    m <- getBM(attributes=c("hgnc_symbol", "hgnc_curated_gene_name", "external_gene_id", "transcript_start", "transcript_end", "chromosome_name", "strand", "refseq_dna", "gene_biotype"), filters = c("refseq_dna"), values = genelist, mart = ensembl54)
               }
               m <- m[m$gene_biotype %in% "protein_coding", ] #Retain only protein coding genes
               index <- which(m$hgnc %in% genelist)
               m <- m[index,]
               m$chromosome_name <- paste("chr", m$chromosome_name, sep="") #chromosome names converted to the form "chr1" because this is how the chromosomes are annotated in external resources (however, not really required for the Infinium annotation as the chromosomes are of the fom 1,2,...X, Y)
               TSS <- sapply(c(1:length(m[,1])), FUN=function(u) {ifelse(m$strand[u]==1, m[u,"transcript_start"], m[u,"transcript_end"])})
               TES <- sapply(c(1:length(m[,1])), FUN=function(u) {ifelse(m$strand[u]==1, m[u,"transcript_end"], m[u,"transcript_start"])})
               genelist.info <- cbind.data.frame(m$hgnc_symbol, m$chromosome_name, TSS, TES, m$strand, m$external_gene_id, stringsAsFactors=F) #table containing the accession nos. and TSS
               colnames(genelist.info) <- c("Gene", "Chromosome", "TSS", "TES", "Strand", "ExternalGeneID")
               #Note: Some TSS and TES (transcription end site) are duplicated. Keep in mind.
               ## Some TSS and TES (transcription end site) are duplicated. Even some of the RefSeq ID (Accession) are duplicated. So retain only genes that have unique TSS (as we are mostly interested in the histone modifications around TSS)
               chr <- unique(genelist.info[,"Chromosome"])
               genelist.info.uniqueTSS <- lapply(chr, FUN=function(i){
                    genelist.info.chr <- genelist.info[genelist.info$Chromosome %in% i, ]
                    genelist.uniqueTSS <-  genelist.info.chr[!(duplicated(genelist.info.chr$TSS)), ]
                    return(genelist.uniqueTSS)
               })
               genelist.info <- do.call(rbind, genelist.info.uniqueTSS)
               # => This is the final genelist.info dataframe that contains necessary information about the given gene list
               assign(paste(genelist.names[i]), genelist.info)
               x[[i]] <- get(genelist.names[[i]])
               #assign and get commands need not be used for assigning genelist.info to the list object x, and instead genelist.info can be directly added to x.  These commands are used here because they ae cool commands to keep in mind! ;-)
          }
          names(x) <- genelist.names
          return(x)
     }
}
#Example: x <- fun.genelist.info(goi.list=genes.quantiles, genelist.names=names(genes.quantiles), goi.list.type="hgnc_symbol")


#Function for getting the gene information about all genes in the genome from Biomart
fun.genelist.info_allGenes_biomaRt <- function(version, hg19UCSCGeneAnnotations=NULL, hg19UCSCGeneAnnotationsPath=NULL){
     #This function downloads the whole transcriptome to get the TSS and other information for each gene
     # Args description
     ## version denotes which version of genome should be used, viz. "hg18", "hg19", or "hg19-UCSC".
     ## hg19UCSCGeneAnnotations is character the name of the text file that contains teh annotations. Defaults to NULL
     ## hg19UCSCGeneAnnotationsPath is the path to the annotation text file downloaded from UCSC. Defaults to NULL
     
     # Download the whole transcriptome to get the TSS
     if(exists("lib.path")){
          # Occasioanlly the JHPCE cluster asks for lib.loc. Take this value from lib.path if exists in the environment. Else, in the case lib.loc is not required, and hence lib.path is not explicitly defined, use default library path.
          library(biomaRt, lib.loc=lib.path)
          library(RCurl, lib.loc=lib.path)  
     } else {
          library(biomaRt)
          library(RCurl)
     }
     
     
     # Set scipen to arbitrary 100 so that scientific notations (like 1.8 e+06) are not printed in the output file. If this happens further functions in BedTools do not work.
     options("scipen"=100, "digits"=4)
     
     # Set chromosome names
     chrom <- c(1:22, "X", "Y")
     
     if(version=="hg18"){
          #listMarts(host="may2009.archive.ensembl.org", path="/biomart/martservice",archive=FALSE)
          ensembl = useMart(host="may2009.archive.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl") #Get the hg18 version in biomart
          
          #listAttributes(ensembl); #listFilters(ensembl)
          m <- getBM(attributes=c("hgnc_symbol", "hgnc_curated_gene_name", "external_gene_id", "transcript_start", "transcript_end", "chromosome_name", "strand", "refseq_dna", "gene_biotype"), filters = c("chromosome_name"), values = chrom, mart = ensembl)
          m <- m[m$gene_biotype %in% "protein_coding", ] #Retain only protein coding genes
          m <- m[!(m[,"refseq_dna"] %in% ""), ] #Retain only genes with a RefSeq ID.
          TSS <- sapply(c(1:length(m[,1])), FUN=function(u) {ifelse(m$strand[u]==1, m[u,"transcript_start"], m[u,"transcript_end"])})
          TES <- sapply(c(1:length(m[,1])), FUN=function(u) {ifelse(m$strand[u]==1, m[u,"transcript_end"], m[u,"transcript_start"])})
          m$chromosome_name <- paste("chr", m$chromosome_name, sep="") #chromosome names converted to the form "chr1" because this is how the chromosomes are annotated in external resources, like CpG island dataframe
          
          genelist.info <- cbind.data.frame(m$hgnc_symbol, m$refseq_dna, m$chromosome_name, TSS, TES, m$strand, m$external_gene_id, stringsAsFactors=F) #table containing the accession nos. and TSS
          colnames(genelist.info) <- c("Gene", "Accession", "Chromosome", "TSS", "TES", "Strand", "ExternalGeneID")
          
          ## Some TSS and TES (transcription end site) are duplicated. Even some of the RefSeq ID (Accession) are duplicated. So retain only genes that have unique TSS (as we are mostly interested in the histone modifications around TSS)
          chr <- unique(genelist.info[,"Chromosome"])
          genelist.info.uniqueTSS <- lapply(chr, FUN=function(i){
               genelist.info.chr <- genelist.info[genelist.info$Chromosome %in% i, ]
               genelist.uniqueTSS <-  genelist.info.chr[!(duplicated(genelist.info.chr$TSS)), ]
               return(genelist.uniqueTSS)
          })
          genelist.info <- do.call(rbind, genelist.info.uniqueTSS)
          # => This is the final genelist.info dataframe that contains necessary information about the whole transcriptome
          return(genelist.info)
     }
     
     if(version=="hg19"){
          #listMarts(host='feb2014.archive.ensembl.org')
          ensembl=useMart(host='feb2014.archive.ensembl.org', biomart='ENSEMBL_MART_ENSEMBL', dataset="hsapiens_gene_ensembl") #Get the hg19 version in biomart
          
          m <- getBM(attributes=c("hgnc_symbol", "external_gene_id", "transcript_start", "transcript_end", "chromosome_name", "strand", "gene_biotype"), filters = c("chromosome_name"), values = chrom, mart = ensembl)
          TSS <- sapply(c(1:length(m[,1])), FUN=function(u) {ifelse(m$strand[u]==1, m[u,"transcript_start"], m[u,"transcript_end"])})
          TES <- sapply(c(1:length(m[,1])), FUN=function(u) {ifelse(m$strand[u]==1, m[u,"transcript_end"], m[u,"transcript_start"])})
          m$chromosome_name <- paste("chr", m$chromosome_name, sep="") #chromosome names converted to the form "chr1" because this is how the chromosomes are annotated in external resources (however, not really required for the Infinium annotation as the chromosomes are of the fom 1,2,...X, Y)
          
          genelist.info <- cbind.data.frame(m$hgnc_symbol, m$chromosome_name, TSS, TES, m$strand, m$external_gene_id, stringsAsFactors=F) #table containing the accession nos. and TSS
          colnames(genelist.info) <- c("Gene", "Chromosome", "TSS", "TES", "Strand", "ExternalGeneID")
          
          ## Some TSS and TES (transcription end site) are duplicated. So retain only genes that have unique TSS
          chr <- unique(genelist.info[,"Chromosome"])
          genelist.info.uniqueTSS <- lapply(chr, FUN=function(i){
               genelist.info.chr <- genelist.info[genelist.info$Chromosome %in% i, ]
               genelist.uniqueTSS <-  genelist.info.chr[!(duplicated(genelist.info.chr$TSS)), ]
               return(genelist.uniqueTSS)
          })
          genelist.info <- do.call(rbind, genelist.info.uniqueTSS)
          # => This is the final genelist.info dataframe that contains necessary information about the whole transcriptome
          
          return(genelist.info)
     }
     
     if(version=="hg19-UCSC"){
          #For hg19-UCSC use the annotation downloaded from UCSC. This is done because the TSS coordinates in grch37 archive in biomart (feb2014.archive.ensembl.org), corresponding to hg19, is slightly shifted from the UCSC-hg19 assembly.
          # STEP-1: Read the UCSC gene annotation.
          hg19UCSCGeneAnnotations <- file.path(hg19UCSCGeneAnnotationsPath, hg19UCSCGeneAnnotations)
          hg19UCSCGeneAnnotations <- read.table(hg19UCSCGeneAnnotations, header=T, sep="\t", skip=1, stringsAsFactors=F)
          colnames(hg19UCSCGeneAnnotations) <- c("UCSCGeneName", "Chromosome", "Strand", "TSS", "TES", "CodingRegionStart", "CodingRegionEnd", "proteinID", "alignID", "geneSymbol")
          head(hg19UCSCGeneAnnotations)
          
          ## Some TSS and TES (transcription end site) are duplicated. So retain only genes that have unique TSS
          chr <- unique(hg19UCSCGeneAnnotations[,"Chromosome"])
          hg19UCSCGeneAnnotations.uniqueTSS <- lapply(chr, FUN=function(i){
               hg19UCSCGeneAnnotations.chr <- hg19UCSCGeneAnnotations[hg19UCSCGeneAnnotations$Chromosome %in% i, ]
               genelist.uniqueTSS <-  hg19UCSCGeneAnnotations.chr[!(duplicated(hg19UCSCGeneAnnotations.chr$TSS)), ]
               return(genelist.uniqueTSS)
          })
          hg19UCSCGeneAnnotations <- do.call(rbind, hg19UCSCGeneAnnotations.uniqueTSS)
          
          #Rename hg19UCSCGeneAnnotations as genelist.info as this was the name used conventionally for this table in lot of other scripts
          ## convert the Strand to numeric form of 1, -1 form
          genelist.info <- hg19UCSCGeneAnnotations
          genelist.info$Strand <- as.numeric(genelist.info$Strand %in% "+") + -as.numeric(genelist.info$Strand %in% "-")
          rm(hg19UCSCGeneAnnotations)
          # => This is the final genelist.info dataframe that contains necessary information about the whole transcriptome
          colnames(genelist.info)[which(colnames(genelist.info) == "geneSymbol")]  <- "Gene"
          return(genelist.info)
     }
}
#Example: genelist.info <- fun.genelist.info_allGenes_biomaRt(version="hg19-UCSC", hg19UCSCGeneAnnotations="hg19UCSCGeneAnnotaions.txt", hg19UCSCGeneAnnotationsPath="/amber2/scratch/baylin/Hari/Required_Files/Annotations/hg19Data")


# Function for plotting coverage data around the TSS as average values and RATIO to input as cartesian plots and heatmaps
fun.average_heat.plots <- function(genelist.info, 
                                   coverage_files, 
                                   input_coverage_file, 
                                   coverage_files.dir, 
                                   RegAroundTSS, 
                                   bin, 
                                   chr.prefix.chromosome, 
                                   plot.Directory, 
                                   whichgenelist="stable.10M",
                                   logHeatmap=F, 
                                   colRangeHeatmap="white-black", 
                                   version,
                                   h3k4.max=1500,
                                   h3k27.max=250,
                                   dnmt.max=50,
                                   ezh2.max=120,
                                   inp.max=50,
                                   h3.max=50
                                   ){
     # Plots the coverage data as raw reads and ratio plots around the TSS as average plot and heatmap; genelist.info is a list of dataframes containing the required info for the genes of interest (the objects in this list should have a name so that the plots can take these names); coverage_files is a vector of file names of the ChIP-seq coverage data; input_coverage_file is the file name of the input coverage file; coverage_files.dir is the directory in which the converage files data are stored; RegAroundTSS is the total region upstream and downstream from TSS for which coverage data will be obtained and plotted (for eg. RegAroundTSS of 10000 will get coverage data 5000 bp up and downstream from the TSS). It should be a multiple of the bin size; bin is the window sizes in which the coverage of ChIP-seq data is summarized; chr.prefix.chromosome is to be passed on to fun.GetGeneCoverage, and takes values T or F for whether or not a "chr" prefix needs to be added to the chromosome values in genelist.info (see fun.GetGeneCoverage for details). plot.Directory is a character vector assigning the directory in which the plots should be stored. logHeatmap is T or F indicating whether or not heatmap should be plotted for log values; colRangeHeatmap is either "white-black" or "red-black", to indicate the color palette to use for the heatmaps, defaults to "white-black".
     library(gplots, lib.loc=lib.path)
     
     # Create directory for saving plots
     dir.create(plot.Directory)
     custom.annotation <- vector()
     
     if(whichgenelist=="stable.10M"){
          load(paste0(dir, "Michelle/BED_files/Coverage_TSS_",bin,"bp_bin/normalizedBED_",bin,"bp_bin/outputdir/c10d_stable.10M_",bin,"bp_ann.bar.Rdata"))
          custom.annotation <- c10d.10M.ann.bar
     }
     if(whichgenelist=="intermediate.10M"){
          load(paste0(dir, "Michelle/BED_files/Coverage_TSS_",bin,"bp_bin/normalizedBED_",bin,"bp_bin/outputdir/c10d_intermediate.10M_",bin,"bp_ann.bar.Rdata"))
          custom.annotation <- c10d.10M.ann.bar
     }
     if(whichgenelist=="stable.10M.new"){
          load(paste0(dir, "Michelle/BED_files/Coverage_TSS_",bin,"bp_bin/normalizedBED_",bin,"bp_bin/outputdir/c10d_new_stable.10M_",bin,"bp_ann.bar.Rdata"))
          custom.annotation <- c10d.10M.ann.bar
     }
     if(whichgenelist=="intermediate.10M.new"){
          load(paste0(dir, "Michelle/BED_files/Coverage_TSS_",bin,"bp_bin/normalizedBED_",bin,"bp_bin/outputdir/c10d_new_intermediate.10M_",bin,"bp_ann.bar.Rdata"))
          custom.annotation <- c10d.10M.ann.bar
     }
     if(whichgenelist=="highly.expressed"){
          load(paste0(dir, "Michelle/BED_files/Coverage_TSS_",bin,"bp_bin/normalizedBED_",bin,"bp_bin/outputdir/c10d_highly.expressed_",bin,"bp_ann.bar.Rdata"))
          custom.annotation <- c10d.10M.ann.bar
     }
     # dim(c10d.10M.ann.bar)
     # dim(custom.annotation)
     
     custom.annotation
     # Define num.bins and int: 
     ## num.bins is the number of coverage bins around the TSS that will be plotted given RegAroundTSS and bin.
     ## If RegAroundTSS/bin is even, add 1 to num.bins to make it odd. This is done because the coverage is plotted for equal number of bins from a central bin. So the total number of bins has to be odd. (%% is the modulo function)
     ## int is the position of bins around the TSS in which coverage is being plotted.
     ## Note that RegAroundTSS should be a multiple of bin (see definition of bin).
     if((RegAroundTSS/bin) %% 2 == 0){
          num.bins <- (RegAroundTSS/bin)+1
          int <- seq(from=-(RegAroundTSS/(bin*2)), to=(RegAroundTSS/(bin*2)), by=1)*bin
     } else {
          num.bins <- (RegAroundTSS/bin)
          int <-  seq(from=-((num.bins-1)/2), to= (num.bins-1)/2, by=1)*bin
     }
     
     # Read the input coverage file
     if(exists("input_coverage_file"))
     {
          input.coverage <- read.table(file.path(coverage_files.dir, input_coverage_file), header=F, sep="\t")
          colnames(input.coverage) <- c("Chromosome", "start", "end", "Seq_tags")
     } else 
     {
          print("INPUT COVERAGE DATA ABSENT")
     }
     
     for(i in coverage_files){
          
          # i <- coverage_files[4] # for troubleshooting
          # coverage_files
          # i
          ############################################################
          ##Load the ChIP_seq coverage data
          ############################################################
          #Read the bed file containing the coverage in bins defined by the bin size using the coverageBed command in BEDtools
          #Retain rows that have bins = bin bp (some bins with are > bin size were generated by coverageBed in the early iterations. This problem has been resolved and thus the scripts to retain bins that are == bin can be deleted)
          chip.coverage <- read.table(file.path(coverage_files.dir, i), header=F, sep="\t")
          #rm.col <- chip.coverage[,3] - chip.coverage[,2]
          #chip.coverage <- chip.coverage[which(rm.col == bin), ]
          colnames(chip.coverage) <- c("Chromosome", "start", "end", "Seq_tags")
          chipCoverage.TSSAverageList <- list() # holds the ChIP average values object for genes in genelist.info
          inputCoverage.TSSAverageList <- list() # holds the Input average values object for genes in genelist.info
          
          for(g in 1:length(genelist.info)){
               g=1 #for troubleshooting
               #Retrieve regions within RegAroundTSS of TSS for each gene for each mark
               chipCoverage.TSS <- fun.GetGeneCoverage(genelist.info=genelist.info[[g]], coverage.data=chip.coverage, RegAroundTSS=RegAroundTSS, bin=bin, chr.prefix.chromosome=chr.prefix.chromosome, version=version)
               inputCoverage.TSS <- fun.GetGeneCoverage(genelist.info=genelist.info[[g]], coverage.data=input.coverage, RegAroundTSS=RegAroundTSS, bin=bin, chr.prefix.chromosome=chr.prefix.chromosome, version=version)
               
               if(dim(t(na.omit(t(chipCoverage.TSS))))[2] >= 2){#there should be coverage data for at least two genes (automatically it means that there is atleast two-gene coverage data for the input also). na.omit removes rows that contain NA in a data frame. t used here because columns containing NAs have to be removed - thus after omitting columns with NA, coverage data should have data for atleast 2 gene  to satisfy the 'if' statement.
                    ## Calculate average ChIP coverage
                    chipCoverage.TSSAverage <- apply(chipCoverage.TSS, 1, mean, na.rm=T)
                    chipCoverage.TSSAverageList[[g]] <- chipCoverage.TSSAverage
                    names(chipCoverage.TSSAverageList)[g] <- paste(i, names(genelist.info)[g], sep="-")
                    
                    ## Calculate average input coverage
                    inputCoverage.TSSAverage <- apply(inputCoverage.TSS, 1, mean, na.rm=T)
                    ## If any of the inputCoverage.TSSAverage vales is 0, replace it with the minimum value
                    inputCoverage.TSSAverage <- replace(inputCoverage.TSSAverage, which(inputCoverage.TSSAverage == 0), min(inputCoverage.TSSAverage[which(inputCoverage.TSSAverage != 0)]))
                    inputCoverage.TSSAverageList[[g]] <- inputCoverage.TSSAverage
                    names(inputCoverage.TSSAverageList)[g] <- paste(i, names(genelist.info)[g], sep="-")
                    
                    # Plot coverage (raw seq reads and ratio to input of seq reads) as heatmap plot
                    ## Create plot name for raw heatmaps
                    heat_plot_raw <- file.path(plot.Directory, paste(i, names(genelist.info)[g], "heat_plot_raw.jpeg", sep="-"))
                    
                    ## Create plot name for RATIO heatmaps
                    heat_plot_RATIO <- file.path(plot.Directory, paste(i, names(genelist.info)[g], "heat_plot_RATIO.jpeg", sep="-"))
                    
                    # Generate matices for heatmaps
                    x.trnposed <- t(chipCoverage.TSS)
                    ratioToInp.x.trnposed <- t(chipCoverage.TSS/inputCoverage.TSSAverage) #ratio of seq reads to average of input for this set of genes
                    if(logHeatmap == T){
                         x.trnposed <- log2(x.trnposed)
                         ratioToInp.x.trnposed <- log2(ratioToInp.x.trnposed)
                    }
                    
                    colnames(x.trnposed) <- colnames(ratioToInp.x.trnposed) <- int
                    rownames(x.trnposed) <- rownames(ratioToInp.x.trnposed) <- genelist.info[[g]]$hgnc_symbol
                    if(colRangeHeatmap == "red-black"){
                         numberOfColors.x.trnposed <- colorRampPalette(c("red", "black"))(round(range(x.trnposed, na.rm=T)[2]))
                         numberOfColors.ratioToInp.x.trnposed <- colorRampPalette(c("red", "black"))(round(range(ratioToInp.x.trnposed, na.rm=T)[2]))
                    }
                    #numberOfColors.x.trnposed <- colorRampPalette(c("black", "white"))(round(range(x.trnposed, na.rm=T)[2]))
                    numberOfColors.x.trnposed <- colorRampPalette(c("white", "black"))(round(range(x.trnposed, na.rm=T)[2]))
                    numberOfColors.ratioToInp.x.trnposed <- colorRampPalette(c("white", "black"))(round(range(ratioToInp.x.trnposed, na.rm=T)[2]))
                    
                    # JPEG plots because pdf sizes are too huge
                    ## Raw heatmap
                    
                    # set breaks depending on range of counts
#                     if(grepl("H3K4", gsub("_R1.*", "", i))){
#                          print("hit: H3K4")
#                          breaks=seq(0, 1500, by=20) #75 values
#                     }
#                     else if(grepl("H3K27", gsub("_R1.*", "", i))){
#                          print("hit: H3K27")
#                          breaks=seq(0, 250, by=2) #125 values
#                     }
#                     else if(grepl("DNMT1", gsub("_R1.*", "", i))){
#                          print("hit: DNMT1")
#                          breaks=seq(0, 50, by=.5) #100 values
#                     }
#                     else if(grepl("EZH2", gsub("_R1.*", "", i))){
#                          print("hit: EZH2")
#                          breaks=seq(0, 120, by=1.5) #80 values
#                     }
#                     else if(grepl("Input", gsub("_R1.*", "", i))){
#                          print("hit: Input")
#                          breaks=seq(0, 50, by=1) #50 values
#                     }
#                     else if(grepl("H3_", gsub("R1.*", "", i))){ #H3
#                          print("hit: H3")
#                          breaks=seq(0, 50, by=.5) #100 values
#                     }
                    h3k4.increment <- h3k4.max/100
                    h3k27.increment <- h3k27.max/100
                    dnmt.increment <- dnmt.max/100
                    ezh2.increment <- ezh2.max/100
                    inp.increment <- inp.max/100
                    h3.increment <- h3.max/100
                    
                    if(grepl("H3K4", gsub("_R1.*", "", i))){
                         print("hit: H3K4")
                         breaks=seq(0, h3k4.max, by=h3k4.increment) #75 values
                    }
                    else if(grepl("H3K27", gsub("_R1.*", "", i))){
                         print("hit: H3K27")
                         breaks=seq(0, h3k27.max, by=h3k27.increment) #125 values
                    }
                    else if(grepl("DNMT1", gsub("_R1.*", "", i))){
                         print("hit: DNMT1")
                         breaks=seq(0, dnmt.max, by=dnmt.increment) #100 values
                    }
                    else if(grepl("EZH2", gsub("_R1.*", "", i))){
                         print("hit: EZH2")
                         breaks=seq(0, ezh2.max, by=ezh2.increment) #80 values
                    }
                    else if(grepl("Input", gsub("_R1.*", "", i))){
                         print("hit: Input")
                         breaks=seq(0, inp.max, by=inp.increment) #50 values
                    }
                    else if(grepl("H3_", gsub("R1.*", "", i))){ #H3
                         print("hit: H3")
                         breaks=seq(0, h3.max, by=h3.increment) #100 values
                    }
                    
                    mycol <- colorpanel(n=length(breaks)-1,low="gray100",mid="gray50",high="gray0")
                    print("raw heatmap")
                    jpeg(heat_plot_raw, height = 500, width = 800, quality=100)
                    heatmap.3(x.trnposed, Rowv=T, Colv=F, col=mycol, RowSideColors=custom.annotation, scale="none", trace="none", dendrogram="row", breaks=breaks, cexRow=1, cexCol=1, key=T, main=paste0("Heatmap of raw seq reads - ", gsub("_R1.*", "", i)), na.rm=TRUE, symkey=FALSE)
                    #         heatmap.3(x.trnposed, Rowv=T, Colv=F, col=numberOfColors.x.trnposed, RowSideColors=custom.annotation, scale="none", trace="none", dendrogram="row", cexRow=0.2, main=paste0("Heatmap of raw seq reads - ", gsub("_R1.*", "", i)), na.rm=TRUE)
                    dev.off()
                    
                    
                    ## Ratio heatmap
                    print("ratio heatmap")
                    jpeg(heat_plot_RATIO, height = 500, width = 800, quality=100)
                    heatmap.3(ratioToInp.x.trnposed, Rowv=T, Colv=F, scale="none", col=numberOfColors.ratioToInp.x.trnposed, trace="none", RowSideColors=custom.annotation, dendrogram="row", cexRow=0.2, main=paste0(main="Heatmap of ratios\nof seq reads to average input - ", gsub("_R1.*", "", i)), na.rm=TRUE, symkey=FALSE)
                    #         heatmap.3(ratioToInp.x.trnposed, Rowv=T, Colv=F, col=mycol, RowSideColors=custom.annotation, scale="none", trace="none", dendrogram="row", breaks=breaks, cexRow=1, cexCol=1, key=T, main=paste0("Heatmap of ratios\nof seq reads to average input - ", gsub("_R1.*", "", i)), na.rm=TRUE)
                    dev.off()
                    
               } else {
                    # These vectors with NA created as a placeholder in the lists for the plot functions below to work smoothly. This is generated if coverage data for at least two genes is not present.
                    chipCoverage.TSSAverageList[[g]] <- rep(NA, num.bins)
                    inputCoverage.TSSAverageList[[g]] <- rep(NA, num.bins)
                    
                    print(paste(names(genelist.info)[g], "coverage data for atleast two genes not present", sep="-"))
                    #Plot dummy heat_plot describing why heatmap is not being plotted
                    heat_plot <- paste(names(genelist.info)[g], "heat_plot.pdf", sep="-")
                    heat_plot <- paste(i, heat_plot, sep="-")
                    heat_plot <- paste(plot.Directory, heat_plot, sep="/")
                    pdf(heat_plot, height = 5, width = 8)
                    plot(0~0)
                    text(x=1, y=0, labels=c("coverage data for atleast two genes not present"))
                    dev.off()
               }
          } # end of for loop to loop through gene lists in genelist.info 
          
          # Plot average coverage and RATIO to input plots
          avAndRATIO_plot <- paste(i, "avAndRATIO_plot.pdf", sep="-")
          avAndRATIO_plot <- paste(plot.Directory, avAndRATIO_plot, sep="/")
          avAndRATIO_plot
          pdf(avAndRATIO_plot, height = 5, width = 8)
          par(mfcol=c(1,1), mar=c(5,5,3,1))
          # Average plot
          plot(0~0, type="n", xlim=range(int), ylim=c(0,range(chipCoverage.TSSAverageList, na.rm=T)[2]), xlab="Dist. from TSS", ylab="Seq Counts", cex.lab=2, cex.axis=2, main="Average Plot")
          plot.col <- vector()
          plot.lty <- vector()
          for(u in 1:length(chipCoverage.TSSAverageList)){
               # Start loop average plot
               plot.col[u] <- u
               if(!any(chipCoverage.TSSAverageList[[u]] %in% NA)){
                    lines(int, chipCoverage.TSSAverageList[[u]], lty=1, lwd=4, col=u)
                    plot.lty[u] <- 1
               } else{lines(int, rep(0, num.bins), lty=3, lwd=4, col=u)
                      text(x=int[1], y=0, labels=c("coverage data for atleast two genes not present"), adj=0, cex=0.8)
                      plot.lty[u] <- 3}
          } # end plotting average plot
          
          # RATIO plot
          ## Get ylim for the RATIO plot
          ratioYLim = matrix(nrow=length(chipCoverage.TSSAverageList), ncol=2)
          for(u in 1:length(chipCoverage.TSSAverageList)){
               ratioYLim[u,] = range(chipCoverage.TSSAverageList[[u]]/inputCoverage.TSSAverageList[[u]])
          }
          ratioYLim = c(0, range(ratioYLim)[2])
          plot(0~0, type="n", xlim=range(int), ylim=ratioYLim, xlab="Dist. from TSS", ylab="Ratio to Input\nof Seq. counts", cex.lab=1, cex.axis=2, main="Ratio to Input Plot")
          plot.col <- vector()
          plot.lty <- vector()
          
          for(u in 1:length(chipCoverage.TSSAverageList)){
               # Start loop ratio plot
               plot.col[u] <- u
               if(!any(chipCoverage.TSSAverageList[[u]] %in% NA)){
                    ratio = chipCoverage.TSSAverageList[[u]]/inputCoverage.TSSAverageList[[u]]
                    lines(int, ratio, lty=1, lwd=4, col=u)
                    plot.lty[u] <- 1
               } else{lines(int, rep(0, num.bins), lty=3, lwd=4, col=u)
                      text(x=int[1], y=0, labels=c("coverage data for atleast two genes not present"), adj=0, cex=0.8)
                      plot.lty[u] <- 3}
          } # end plotting ratio plot
          
          # Plot dummy plot showing the figure legend
          par(mfcol=c(1,1), mar=c(5,5,3,1))
          plot(0~0, type="n", xlim=range(int), ylim=c(0,range(chipCoverage.TSSAverageList, na.rm=T)[2]), xlab="Dist. from TSS", ylab="Seq Counts", cex.lab=2, cex.axis=2)
          legend("topright", legend=names(genelist.info), col=plot.col, lty=plot.lty, lwd=4, cex=1)
          dev.off()
          
          # Write the average coverage values as a table in the same folder where plots are saved
          output.table <- do.call(cbind, chipCoverage.TSSAverageList)
          output.table.name <- paste("average_values", "txt", sep=".")
          output.table.name <- paste(i, output.table.name, sep="-")
          output.table.name <- paste(plot.Directory, output.table.name, sep="/")
          write.table(output.table, output.table.name, sep="\t", quote=F, row.names=F)
     }
}


#######################################################
# Function to plot average enrichment profile as a ratio to the average of input seq read profile
#######################################################
# files.path <- plot.Directory
fun.Plot_RatioToInput <- function(files.path, RegAroundTSS, bin, plot.0to20.81to100){
     # files.path defines the folder/file names of average coverage data as *.txt files. Usually it is the same as plot.Directory value in function fun.average_heat.plots;  RegAroundTSS is the region around TSS used to create the average coverage data; bin is the windows in which seq tags are binned in the coverage data used to create the average coverage data around gene TSS; plot.0to20.81to100 is T or F indicating if only the 0to20 and 81to100 expression groups be plotted or not (this is done because plotting all 5 expression percentile groups with the 10bp binned data jumbles the plot.  So only the top and bottom groups selected).
     
     
     # Define num.bins and int for plotting (this is taken from (fun.average_heat.plots): 
     ## num.bins is the number of coverage bins around the TSS that will be plotted given RegAroundTSS and bin.
     ## If RegAroundTSS/bin is even, add 1 to num.bins to make it odd. This is done because the coverage is plotted for equal number of bins from a central bin. So the total number of bins has to be odd. (%% is the modulo function)
     ## int is the position of bins around the TSS in which coverage is being plotted.
     ## Note that RegAroundTSS should be a multiple of bin (see definition of bin).
     
     if((RegAroundTSS/bin) %% 2 == 0){
          num.bins <- (RegAroundTSS/bin)+1
          int <- seq(from=-(RegAroundTSS/(bin*2)), to=(RegAroundTSS/(bin*2)), by=1)*bin
     } else {
          num.bins <- (RegAroundTSS/bin)
          int <-  seq(from=-((num.bins-1)/2), to= (num.bins-1)/2, by=1)*bin
     }
     length(int)
     
     write.table(int, file=paste0(plot.Directory, "/int.txt"), row.names=FALSE, col.names=FALSE, sep="\t")
     
     # Are the target *.txt files in one folder? 
     ## If TRUE, TargetFilesInOneFolder=T, if FALSE TargetFilesInOneFolder=F
     ## This is done because in certain cases, like the CpG island plots coverage data, the target files are in separate folders within the files.path.
     ifelse(sum(grep(".txt", dir(files.path))) >= 1, TargetFilesInOneFolder <- T, TargetFilesInOneFolder <- F)
     
     ifelse(TargetFilesInOneFolder == T, folders <- basename(files.path), folders <- dir(files.path))
     
     # Plot ratio of average histone-mark enrichment wrt Input enrichment in windows of bin size
     for(i in folders){
          if(TargetFilesInOneFolder == T){#get all average data file names
               fls <- basename(list.files(path=file.path(dirname(files.path), i), ".txt$", full=TRUE))} else {fls <- basename(list.files(path=file.path(files.path, i), ".txt$", full=TRUE))
               }
          INP <- grep("INP", fls, ignore.case=T, value=T) #get file name for input average data
          fls.FOI <- setdiff(fls, INP) #isolate files of interest (all file names except input)
          sapply(fls.FOI, FUN=function(file.x){
               if(TargetFilesInOneFolder == T){
                    av <- file.path(dirname(files.path), i, file.x)
                    INP <- file.path(dirname(files.path), i, INP)
               } else {
                    av <- file.path(files.path, i, file.x)
                    INP <- file.path(files.path, i, INP)
               }
               av <- read.table(av, header=T, sep="\t")
               INP <- read.table(INP, header=T, sep="\t") #INP is the common input file for all the histone-mark average data in this sapply interation
               if(exists("plot.0to20.81to100")){
                    if (plot.0to20.81to100 == T){
                         av <- av[, grep(paste(c("0to20", "81to100"), collapse="|"), colnames(av), value=T)]
                         INP <- INP[, grep(paste(c("0to20", "81to100"), collapse="|"), colnames(INP), value=T)]
                    }
               }
               
               ratio <- av/INP
               
               # Ratio plot
               if(TargetFilesInOneFolder == T){plot.Directory <- file.path(dirname(files.path), i)} else {plot.Directory <- file.path(files.path, i)}
               av_RATIO_plot <- paste(file.x, "av_RATIO_plot.pdf", sep="-")
               av_RATIO_plot <- file.path(plot.Directory, av_RATIO_plot)
               pdf(av_RATIO_plot, height = 10, width = 8)
               par(mfcol=c(2,1), mar=c(5,5,3,1))
               plot(0~0, type="n", xlim=range(int), ylim=c(0,max(ratio)), xlab="Dist. from TSS", ylab="Ratio (to input) of Seq Counts", cex.lab=2, cex.axis=2)
               plot.col <- vector()
               for(u in 1:ncol(ratio)){
                    plot.col[u] <- u
                    lines(int, ratio[,u], lty=1, lwd=4, col=u)
               }
               title(main="Ratio to Input")
               # Plot legend separately (below the above plot)
               plot(0~0, type="n", xlim=range(int), ylim=c(0,max(ratio)), xlab="Dist. from TSS", ylab="Seq Counts", cex.lab=2, cex.axis=2)
               legend("topright", legend=names(ratio), col=plot.col, lty=1, lwd=4, cex=1)
               dev.off()
          })
     }
     
}
#Example: fun.Plot_RatioToInput(files.path="/Users/harieaswaran/Work Related Files/Cluster scripts/Subodh-ChIP_seq/Plots/TSS_Plots_200bpBin_NCCIT-chromatin-genes_10000bp-new", RegAroundTSS=10000, bin=200, plot.0to20.81to100=F)


############################################################################
#Function to retrieve ChIP-seq coverage in regions flanking a genomic element (GE; eg. CpG island)
fun.GetGECoverage <- function(GE.info, flank.coverage.data, central.coverage.data, FlankRegion=2000, bin=200, norm.bins=10){
     #This function extracts the coverage data in a defined region flanking a genomic element (GE; eg. CpG island) (FlankRegion); GE.info is the dataframe consisting the required information (Chromosome, Start, End) about the GE; flank.coverage.data is the coverage data in bed format obtained for the flanking regions using coverageBed in BedTools; central.coverage.data is the coverage data in bed format for the central region of the GE (normalized into bins defined by norm.bins); bin is the window size in which the coverage data is binned (default is 200 bp); norm.bins is the number of bins into which the actual GE (central region) is divided to normalize all GE into a similar size group; FlankRegion is the region flanking up- and down-stream from the edges of the genomic element for which the coverage in the defined bin size has to be extracted (default FlankRegion is 2000bp, i.e. 2000 bp up- and down-stream from the edges of the GE). Output is an object of class list containing two data frames as elements of the list, holding the coverage data in the flanking region upstream and downstream of the GE.
     
     num.bins <- (FlankRegion/bin) #this is the number of coverage bins on each side flanking the genomic element that should be obtained unless the genomic element is at the edge of a chromosome
     
     # split genomic element (GE, eg. CpG island) information by chromosomes
     # split coverage data by chromosomes
     # then get coverage in each chromosome around each GE
     chr <- paste("chr", c(1:22, "X", "Y"), sep="")
     GE.info.by.chr <- lapply(chr, FUN=function(u){GE.info[GE.info[,"Chromosome"] %in% u, ]}) # splits GE.info by chromosomes and creates a list (each dataframe object in list holds information about GE.info in only one chromosome)
     names(GE.info.by.chr) <- chr
     flank.coverage.data <- lapply(chr, FUN=function(u){flank.coverage.data[flank.coverage.data[,"Chromosome"] %in% u, ]}) # splits coverage data by chromosomes and creates a list (each dataframe object in list holds coverage data for only one chromosome)
     names(flank.coverage.data) <- chr
     
     central.coverage.data <- lapply(chr, FUN=function(u){central.coverage.data[central.coverage.data[,"Chromosome"] %in% u, ]}) # splits coverage data by chromosomes and creates a list (each dataframe object in list holds coverage data for only one chromosome)
     names(central.coverage.data) <- chr
     
     newlist <- lapply(chr, FUN=function(u){
          GE.info.by.chr.u <- GE.info.by.chr[[u]]
          
          # Get coverage data only if there is a genomic element (eg. CpG islands) in a particular iteration of a chromosome. 'if' condition used here so that the script does not stop here with an error if there is no genomic element in the query belonging to a particular chromosome.
          if(dim(GE.info.by.chr.u)[1] > 0){
               flank.coverage.data.u <- flank.coverage.data[[u]]
               flank.coverage.data.u <- flank.coverage.data.u[!duplicated(flank.coverage.data.u$start), ] # Some of the binned coverage data is duplicated because if the same flanking regions is part of two or more different genomic elements, i.e. the GE overlap partially, then the coverage data will be obtained for each instance of the flanking regions. Remove these duplicates as it causes problems with "if(length(GE.coverage.up) != num.bins & length(GE.coverage.down) != num.bins)".
               central.coverage.data.u <- central.coverage.data[[u]]
               central.coverage.data.u <- central.coverage.data.u[!duplicated(central.coverage.data.u$start), ] # remove duplicates as above, just in case.
               
               cov.chr <- lapply(1:nrow(GE.info.by.chr.u), FUN=function(x){
                    GE.info.start <- GE.info.by.chr.u[x, "Start"]
                    GE.info.end <- GE.info.by.chr.u[x, "End"]
                    up <- GE.info.start - FlankRegion
                    down <- GE.info.end + FlankRegion #pick the flanking region up and downstream from the edges of the genomic element being queried (eg. CpG island) (defined by the FlankRegion parameter).  The nucleotide position of the edge will be considered as the starting  bin (if bin > 1; if bin=1, i.e. nucleotide-level resolution, then the term bin here would mean nucelotides) and coverage in the required number of bins (given by FlankRegion/bin) up- and down-stream from this starting bin will be extracted
                    ## Since the coverage for the same 200 bp regions could be obtained in regions where the GE are within the FlankRegion (eg. 2000 bp), for. example the coverage will be obtained on 200-400, 400-600 as well as 199-399, 399-599...., it is better to exactly define the windows in which the coverage was obtained and then extract the coverage data by matching to the windows. Hence the windows.up and windows.down are generated here.
                    chr.no <- GE.info.by.chr.u[x, "Chromosome"]
                    # Create the windows for the up flank
                    a1 <- seq(from=up, to=GE.info.start,  by=bin)
                    b1 <- seq(from=up + bin, to=GE.info.start,  by=bin)
                    u1=min(length(a1), length(b1))
                    windows.up <- data.frame("Chromosome"=rep(chr.no, u1), "Start"=a1[1:u1], "End"=b1[1:u1]) #colnames defined so that rbind works below
                    # Create the windows for the down flank
                    a2 <- seq(from=GE.info.end, to=down,  by=bin)
                    b2 <- seq(from=GE.info.end + bin, to=down,  by=bin)
                    u2=min(length(a2), length(b2))
                    windows.down <- data.frame("Chromosome"=rep(chr.no, u2), "Start"=a2[1:u2], "End"=b2[1:u2])
                    GE.coverage.up <- flank.coverage.data.u[flank.coverage.data.u$start %in% windows.up$Start & flank.coverage.data.u$end %in% windows.up$End , "Seq_tags"] #gets the coverage in the upstream flanking regions
                    GE.coverage.down <- flank.coverage.data.u[flank.coverage.data.u$start %in% windows.down$Start & flank.coverage.data.u$end %in% windows.down$End, "Seq_tags"] #gets the coverage in the downstream flanking regions
                    GE.coverage.central <- central.coverage.data.u[central.coverage.data.u$start >= GE.info.start &  central.coverage.data.u$end <= GE.info.end, "Seq_tags"] #gets the coverage in the central region
                    # If the queried genomic element is at an end of a chromosome (or the data for that region is lacking), then report these as NA's so that the plotting functions later on work. This has happened with a central region (Enhancer element chr16:14340722-14343849). Don't know why! 
                    if(any(length(GE.coverage.up) != num.bins | length(GE.coverage.down) != num.bins | length(GE.coverage.central) != norm.bins)) {
                         GE.coverage.up=rep(NA,num.bins)
                         GE.coverage.down=rep(NA,num.bins)
                         GE.coverage.central=rep(NA, norm.bins)
                    }
                    
                    return(list("GE.coverage.up"=GE.coverage.up, "GE.coverage.central"=GE.coverage.central, "GE.coverage.down"=GE.coverage.down))
               })
               cov.chr.upGE <- sapply(1:length(cov.chr), FUN=function(i){cov.chr[[i]]["GE.coverage.up"]}) #upGE for upstream queried genomic element
               cov.chr.downGE <- sapply(1:length(cov.chr), FUN=function(i){cov.chr[[i]]["GE.coverage.down"]}) #downGE for downstream queried genomic element
               cov.chr.centralGE <- sapply(1:length(cov.chr), FUN=function(i){cov.chr[[i]]["GE.coverage.central"]}) #downGE for downstream queried genomic element
               cov.chr <- list("cov.chr.upGE"=do.call(cbind, cov.chr.upGE), "cov.chr.centralGE"=do.call(cbind, cov.chr.centralGE), "cov.chr.downGE"=do.call(cbind, cov.chr.downGE)) #creates a list of three dataframes - holding the coverage data for the upstream flank, normalized central region and downstream flank of each queried genomic element.
               return(cov.chr)
          }
     })
     names(newlist) <- chr #the oblects in this list named for convenience, which is used in the following sapply function
     coverage.up <- sapply(chr, FUN=function(i){newlist[[i]]["cov.chr.upGE"][[1]]})
     coverage.up <- do.call(cbind, coverage.up)
     coverage.down <- sapply(chr, FUN=function(i){newlist[[i]]["cov.chr.downGE"][[1]]})
     coverage.down <- do.call(cbind, coverage.down)
     coverage.central <- sapply(chr, FUN=function(i){newlist[[i]]["cov.chr.centralGE"][[1]]})
     coverage.central <- do.call(cbind, coverage.central)
     
     return(list("coverage.up"=coverage.up, "coverage.central"=coverage.central, "coverage.down"=coverage.down))
}


# Create windows file to generate the coverage data for the flanking region and the exact region of overlap (central region)  
fun.WindowsGE_CentralFlank <- function(GE.info, GE.info.name, GE.info.path, norm.bins=10, ouput.Directory, flank.bin.size=200, flank.region.size=2000, ...){
     # This function creates windows file to create the coverage data for the flanking region and the exact region of overlap (central region)  
     ## GE.info is a dataframe containing the required info (with column names "Chromosome", "Start", "End") for the genomic elements of interest; GE.info.name is a character name describing the genomic elements in GE.info (eg. "CpG-island", "H3K4me3_peaks", etc.); GE.info.path is the path of the directory in which the GE.info dataframe is stored. If this dataframe is rather passed directly as an R object, then select GE.info.path to be NULL ; flank.region.size is the region flanking up- and down-stream from the edges of the genomic element for which the coverage in the defined bin size has to be extracted; flank.bin.size is the window size in which the coverage data in the flanks should be binned; ouput.Directory is a character vector assigning the directory in which the windows output table should be stored; norm.bins (normalization bins), the number of bins in which the coverage data is to be computed for a defined genomic region; ... arguments to be passed to fun.read.table.path (read.table)
     options("scipen"=100, "digits"=4) #scipen set to arbitrary 100 so that scientific notations (like 1.8 e+06) are not printed in the output file below. If this happens further functions in BedTools do not work.
     
     ## 0) Create the output directory (if it does not already exist) where windows file will be stored.  
     dir.create(ouput.Directory)
     
     ## 1) Create windows of defined size for each genomic element such that the window sizes divides the genomic element into a defined number of windows (norm.bins).  This is done to compute coverage data (similar to fun.windows, but for a defined genomic element).  The goal is to compute coverage data in these window sizes such that the number of bins is constant for all genomic elements and then standardize it to a size of 200 bp, so that all the genomic elements can be put on the same plot.
     if(!(is.null(GE.info.path))){
          GE <- fun.read.table.path(file=GE.info, path=GE.info.path, ...)
     } else{
          GE <- GE.info
     }
     ### Remove rows that have strange overlap values (has "i" and "j" values in the Start and End columns respectively)
     #GE <- GE[-grep("i", GE$Start),]
     GE <- transform(GE, Start=as.numeric(Start), End=as.numeric(End))
     GE <- na.omit(GE)
     GE.dif <- GE$End - GE$Start
     GE <- cbind.data.frame(GE, "Overlap.Size" = GE.dif) # name Overlap.Size does not have any real meaning. It was used initially for the peak overlaps; basically a reference to size of the genomic elements.
     head(GE)
     #plot(density(GE.dif), xlim=c(0,5000))
     
     ### Percent peaks that are greater than or lesser than 2000 bp
     #dim(GE[GE.dif <= 2000, ])[1]/nrow(GE)
     #dim(GE[GE.dif >= 2000, ])[1]/nrow(GE)
     
     c1 <- list()
     for(i in 1:nrow(GE)){
          size.bin <- floor(GE$Overlap.Size[i]/norm.bins) # size in bp of each normalized bin
          chr.no <- GE$Chromosome[i]
          a <- seq(from=GE$Start[i], to=GE$End[i],  by=size.bin)
          b <- seq(from=GE$Start[i]+size.bin, to=GE$End[i],  by=size.bin)
          u=min(length(a), length(b))
          c1[[i]] <- data.frame(rep(chr.no, u), a[1:u], b[1:u])
     }
     windows_GE_central <- do.call(rbind, c1)
     colnames(windows_GE_central) <- c("Chromosome", "Start", "End")
     windows_GE_central_file.name <- file.path(ouput.Directory, paste(GE.info.name, "windows_GE_central.bed", sep="-"))
     #return(windows_GE_central)
     write.table(windows_GE_central, file=windows_GE_central_file.name, sep="\t", quote=F, row.names=F, col.names=F)
     
     ##2) Create windows of a defined size (flank.bin.size, usually 200 bp) flanking each genomic element extending up to a region of defined size (flank.region.size)
     c1 <- list()
     for(i in 1:nrow(GE)){
          up.flank <- GE$Start[i] - flank.region.size
          down.flank <- GE$End[i] + flank.region.size
          chr.no <- GE$Chromosome[i]
          # Create the windows for the up flank
          a1 <- seq(from=up.flank, to=GE$Start[i],  by=flank.bin.size)
          b1 <- seq(from=up.flank + flank.bin.size, to=GE$Start[i],  by=flank.bin.size)
          u1=min(length(a1), length(b1))
          A <- data.frame("Chromosome"=rep(chr.no, u1), "Start"=a1[1:u1], "End"=b1[1:u1]) #colnames defined so that rbind works below
          # Create the windows for the down flank
          a2 <- seq(from=GE$End[i], to=down.flank,  by=flank.bin.size)
          b2 <- seq(from=GE$End[i] + flank.bin.size, to=down.flank,  by=flank.bin.size)
          u2=min(length(a2), length(b2))
          B <- data.frame("Chromosome"=rep(chr.no, u2), "Start"=a2[1:u2], "End"=b2[1:u2])
          # rbind A and B
          c1[[i]] <- rbind.data.frame(A, B)
     }
     windows_GE_flank <- do.call(rbind, c1)
     #return(windows_GE_flank)
     windows_GE_flank_file.name <- file.path(ouput.Directory, paste(GE.info.name, "windows_GE_flank.bed", sep="-"))
     write.table(windows_GE_flank, file=windows_GE_flank_file.name, sep="\t", quote=F, row.names=F, col.names=F)
}
#Example: fun.WindowsGE_CentralFlank(GE.info="K4K27_peak_overlap.bed", GE.info.name="K4K27_peak_overlap", GE.info.path="/Users/harieaswaran/Work Related Files/Cluster scripts/Subodh-ChIP_seq/Plots/Sicer_calls/Sicer_data_new", norm.bins=10, ouput.Directory="/Users/harieaswaran/Work Related Files/Cluster scripts/Subodh-ChIP_seq/Plots/Sicer_calls/Sicer_data_new", flank.bin.size=200, flank.region.size=5000, header=T, sep="\t", stringsAsFactors=F)


#Function for plotting coverage data flanking the genomic elements of interest as average values and heatmap
###change flank to FlankRegion
fun.GE_average_heat.plots <- function(GE.info, coverage_files, FlankRegion, bin, GE="CpGi", PlotHeatmap=T, plot.Directory){
     # Plots the coverage data flanking the genomic elements of interest as average plot and heatmap; GE.info is a list of dataframes containing the required info (with column names "Chromosome", "Start", "End") for the genomic elements of interest (the objects in this list should have a name so that the plots can take these names); coverage_files is a vector of file names of the ChIP-seq coverage data; FlankRegion is the region flanking up- and down-stream from the edges of the genomic element for which the coverage in the defined bin size has to be extracted; bin is the window size in which the coverage data is binned; GE is character vector defining the nature of the genomic element around which the coverage is being plotted (eg. "CpGi" for CpG island); PlotHeatmap controls if heatmaps should be plotted (T or F); plot.Directory is a character vector assigning the directory in which the plots should be stored.
     
     library(gplots, lib.loc=lib.path)
     num.bins <- (FlankRegion/bin)+1 #this is the number of coverage bins on each side flanking the genomic element that should be obtained unless the genomic element is at the edge of a chromosome
     int <- c(seq(from=-(FlankRegion/bin), to=0, by=1)*bin, seq(from=0, to=(FlankRegion/bin), by=1)*bin)#int is the position of bins flanking the GE in which coverage is being plotted.
     int[num.bins:(num.bins+1)] <- GE
     
     for(i in coverage_files){
          ############################################################
          ##Load the ChIP_seq coverage data
          ############################################################
          #Read the bed file containing the coverage in bins defined by the bin size using the coverageBed command in BEDtools.
          chip.coverage <- read.table(i, header=F, sep="\t")
          colnames(chip.coverage) <- c("Chromosome", "start", "end", "Seq_tags")
          GE.info.x.up.av <- list() #create list to store the average coverage values upstream
          GE.info.x.down.av <- list() #create list to store the average coverage values downstream
          for(g in 1:length(GE.info)){
               #Retrieve coverage data for regions flanking each GE for each mark
               x <- fun.GetGECoverage(GE.info=GE.info[[g]], coverage.data=chip.coverage, FlankRegion=FlankRegion, bin=bin)
               x.up <- x[["coverage.up"]]
               x.down <- x[["coverage.down"]]
               if(dim(t(na.omit(t(x.up))))[2] >= 2 & dim(t(na.omit(t(x.down))))[2] >= 2){#there should be coverage data for at least two GE on both flanking sides of the GE. na.omit removes rows that contain NA in a data frame. t used here because columns containing NAs have to be removed - thus after omitting columns with NA, coverage data should have data for atleast 2 genes to satisfy the 'if' statement.
                    x.average.up <- apply(x.up, 1, mean, na.rm=T)
                    x.average.down <- apply(x.down, 1, mean, na.rm=T)
                    GE.info.x.up.av[[g]] <- x.average.up
                    GE.info.x.down.av[[g]] <- x.average.down
                    names(GE.info.x.up.av)[g] <- paste(i, names(GE.info)[g], sep="-")
                    names(GE.info.x.down.av)[g] <- paste(i, names(GE.info)[g], sep="-")        
                    
                    if(PlotHeatmap==T){#Plot coverage as heatmap plot
                         heat_plot <- paste(names(GE.info)[g], "heat_plot.pdf", sep="-")
                         heat_plot <- paste(i, heat_plot, sep="-")
                         heat_plot <- paste(plot.Directory, heat_plot, sep="/")
                         pdf(heat_plot, height = 5, width = 8)
                         x.trnposed <- cbind(t(x.up), t(x.down))
                         colnames(x.trnposed) <- int
                         rownames(x.trnposed) <- c(1:nrow(x.trnposed)) #rownames changed because otherwise by default it gets the name "GE.coverage.up"
                         heatmap.2(x.trnposed, Rowv=F, Colv=F, scale="none", col=greenred(50), trace="none", dendrogram="none")
                         dev.off()
                    }
               } else{
                    GE.info.x.up.av[[g]] <- rep(NA, num.bins) #this vector created as a placeholder in the list for the plot functions below to work smoothly
                    GE.info.x.down.av[[g]] <- rep(NA, num.bins) #this vector created as a placeholder in the list for the plot functions below to work smoothly
                    print(paste(names(GE.info)[g], "coverage data for atleast two genes not present", sep="-"))
                    if(PlotHeatmap==T){
                         #Plot dummy heat_plot describing why heatmap is not being plotted
                         heat_plot <- paste(names(GE.info)[g], "heat_plot.pdf", sep="-")
                         heat_plot <- paste(i, heat_plot, sep="-")
                         heat_plot <- paste(plot.Directory, heat_plot, sep="/")
                         pdf(heat_plot, height = 5, width = 8)
                         plot(0~0)
                         text(x=1, y=0, labels=c("coverage data for atleast two GE not present"))
                         dev.off()
                    }
               }
          }
          
          #Plot average coverage plot
          av_plot <- paste(i, "av_plot.pdf", sep="-")
          av_plot <- paste(plot.Directory, av_plot, sep="/")
          pdf(av_plot, height = 5, width = 8)
          par(mfcol=c(1,1), mar=c(5,5,3,1))
          par(las=2)
          plot(rep(0, num.bins*2), xaxt = "n", type="n", ylim=c(0,range(GE.info.x.up.av, GE.info.x.down.av, na.rm=T)[2]), xlab="Regions flanking", ylab="Seq Counts", cex.lab=2, cex.axis=2) # Null plot made using 0s just to make the plotting area
          plot.col <- vector()
          plot.lty = vector()
          for(i in 1:length(GE.info.x.up.av)){
               plot.col[i] <- i
               if(!any(GE.info.x.up.av[[i]] %in% NA) & !any(GE.info.x.down.av[[i]] %in% NA)){
                    y=c(GE.info.x.up.av[[i]], GE.info.x.down.av[[i]])
                    lines(y, type="b", lty=1, lwd=4, col=i)
                    polygon(c(num.bins, num.bins, num.bins+1, num.bins+1), c(0, max(y)*1.5, max(y)*1.5, 0), col=rgb(red=0.7, green=0.7, blue=0.7, alpha=0.6, maxColorValue=1))
                    axis(1, at=1:length(int), labels=int)
                    plot.lty[i] <- 1
               } else{lines(rep(0, num.bins*2), lty=3, lwd=4, col=i)
                      text(x=1, y=0, labels=c("coverage data for atleast two GE not present"), adj=0, cex=0.8)
                      plot.lty[i] <- 3
               }
          }
          legend("topright", legend=names(GE.info), col=plot.col, lty=plot.lty, lwd=4, cex=1)
          dev.off()
     }
}
#fun.GE_average_heat.plots(GE.info=CpGi.list, coverage_files=c("coverage_chr22.bed"), FlankRegion=2000, bin=200, GE="CpGi", PlotHeatmap=F)


#Function for QC plotting of duplicate reads
QCplot.fun <- function(chromInfo, chromInfo.file_path, bam.files, bam.files_path, axis.labels){
     #plots the duplicate reads as ratio and proportion of unique or totla read respectively. chromInfo is the character vector describing the filename that contains the genome information (chromosomes and their sizes) downloaded from USCSC; chromInfo.file_path is the path to chromInfo; bam.files is a character vector of bam files; bam.files_path is the path to the directory in which the bam files and the bam index files (.bai) are stored; axis.labels is character vector of abbreviated names of bam.files, in the SAME ORDER as bam.files, for labeling the X-axis of the plots.
     library(Rsamtools)
     
     # Read and setup the chromInfo table that contains information about chromosomes
     chromInfo <- fun.read.table.path(file=chromInfo, path=chromInfo.file_path)
     chromInfo <- chromInfo[-grep("_", chromInfo[,1]), ] #remove rows that contain characters of the form "chr1_random".
     colnames(chromInfo) <- c("Chromosome", "Chromosome_Length", "source")
     chromInfo <- subset(chromInfo, chromInfo$Chromosome != "chrM") #remove mitochondrial chromosome
     rownames(chromInfo) <- chromInfo$Chromosome
     
     #Create a list of ratio and proportion of duplicated sequence tags for each cromosome in each bam file
     ratio <- vector("list", nrow(chromInfo)) #list to collect the ratio of duplicated seq tags to unique seq tags
     proportion <- vector("list", nrow(chromInfo)) #list to collect the proportion of duplicated seq tags to total number of seq tags (duplicated + unique tags)
     pos.frequency <- vector("list", length(bam.files)) #list oblect to collect the frequency of occurrence of position in each chromosome. The frequency from all chromosomes will be appended to a vector for each ChIP-seq, hence the length of this list is equal to the number of ChIP-seq samples/bam files
     for(a in 1:nrow(chromInfo)){#for each chromosome (a), get the ratio and proportion of duplicated sequence tags for each bam file
          which <- RangesList(IRanges(1, chromInfo[rownames(chromInfo)[a], "Chromosome_Length"])) #defines which chromosome and what region is to be obtained
          names(which) <- rownames(chromInfo)[a]
          #scanBamWhat() #this list the fields that can be obtained from the bam file using the what command
          what <- c("rname", "strand", "pos", "qwidth") #defines what fields from the bam file should be obtained
          param <- ScanBamParam(which=which, what=what) #creates instance of class param to be interpreted by scanBam
          for(i in 1:length(bam.files)){#for each bam file (i), compute the ratio and proportion of duplicated sequence tags
               bamFile <- paste(bam.files_path, bam.files[i], sep="")
               bam <- scanBam(bamFile, param=param)
               #Collapse bam, which is a list-of-lists, into a single list
               ## Source: http://www.biostat.jhsph.edu/~lcollado/Presentations/SNB_BiocHTS.R
               lst <- lapply(names(bam[[1]]), function(elt) {
                    do.call(c, unname(lapply(bam, "[[", elt)))
               })
               names(lst) <- names(bam[[1]])
               #Convert lst to dataframe object
               bam <- do.call("DataFrame", lst)
               unique.start <- sum(!duplicated(bam$pos)) #number of uniqe seq tags
               duplicated.start <- sum(duplicated(bam$pos)) # #number of duplicated seq tags
               ratio[[a]][i] <- duplicated.start/unique.start
               proportion[[a]][i] <- duplicated.start/(unique.start + duplicated.start)
               pos.frequency[[i]] <- c(pos.frequency[[i]], as.vector(table(bam$pos))) #counts the number of times sequence reads are repeated in each chromsome for each bam file. The frequencies are collected as an object for each bam file (i).
          }
          #names(ratio)[a] <- rownames(chromInfo)[a]
          names(proportion)[a] <- rownames(chromInfo)[a]
     }
     ratio <- do.call(rbind, ratio)
     proportion <- do.call(rbind, proportion)
     colnames(ratio) <- axis.labels
     colnames(proportion) <- axis.labels  
     #Plot the ratio and proportion stats
     pdf("QC_duplicates.pdf", height = 8, width = 5)
     par(mfcol=c(2,1), mar=c(5,5,3,1))
     boxplot(ratio)
     title(main="Distribution of Ratio (duplicate reads:unique reads)\nacross all chromosomes", sub="ChIP-seq samples")
     boxplot(proportion)
     title(main="Distribution of Proportion (duplicates reads to total reads)\nacross all chromosomes", sub="ChIP-seq samples")
     dev.off()
     #Plot the histogram of frequency at which reads are repeated (plots frequency reads with 1, 2, 3.... number of repeats; for clarity, different inervals of read repetitions are plotted)
     for(i in 1:length(axis.labels)){
          freq_hist <- paste("QC_freq_hist", i, sep="-")
          freq_hist <- paste(freq_hist, "pdf", sep=".")
          pdf(freq_hist, height = 5, width = 8)
          par(mfcol=c(2,2), mar=c(5,5,3,1))
          hist(pos.frequency[[i]], breaks=c(0:max(pos.frequency[[i]])), main="", col="black")
          title(main=paste("Frequency of seq. reads repetition\n", axis.labels[i], sep=""), sub="Number of times reads are repeated")
          #Plot histogram of seq reads repeated 0:10 times
{rep.range=c(0:10)
 x <- pos.frequency[[i]][pos.frequency[[i]] %in% rep.range] #retain read positions that occur only 0:10 times
 hist(x, breaks=rep.range, main="", col="black")
 title(main=paste("Frequency of seq. reads repetition\n", axis.labels[i], sep=""), sub="Number of times reads are repeated (range 0:10)")}
#Plot histogram of seq reads repeated 0:25 times
{rep.range=c(0:25)
 x <- pos.frequency[[i]][pos.frequency[[i]] %in% rep.range] #retain read positions that occur only 0:25 times
 hist(x, breaks=rep.range, main="", col="black")
 title(main=paste("Frequency of seq. reads repetition\n", axis.labels[i], sep=""), sub="Number of times reads are repeated (range 0:25)")}
#Plot histogram of seq reads repeated 0:50 times
{rep.range=c(0:50)
 x <- pos.frequency[[i]][pos.frequency[[i]] %in% rep.range] #retain read positions that occur only 0:50 times
 hist(x, breaks=rep.range, main="", col="black")
 title(main=paste("Frequency of seq. reads repetition\n", axis.labels[i], sep=""), sub="Number of times reads are repeated (range 0:50)")}
dev.off()
     }
}
#QCplot.fun(chromInfo="hg18_chromInfo.txt", chromInfo.file_path="/amber2/scratch/baylin/Hari/JHUSB01002/JHUSB01002_000_analysis/ExtraFiles_for_Analysis", bam.files=c("JHUSB01002_002_K4_SS02.sorted.bam", "JHUSB01002_003_K27_SS03.sorted.bam", "JHUSB01002_004_2o_SS04.sorted.bam", "JHUSB01002_005_H2AZ_SS05.sorted.bam", "JHUSB01002_001_INP_SS01.sorted.bam"), bam.files_path="/amber2/scratch/baylin/Hari/JHUSB01002/JHUSB01002_000_analysis/sorted_bam/", axis.labels=c("K4", "K27", "Biv", "H2AZ", "Input"))



#######################################################################################
# Function for identifying genes with peaks within 5000bp/spanning TSS of genes
## This is modification of the functions (fun.getpeakswithin5000bpTSS.1 and fun.getpeakswithin5000bpTSS.2) used in the MSC/Osteo Genome Research paper. This picks peaks whose:
### 1) Start OR End is within 5000 bp from TSS or 
### 2) the Start AND End is greater than 5000 bp from TSS (thus will pick peaks that are very broad and span large regions across the TSS)
## The modification here is that:
### 1) fun.getpeakswithin5000bpTSS.1 and fun.getpeakswithin5000bpTSS.2 were combined to have a single function called fun.getGenesWithPeaksWithinTSS
### 2) The regions around the TSS (eg. 5000bp) can be passed as an argument, regionaroundTSS, in the function.
## Output is the table containing genes and the info
######################################################################################
fun.getGenesWithPeaksWithinTSS <- function(genelist.info, dat.name, dat.path, dat.workspace=NULL, SICER.control.lib, regionaroundTSS, output.Directory, save.output.dataframe, ...){
     # This function generates a dataframe object that contains the gene names and associated number of peaks (from SICER peak call output) around its TSS, along with the TSS, Chromosome; genelist.info is the dataframe consisting the required information about the genes (typically info for all genes obtained from biomaRt); dat.name is name of the SICER peak call output file that contains the peak call information; dat.path is the path to the peak call file; dat.workspace is NULL if the peak data is being read from file on disk, or T if it is being referenced from the workspace; SICER.control.lib is TRUE or FALSE values indicating if a control library was used for peak calling; regionaroundTSS is the length in basepairs upstream and downstream (for eg. 2500 bp upstream and 2500 bp downstream) in which presence /absence of peaks will be tested; output.Directory is the directory in which the histogram plot of peak calls and dataframe containing the gene names and peak numbers has to be saved;  save.output.dataframe T or F value indicating if dataframe containing the gene names and peak numbers be saved or not; ... arguments to be passed to fun.read.table.path (read.table)
     
     # Load the SICER peak call data
     # Depending on the SICER peak call data had a control library (input) or not, the column names differ. Hence this has to be defined by SICER.control.lib
     if (is.null(dat.workspace)){# whether or not the peak data file should be read from file o disk or from R workspace
          dat <- fun.read.table.path(file=dat.name, path=dat.path, ...)
          if(SICER.control.lib==T){
               colnames(dat) <- c("Chromosome", "Start", "End", "ChIP_island_read_count", "CONTROL_island_read_count", "p_value", "fold_change", "FDR_threshold")
               dat <- cbind.data.frame(dat, "Length"=dat[,"End"] - dat[,"Start"])
          } else {
               colnames(dat) <- c("Chromosome", "Start", "End", "score")
               dat <- cbind.data.frame(dat, "Length"=dat[,"End"] - dat[,"Start"])
          }
     } else { # else peak data is used from R workspace
          dat <- dat.name # set the appropriate column names in the referenced dat.name
     }
     
     # Plot a histogram of the peak sizes
     hist_plot <- paste(dat.name, "peak_call_hist.pdf", sep="-")
     hist_plot <- paste(output.Directory, hist_plot, sep="/")
     pdf(hist_plot, height = 5, width = 8)
     par(mfcol=c(1,1), mar=c(5,5,3,1))
     hist(log2(dat[,"Length"]), xlab="Peak Lengths (Log2)", col="grey40", ylab="Frequency", cex.lab=2, cex.axis=2, main=paste("Histogram of peak sizes for", dat.name, sep=" "))
     dev.off()
     
     # Create chromosome name objects
     chr <- c("chr1", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr2", "chr20", "chr21", "chr22", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chrX", "chrY")
     
     # The chromsome names in the peak call dataframe and the genelist.info are respectively of the form "chr1" and "1".  To keep them consistent, change the chromosome names in the dat object to the form in genelist.info.
     for(i in 1:length(chr)) {
          dat[dat[,1] %in% chr[i], "Chromosome"] <- chr[i]
     }
     
     X <- lapply(chr, FUN=function(v){# this function outputs a dataframe containing gene names, TSS,  chromosome name and info about peak calls within the regionaroundTSS for genes that have at least  one peak call within the regionaroundTSS.
          # Get the rows from genelist.info and dat that match a particular chromosome (v).
          genelist.chr <- genelist.info[genelist.info[,"Chromosome"] %in% v, ]
          dat.chr <- dat[dat[,"Chromosome"] %in% v, ]
          peaks.table <- sapply(c(1:nrow(genelist.chr)), FUN=function(u){# this function was formerly fun.getpeakswithin5000bpTSS.1
               TSS <- genelist.chr[u, "TSS"]
               # Set the  boundary for the broad region around the TSS within which peaks will be analyzed for conditions 1 and 2 (below). This is obtained as the largest peak size + regionaroundTSS to ensure that any peak and even the largest peak will bordering the  regionaroundTSS will be considered. By taking the largest peak size, automatically every peak that can potentially overlap the regionaroundTSS will be considered
               bound <- max(dat.chr[,"Length"]) + regionaroundTSS
               # Get the peak calls that +/- from the TSS the number of base pairs defined by "bound" (boundary). This is done so that the subsequent steps to identify if a peak is inclusive in a defined region around the TSS (regionaroundTSS) are executed faster.
               
               ## Various configuration of largest peak wrt TSS +/- bound that will be considered
               ## Largest peak= ~~~~~~~~~~~~~~~~
               
               ##        bound          -regionaroundTSS   +regionaroundTSS          bound
               ##        ~~~~~~~~~~~~~~~~«««««««««««««««TSS»»»»»»»»»»»»»»»~~~~~~~~~~~~~~~~
               ##        ~~~~~~~~~~~~~~~~
               ##             ~~~~~~~~~~~~~~~~
               ##                                                   ~~~~~~~~~~~~~~~~
               ##                                                         ~~~~~~~~~~~~~~~~
               
               Sub <- dat.chr[dat.chr[,"Start"] >= (TSS - bound) & dat.chr[,"End"] <= (TSS + bound), c("Start", "End")]
               Peaks.Info <- vector("list", length=nrow(Sub)) # for grabing genelist info for the current TSS iteration and Peak coordinates
               for(i in 1:length(Sub[,1])){
                    a <- Sub[i,"Start"]; b <- Sub[i,"End"] #Start and End postion of i-th peak in Sub
                    # The following conditions for peaks overlapping the regionaroundTSS does not take into account the strand on which the gene is present because it does not matter.  We only want the peaks that are in the vicinity of the TSS, any asymmetric nature of the peaks are not considered:
                    # Condition-1 for peaks overalapping regionaroundTSS: (a)Pick peaks whose Start OR End are within the length defined by regionaroundTSS. If the Start OR End of the i-th peak in Sub is upstream or downstream from the TSS within the length defined by regionaroundTSS, return 1, else NA.  This will return 1 for all peak calls satisfying the condtion along with the peak information.
                    
                    ##                       -regionaroundTSS   +regionaroundTSS
                    ##        ________________«««««««««««««««TSS»»»»»»»»»»»»»»»________________
                    ##                               a~~~~b
                    ##                a~~~~~~~~~~~b
                    ##        a~~~~~~~~~~~~~~~b
                    ##                                a~~~~~~~~~~~~~~b
                    ##                                                       a~~~~~~b
                    ##                               a~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~b
                    
                    ifelse(abs(TSS-a) <= regionaroundTSS | abs(TSS-b) <= regionaroundTSS, Peaks.Info[[i]] <- cbind(genelist.chr[u,], Sub[i,]) , Peaks.Info[[i]] <- NA)
                    
                    # Condition-2 for peaks overalapping regionaroundTSS: (a) Picks peaks whose Start AND End are at a distance greater than the length defined by regionaroundTSS. If the Start AND End of the i-th peak in Sub is, respectively, upstream and downstream from the TSS by greater than the length defined by regionaroundTSS, return 1, else NA.
                    
                    ##                       -regionaroundTSS   regionaroundTSS
                    ##        ________________«««««««««««««««TSS»»»»»»»»»»»»»»»________________
                    ##                a~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~b
                    
                    ifelse(TSS-a >= regionaroundTSS & TSS-b <= -regionaroundTSS, Peaks.Info[[i]] <- cbind(genelist.chr[u,], Sub[i,]) , Peaks.Info[[i]] <- Peaks.Info[[i]])
               }
               # Get total number of peaks associated with a gene TSS
               Peaks.Info <- na.omit(do.call(rbind.data.frame, Peaks.Info))
               return(Peaks.Info)
          }
          ) # this creates the peaks.table for every gene u in genelist.chr
          
          # Combine peaks.table for all genes in genelist.chr
          peaks.table <- do.call(rbind, peaks.table)
          return(peaks.table)
     })
     X <- do.call(rbind, X)
     if (save.output.dataframe==T){
          output.name <- paste("GenesWithPeaks", dat.name, sep="-") #give a name to the output file
          write.table(X, file=paste(output.Directory, output.name, sep="/"), quote=F, row.names=F, sep="\t")
     }
     return(X)
}


# Function to get the regions of peak overlaps from two different peak call files
fun.getPeakOverlap <- function(peak.dat1.name, peak.dat2.name, peak.dat1.path, peak.dat2.path, ...){
     #This function gets the regions of peak overlaps from two different peak call files; peak.dat1.name and peak.dat2.name are the two peak call files; peak.dat1.path and peak.dat2.path are the directories in which the two peak call files are located; ... arguments to be passed to fun.read.table.path (read.table)
     # Get regions of overlap
     ## Get peak call data for each chromosome
     ## For a peak region in peak.dat1, get the coordinates that exactly overlap with the peaks in peak.dat2.
     peak.dat1 <- fun.read.table.path(file=peak.dat1.name, path=peak.dat1.path, ...)
     peak.dat2 <- fun.read.table.path(file=peak.dat2.name, path=peak.dat2.path, ...)
     colnames(peak.dat1) <- c("Chromosome", "Start", "End", "ChIP_island_read_count", "CONTROL_island_read_count", "p_value", "fold_change", "FDR_threshold")
     colnames(peak.dat2) <- c("Chromosome", "Start", "End", "ChIP_island_read_count", "CONTROL_island_read_count", "p_value", "fold_change", "FDR_threshold")
     peak.dat1 <- transform(peak.dat1, Start=as.numeric(Start), End=as.numeric(End))
     peak.dat2 <- transform(peak.dat2, Start=as.numeric(Start), End=as.numeric(End))
     
     # Get the chromosome names in the peak call data so as to get overlaps only in those chromosomes present in the two peak call data sets.
     chr <- unique(c(peak.dat1[,"Chromosome"], peak.dat2[,"Chromosome"]))
     
     # AB is an example of peak call in peak.dat1
     # Various ab are possible peak positions of every peak in peak.dat2
     ##       _____________________A===========B__________________________
     ## cond.1                        a~~~~b
     ## cond.2                a~~~~~~~b
     ## cond.3                             a~~~~~~~b
     ## cond.4                 a~~~~~~~~~~~~~~~~~~~~~~~~b
     ## cond.5                     a~~~~~~~~~~~b
     ## cond.6                     a~~~~~~~~~~~~~~~~~~b
     ## cond.7             a~~~~~~~~~~~~~~~~~~~b
     ## cond.8       a~~~~~~~b
     ## cond.9                                             a~~~~~~~b
     
     cond.1 <- c(-1, -1, 1, 1) #overlap = ab
     cond.2 <- c(1, -1, 1, 1) #overlap = Ab
     cond.3 <- c(-1, -1, 1, -1) #overlap = aB
     cond.4 <- c(1, -1, 1, -1) #overlap = AB
     cond.5 <- c(NaN, -1, 1, NaN) #overlap = AB=ab
     cond.6 <- c(NaN, -1, 1, -1) #overlap = AB
     cond.7 <- c(1, -1, 1, NaN) #overlap = AB
     cond.8 <- c(1, 1, 1, 1) # No overlap
     cond.9 <- c(-1, -1, -1, -1) # No overlap
     
     overlap <- lapply(chr, FUN=function(v){
          peak.dat1.chr <- peak.dat1[peak.dat1$Chromosome %in% v, ]
          peak.dat2.chr <- peak.dat2[peak.dat2$Chromosome %in% v, ]
          # For every peak in peak.dat1.chr (rows I), get all the peaks in peak.dat2.chr (rows J)
          overlap.IJ <- lapply(1:nrow(peak.dat1.chr), FUN=function(i){#
               A <- peak.dat1.chr[i, "Start"]; B <- peak.dat1.chr[i, "End"]
               overlap.iJ <- lapply(1:nrow(peak.dat2.chr), FUN=function(j){##
                    a <- peak.dat2.chr[j, "Start"]; b <- peak.dat2.chr[j, "End"]
                    distances <- c((A-a)/abs(A-a), (A-b)/abs(A-b), (B-a)/abs(B-a), (B-b)/abs(B-b)) # divided by absolute values to get values as 1 or -1
                    if(identical(distances, cond.8) | identical(distances, cond.9)) {overlap.ij.x <- c(v,0,0)} else {
                         if(identical(distances, cond.1)) {overlap.ij.x <- c(v,a,b)} else {
                              if(identical(distances, cond.2)) {overlap.ij.x <- c(v,A,b)} else {
                                   if(identical(distances, cond.3)) {overlap.ij.x <- c(v,a,B)} else {
                                        if(identical(distances, cond.4)) {overlap.ij.x <- c(v,A,B)} else {
                                             if (identical(distances, cond.5) | identical(distances, cond.6) | identical(distances, cond.7)) {overlap.ij.x <- c(v,A,B)} else {
                                                  overlap.ij.x <- c(v, paste("i", i, sep="="), paste("j", j, sep="=")) # this gets the rows in the two dataframes that for some reason do not match any of the conditions
                                             }
                                        }
                                   }
                              }
                         }
                    }
                    return(overlap.ij.x) # overlaps between i-th row in peak.dat1.chr and j-th row in peak.dat2.chr
               })##
               overlap.iJ <- do.call(rbind, overlap.iJ)
               return(overlap.iJ) # all overlaps between i-th row in peak.dat1.chr and all j rows in peak.dat2.chr for chromosome v
          })#
          overlap.IJ <- data.frame(do.call(rbind, overlap.IJ)) # all overlaps between all i rows in peak.dat1.chr and all j rows in peak.dat2.chr for chromosome v
          colnames(overlap.IJ) <- c("Chromosome", "Start", "End")
          overlap.IJ <- overlap.IJ[!(overlap.IJ$Start == 0  & overlap.IJ$End == 0), ]
          return(overlap.IJ)
     })
     overlap <- data.frame(do.call(rbind, overlap)) # all overlaps between all i rows in peak.dat1.chr and all j rows in peak.dat2.chr for all chromosomes
     colnames(overlap) <- c("Chromosome", "Start", "End")
     return(overlap)
}


# Function to plot the distribution of the distances of centre of regions of interest, like peak calls or CpG-islands, relative to the TSS
fun.plot.ROIsDistributionAroundTSS <- function(ROI, genelist.info, regionaroundTSS, plot.name, plot.Directory){
     #ROIs are the regions of interest, like peaks calls, CpG islands, etc., the centre of whose distribution has to be plotted around the TSS. The ROIs is a named list object containing the list of dataframes that has the minimal information about the regions ("Chromosome", "Start", "End"). This function will plot the central position of the regions by calcualting the centre of the region from its Start and End positions; genelist.info is the dataframe consisting the required information ("Gene", "Chromosome", "TSS", "Strand") about the genes (typically obtained from biomaRt; regionaroundTSS is the region up/downstream from TSS within which the ROI centre shold be to be considered for the plot; plot.name is the name of the output plot file;  plot.Directory is the directory in which plot should be saved)
     
     # Approach: For each TSS in genelist.info:
     ## get the centre of ROIs that are within -regionaroundTSS to +regionaroundTSS from the TSS (note that the +/-2500bp from TSS was used to assign peaks to the corresponding TSS but using  regionaroundTSS=2500 in this function is a bit problematic because, for e.g. for the K4-exclusive genes, only the K4-peak centres will be within +/-2500bp from TSS and it is very unlikely that (K27-peak or K4K27-overlap centres) will fall in that range. Hence, there won't be any of the latter peaks identified to fall in that window and so the plot will ultimately have only the distribution of K4-peak centres. Thus, in general it makes sense to plot a window a much larger than +/-2500bp from TSS because there is a higher chance of picking distrbution of other peaks in that window).
     ## compute the distance of the centre of each ROI from TSS to which it belongs.
     ## plot the density of these distances which will give the probabaility distribution of the distances of the ROIs/peaks/overlaps from the TSS.
     
     # Split genelist.info by chromosomes and create a list (each dataframe object in list holds information about genes in one chromosome)
     chr <- paste("chr", c(1:22, "X", "Y"), sep="")
     #genelist.info$Chromosome <- paste("chr", genelist.info$Chromosome, sep="") # this has to be  done because the chromaosmes are of the form 1,2,3... in the output from Biomart
     genelist.info <- transform(genelist.info, TSS=as.numeric(TSS)) # convert TSS to numeric
     genelist.info <- lapply(chr, FUN=function(u){genelist.info[genelist.info[,"Chromosome"] %in% u, ]}) # splits genelist.info by chromosomes
     names(genelist.info) <- chr
     
     ROI.PeakSE.Dist.From.TSS <- lapply(ROI, FUN=function(ROI.x){
          # If Chromosome in ROI are of form "1", change them to the form "chr1". Just to make sure that it is consistent with genelist.info (This happens due to creation of these files in various ways).
          chr.num=c(1:22, "X", "Y")
          for(w in 1:length(chr.num)){ROI.x[ROI.x$Chromosome %in% chr.num[w], "Chromosome"] <- chr[w]}
          # Convert the Start and End values to numeric (these columns could be integer or character depending on how it was created, loaded, etc).
          ROI.x <- transform(ROI.x, Start=as.numeric(Start), End=as.numeric(End))
          
          # Compute the centre position of the ROIs
          ROI.x$centre.ROI <- apply(ROI.x[,c("Start", "End")], 1, mean)
          # Split ROI.x by chromosomes
          ROI.x <- lapply(chr, FUN=function(u){ROI.x[ROI.x[,"Chromosome"] %in% u, ]})
          names(ROI.x) <- chr
          
          # Compute the distance from TSS of the centre position of the ROIs and get teh Start and End of corresponding peaks
          PeakSE.Dist.From.TSS <- lapply(chr, FUN=function(u){
               genelist.info.u <- genelist.info[[u]]
               ROI.x.u <- ROI.x[[u]]
               if(dim(genelist.info.u)[1] > 0){# if there are any genes in this iteration of chromosome (otherwise gives an error)
                    dist.from.TSS <- vector("list", length=nrow(genelist.info.u)) #for holding deistance of peak-centre from TSS
                    Peak.Start.End <- vector("list", length=nrow(genelist.info.u)) #for holdng peak start and end
                    for(i in 1:nrow(genelist.info.u)){if(genelist.info.u$Strand[i] == 1) {
                         a <- which(abs(ROI.x.u$centre.ROI-genelist.info.u$TSS[i]) <= regionaroundTSS)
                         dist.from.TSS[[i]] <- ROI.x.u[a, "centre.ROI"]-genelist.info.u$TSS[i]
                         Peak.Start.End[[i]] <- ROI.x.u[a, c("Start", "End")] - genelist.info.u$TSS[i]
                    } else {
                         a <- which(abs(ROI.x.u$centre.ROI-genelist.info.u$TSS[i]) <= regionaroundTSS)
                         dist.from.TSS[[i]] <- genelist.info.u$TSS[i]-ROI.x.u[a, "centre.ROI"]
                         Peak.Start.End[[i]] <- genelist.info.u$TSS[i]-ROI.x.u[a, c("Start", "End")]
                    }
                    }
                    dist.from.TSS <- unlist(dist.from.TSS)
                    Peak.Start.End <- do.call(rbind, Peak.Start.End)
                    #return(dist.from.TSS)
                    return(list(dist.from.TSS=dist.from.TSS, Peak.Start.End=Peak.Start.End))
               }
          })
          Dist.From.TSS <- unlist(lapply(c(1:length(PeakSE.Dist.From.TSS)), FUN=function(p){return(PeakSE.Dist.From.TSS[[p]]$dist.from.TSS)}))
          PeakSE <- lapply(c(1:length(PeakSE.Dist.From.TSS)), FUN=function(p){return(PeakSE.Dist.From.TSS[[p]]$Peak.Start.End)})
          PeakSE <- do.call(rbind, PeakSE)
          return(list(Dist.From.TSS=Dist.From.TSS, PeakSE=PeakSE))
     })
     # For each ROI, create separate lists with the Dist.From.TSS and PeakSE
     ROI.Dist.From.TSS <- lapply(ROI.PeakSE.Dist.From.TSS, FUN=function(L){return(L$Dist.From.TSS)})
     names(ROI.Dist.From.TSS) <- names(ROI)
     ROI.PeakSE <- lapply(ROI.PeakSE.Dist.From.TSS, FUN=function(L){return(L$PeakSE)})
     names(ROI.PeakSE) <- names(ROI)
     
     # Plot density distributon of the distances of the ROIs from the TSS
     plot.col <- c(1:length(ROI.Dist.From.TSS))
     plot.lty <- c(1:length(ROI.Dist.From.TSS))
     plot.name.density <- file.path(plot.Directory,  paste(plot.name, "ROI_DensityDistr_wrt_TSS.pdf", sep="_"))
     ## Get the maximum of the density distributions for the ylim
     yy <- lapply(ROI.Dist.From.TSS, FUN=function(x){
          if(length(x) > 2){(density(x)$y)*length(x)}# need at least 2 points to select a bandwidth automatically
     })
     #xx <- lapply(ROI.Dist.From.TSS, FUN=function(x){(density(x)$x)})
     y.lim=max(range(yy))
     #x.lim=range(xx)
     x.lim= (regionaroundTSS + 500)*c(-1,1) #set the x-limit to the range where we are interested (+ an additional 500bp so taht teh regionaroundTSS boundary shows in the plot)
     pdf(plot.name.density, height = 5, width = 16)
     par(mfcol=c(1,2))
     ## Plot-1: Density plot
     plot(0~0, type="n", main="Density distributon of ROI-centre distances from TSS", ylim=c(0, y.lim), xlim=x.lim, cex.axis=1.5, axes=F, xlab="Distance of peak-center to TSS (Kilobases)", ylab="Frequency")
     legend("topright", legend=names(ROI.Dist.From.TSS), col=plot.col, lty=plot.lty, lwd=4, cex=1)
     axis(side=1, at=seq(-regionaroundTSS, regionaroundTSS, by=regionaroundTSS/10), labels=seq(-regionaroundTSS/1000, regionaroundTSS/1000, by=regionaroundTSS/10000), pos=0, lty=1, cex.axis=1.5)
     axis(2)
     for(i in 1:length(ROI.Dist.From.TSS)){
          if(length(ROI.Dist.From.TSS[[i]] > 2)){
               x.points <- density(ROI.Dist.From.TSS[[i]])$x # distances from TSS
               y.points <- density(ROI.Dist.From.TSS[[i]])$y*length(ROI.Dist.From.TSS[[i]]) #multiplied by the number of observations so that the scale of the plot gives a relative frequency among peaks/ROIs for the different histone marks.  Just plotting the probablilities (density distribution) masks the relative frequencies.
               lines(y.points~x.points, col=plot.col[i], lty=plot.lty[i], lwd=4, ylim=c(0, y.lim))
          }
     }
     
     ## Plot-2: segment plot
     x.lim=(regionaroundTSS*3 + 500)*c(-1,1)
     plot(0~0, type="n", xlim=x.lim, ylim=c(0,4), main="Segment plot of peak length distribution", axes=F, xlab="Distance from TSS (Kilobases)")
     axis(side=1, at=seq(-regionaroundTSS*3, regionaroundTSS*3, by=regionaroundTSS*3/10), labels=seq(-regionaroundTSS*3/1000, regionaroundTSS*3/1000, by=regionaroundTSS*3/10000), pos=0, lty=1, cex.axis=1.5)
     axis(2)
     for(i in 1:length(ROI.PeakSE)){
          d <- ROI.PeakSE[[i]]
          if(length(nrow(d) > 1)){ # i.e. if there are any peaks in this iteration
               for(q in 1:nrow(d)){segments(x0=d[q,"Start"], y0=i, x1=d[q,"End"], y1=i, lwd=20, col=rgb(t(col2rgb(plot.col[i], alpha=F)), alpha=1, maxColorValue=255))}
          }
     }
     
     dev.off()
     
     # Plot histogram of frequency distributon of the distances of the ROIs from the TSS
     plot.name.hist <- file.path(plot.Directory,  paste(plot.name, "ROI_HistDistr_wrt_TSS.pdf", sep="_"))
     pdf(plot.name.hist, height = 15, width = 8)
     par(mfcol=c(3,1))
     for(i in 1:length(ROI.Dist.From.TSS)){
          if(length(ROI.Dist.From.TSS[[i]] > 2)){
               hist(ROI.Dist.From.TSS[[i]], freq=T, main=names(ROI.Dist.From.TSS)[i], col=i, cex.axis=1.5)}
     }
     dev.off()
}


# Function for loading the minimal information about ROIs.
fun.getROI <- function(get_CpG.Island=F){
     #This function loads the the minimal information ("Chromosome", "Start", "End") about the ROIs (regions of interest), like peaks calls, CpG islands, etc., whose distribution has to be plotted around the TSS. The output is used in fun.plot.ROIsDistributionAroundTSS.
     #get_CpG.Island is a T or F indicating if the information for CpG-islands should be loaded. By defaultit is F, and the function loads the SICER H3K4Me3 and H3K27Me3-peak calls, and the regions of exact overlaps of H3K4Me3 and H3K27Me3. Check the function for the files that are loaded (incorporate the following variables as function arguments}.
     # Read the H3K4Me3 and H3K2Me3 peaks call data
     peak.dat1.name="K4-W200-G600-islands-summary-FDR0.01"
     peak.dat2.name="K27-W200-G1000-islands-summary-FDR0.01"
     peak.dat1.path="/amber2/scratch/baylin/Hari/JHUSB01002/JHUSB01002_000_analysis/SICER_calls/K4"
     peak.dat2.path="/amber2/scratch/baylin/Hari/JHUSB01002/JHUSB01002_000_analysis/SICER_calls/K27"
     
     peak.dat1 <- fun.read.table.path(file=peak.dat1.name, path=peak.dat1.path, header=F, sep="\t", stringsAsFactors=F)
     peak.dat2 <- fun.read.table.path(file=peak.dat2.name, path=peak.dat2.path, header=F, sep="\t", stringsAsFactors=F)
     colnames(peak.dat1) <- c("Chromosome", "Start", "End", "ChIP_island_read_count", "CONTROL_island_read_count", "p_value", "fold_change", "FDR_threshold")
     colnames(peak.dat2) <- c("Chromosome", "Start", "End", "ChIP_island_read_count", "CONTROL_island_read_count", "p_value", "fold_change", "FDR_threshold")
     peak.dat1 <- transform(peak.dat1, Start=as.numeric(Start), End=as.numeric(End))
     peak.dat2 <- transform(peak.dat2, Start=as.numeric(Start), End=as.numeric(End))
     
     # Read the exact H3K4Me3 and H3K27Me3 peak overlap data (generated using fun.getPeakOverlap in the R script Get_peak_overlaps.R)
     K4K27_peak_overlap <- read.table("/amber2/scratch/baylin/Hari/JHUSB01002/JHUSB01002_000_analysis/SICER_calls/Peak_Overlap/K4K27_peak_overlap.txt", header=T, sep="\t", stringsAsFactors=F)
     ## Remove rows that have strange overlap values (has "i" and "j" values in the Start and End columns respectively)
     K4K27_peak_overlap <- K4K27_peak_overlap[-grep("i", K4K27_peak_overlap[,2]),]
     K4K27_peak_overlap <- transform(K4K27_peak_overlap, Start=as.numeric(Start), End=as.numeric(End))
     
     if(get_CpG.Island == T){
          CpG.I <- read.table("/amber2/scratch/baylin/Hari/JHUSB01002/JHUSB01002_000_analysis/ExtraFiles_for_Analysis/hg18_cpgIslandExt.txt", header=F, sep="\t")
          #CpG.I <- read.table("/Users/harieaswaran/Work Related Files/Cluster scripts/Subodh-ChIP_seq/Required Data/hg18_data_fromUCSC/hg18_cpgIslandExt.txt", header=F, sep="\t")
          colnames(CpG.I) <- c("chrom", "chromStart", "chromEnd",  "name",  "length",  "cpgNum",  "gcNum",  "perCpg",  "perGc",  "obsExp")
          CpG.I <- data.frame(CpG.I[,c("chrom", "chromStart", "chromEnd")])
          colnames(CpG.I) <- c("Chromosome", "Start", "End")
     }
     
     # Create a list of the peak/island data (Regions of interest) to pass to fun.plot.ROIsDistributionAroundTSS
     if(get_CpG.Island == T){
          ROI <- list(peak.dat1, peak.dat2, K4K27_peak_overlap, CpG.I)
          names(ROI) <- c("H3K4Me3_peaks", "H3K27Me3_peaks", "H3K4Me3+H3K27Me3_overlaps", "CpG-islands")
     } else {
          ROI <- list(peak.dat1, peak.dat2, K4K27_peak_overlap)
          names(ROI) <- c("H3K4Me3_peaks", "H3K27Me3_peaks", "H3K4Me3+H3K27Me3_overlaps")
     }
     return(ROI)
}


# Function for getting enhancers that overlap completely or partially with CpG-islands, and thhose that don't.
fun.getEnhancerOverlappingCpGI <- function(enhancerData, cpgIData){
     # enhancerData is the enhancer dataframe (with minimal columns "Chromosome", "Start", "End", "name"); cpgIData is the CpG islands data (with minimal columns "Chromosome", "Start", "End").
     
     # Split enhancerData by chromosomes
     enhancerDataByChr <- lapply(unique(enhancerData$Chromosome), FUN=function(i){
          x <- enhancerData[enhancerData$Chromosome %in% i, ]
          return(x)
     })
     names(enhancerDataByChr) <- unique(enhancerData$Chromosome)
     
     # Split CpG.I by Chromosome
     cpgIByChr <- lapply(unique(cpgI$Chromosome), FUN=function(i){
          x <- cpgI[cpgI$Chromosome %in% i, ]
          return(x)
     })
     names(cpgIByChr) <- unique(cpgI$Chromosome)
     
     enhancerDataWithCpGI <- vector("list", length(enhancerDataByChr))
     enhancerDataWithoutCpGI <- vector("list", length(enhancerDataByChr))
     
     for(chr in names(enhancerDataByChr))
     {
          # get CpG-I and Enhancers for each chromosome
          cpgI.chr <- cpgIByChr[[chr]]
          enhancerData.chr <- enhancerDataByChr[[chr]]
          
          # which enhancer in this iteration of chromsome has a CpG island within (cond1) or partailly overlapping (left is cond2, right is cond3) or completely encompassing (cond4) the enhancer Start and End
          # Get enhancer rows with CpG islands
          enhancerData.chrWithCpGI <- sapply(1:nrow(enhancerData.chr), FUN=function(e){
               enhStart <- enhancerData.chr$Start[e] # enhancer Start
               enhEnd <- enhancerData.chr$End[e] # enhancer End
               cond1 <- which( (cpgI.chr$Start >=  enhStart & cpgI.chr$End <=  enhEnd) )
               cond2 <- which( (cpgI.chr$Start <=  enhStart & cpgI.chr$End >= enhStart) )
               cond3 <- which( (cpgI.chr$Start <=  enhEnd & cpgI.chr$End >= enhEnd) )
               cond4 <- which( (cpgI.chr$Start <=  enhStart & cpgI.chr$End >= enhEnd) )
               if(any(c(cond1, cond2, cond3, cond4) >= 1) == T){
                    return(e)
               } else{ return(NA) }
          })
          enhancerData.chrWithCpGI <- enhancerData.chrWithCpGI[complete.cases(enhancerData.chrWithCpGI)]
          enhancerData.chrWithoutCpGI <- setdiff(c(1:nrow(enhancerData.chr)), enhancerData.chrWithCpGI) #Enhancer rows without CpG islands.
          
          enhancerDataWithCpGI[[chr]] <- enhancerData.chr[enhancerData.chrWithCpGI, ]
          enhancerDataWithoutCpGI[[chr]] <- enhancerData.chr[enhancerData.chrWithoutCpGI, ]
     }# end for loop looping thrugh chromosomes
     enhancerDataWithCpGI <- do.call(rbind, enhancerDataWithCpGI)
     enhancerDataWithoutCpGI <- do.call(rbind, enhancerDataWithoutCpGI)
     return( list(enhancerDataWithCpGI=enhancerDataWithCpGI, 
                  enhancerDataWithoutCpGI=enhancerDataWithoutCpGI) )
}
# Example
## x <- fun.getEnhancerOverlappingCpGI(enhancerData=sampleE100, cpgIData=cpgI)


# Function for reading Bam files and making QC plots
fun.ReadBam.QC_Plots <- function(bam.files, bam.files_path, chr.to.analyze=NULL, plot.dir){
     #chr.to.analyze is NULL if all chromosomes to be analyzed, else specify the chromosome as a vector, eg. c("chr1", "chr22"); bam.files is vector of bam files to be read; bam.files_path is path to bam files; plot.dir is directory where plots have to be saved.
     ###
     library(ShortRead)
     #library(Rsamtools)
     #library(GenomicFeatures)
     #library(chipseq)
     #library(BayesPeak)
     library(BSgenome.Hsapiens.UCSC.hg18)
     #library(rtracklayer)
     
     # Get the hg18 chromosome information
     chromInfo <- seqlengths(Hsapiens)
     chromInfo <- chromInfo[-grep(paste(c("_", "chrM"), collapse="|"), names(chromInfo))] #remove elements that contain characters of the form "chr1_random" and remove mitochondrial chromosome
     
     # CREATE DIRECTORY FOR QC PLOTS
     dir.create(file.path(plot.dir, "QC_Plots"))
     qc.plot.dir <- file.path(plot.dir, "QC_Plots")
     
     if(!(is.null(chr.to.analyze))) {chr.to.analyze <- chr.to.analyze} else{chr.to.analyze <- names(chromInfo)} #if a specific chromosome is not to be analyzed, then analyze all chromosomes
     # Read data from all bam files and store each as an element in a list object
     cov.islands.CseqData <- seqapply(bam.files, FUN=function(bamf, param){ #Gets the coverage data for a given chromosome
          bamFile <- file.path(bam.files_path, bamf) #set the path to to the bam file that should be read.
          aln <- seqapply(chr.to.analyze, FUN=function(a){ #for each chromosome (a), set the param object to read the corresponding region in the bam file and get the coverage vector
               which <- RangesList(IRanges(1, chromInfo[a])); names(which) <- a #defines what region is to be obtained and which chromosome
               #which <- GRanges(a, IRanges(1, chromInfo[a])) #alternative of defining the chromosome and region in one step (try this once if it works).
               #scanBamWhat() #this lists the fields that can be obtained from the bam file using the what command (check http://samtools.sourceforge.net/SAMv1.pdf for details of field description)
               what <- c("rname", "strand", "pos", "qwidth") #defines what fields from the bam file should be obtained
               param <- ScanBamParam(which=which, what=what) #creates instance of class param to be interpreted by scanBam
               aln <- readGappedAlignments(bamFile, format="BAM", param=param) #Read bam file as GappedAlignments object. See help readGappedAlignments for details; if param is not defined, it runs into memory issues.
               seqlevels(aln) <- names(bamWhich(param)) #seqlevels are adjusted to contain just the levels (chromosome name) we are interested in, rather than all levels in the BAM file (default gives all chromosome names)
               aln <- as(aln, "GRanges") #Coerce to a GRanges object as a convenient way to retrieve minimal amount of data (chromosome name, start, end, strand)
               aln <- resize(aln, width=150) #Extend the reads by length of the nucleosomal DNA to get 150 bp reads.
               return(aln)
          })
          names(aln) <- chr.to.analyze
          
          # calculate the coverage vector (how many times each base in the genome is covered by each of the resized intervals in aln)
          cov.islands <- slice(coverage(aln), lower=1) #get regions consisting of contiguous segments of non-zero coverage for all chromosomes.
          
          # QC Plots
          ## QC Plots-1: PLOT NUMBER OF TIMES SAME TAG IS SEQUENCED
          ### Pending: Compute nonredundant fraction (NRF): the ratio between the number of positions in the genome that uniquely mappable reads map to and the total number of uniquely mappable reads ({Landt, 2012 #1711}).
          start.occurrence <- unlist(sapply(chr.to.analyze, FUN=function(a){
               runLength(Rle(start(aln[[a]]))) # Number of times each start postion of the read occurs in the current interation of chromosome
          }))
          plot.name <- file.path(qc.plot.dir, paste(c("StartOccurrence-", bamf, ".jpg"), collapse=""))
          pdf(plot.name)
          hist(start.occurrence, freq=FALSE, breaks=range(start.occurrence)[2], col=c("green"), ylab="Frequency", xlab="Number of times same tag is sequenced", main="")
          title(main="Frequency distribution of seq reads repetition", sub=bamf)
          dev.off()
          rm(start.occurrence, aln)
          gc()
          
          ## QC Plots-2: PLOT WIDTHS OF CONTIGUOUS SEGMENTS OF NON-ZERO COVERAGE (CALLED ISLANDS)
          x <- table(width(cov.islands)) #get the width of islands; colnames here are the widths and rows are the frequencies of the widths.
          x <- colSums(x)
          X <- as.numeric(names(x)) #X-axis values; converted to numeric so that log option can be used in plotting function below
          Y <- as.vector(x) #Y-axis values
          plot.name <- file.path(qc.plot.dir, paste(c("Island_Widths-", bamf, ".jpg"), collapse=""))
          pdf(plot.name, width=8, height=8)
          #par(ylog=T)
          plot(Y~X, log="xy", ylab="Frequency", xlab="CONTIGUOUS SEGMENTS (bp) OF NON-ZERO COVERAGE", col=rgb(1,0,0, alpha=0.5), bg=rgb(1,0,0, alpha=0.5), cex=0.5, pch=19) 
          title(main="FREQUENCY OF WIDTHS OF CONTIGUOUS SEGMENTS (bp)\nOF NON-ZERO COVERAGE", sub=bamf)
          dev.off()
          rm(x,X,Y)
          gc()
          
          ## QC Plots-3: PLOTS OF READ DEPTH & MAX. COVERAGE AS A FUNCTION OF ISLAND WIDTH
          X <- unlist(width(cov.islands)) #get the width of islands; colnames here are the widths and rows are the frequencies of the widths.
          Y <- unlist(viewSums(cov.islands)) #total number of times all nucleotides in a island are read (depth)
          Z <- unlist(viewMaxs(cov.islands)) #maximum number of reads in a given island
          plot.name <- file.path(qc.plot.dir, paste(c("Depth_maxCov-", bamf, ".jpg"), collapse=""))
          pdf(plot.name, width=8, height=8)
          par(mfcol=c(1,2))
          #plot(Y~X, log="y", ylab="Total coverage (Seq. Depth)", xlab="Island width (bp)", col=rgb(0,1,0, alpha=0.4), pch=19, cex=0.3)
          scatter.smooth(x=X, y=Y, log="y", family="gaussian", ylab="Total coverage (Seq. Depth)", xlab="Island width (bp)", col=rgb(0,1,0, alpha=0.4), pch=19, cex=0.3)
          title(main="READ DEPTH AS A FUNCTION OF ISLAND WIDTH (bp)", sub=bamf)
          #plot(Z~X, ylab="Max. Coverage", xlab="Island width (bp)", col=rgb(0,0,1, alpha=0.4), pch=19, cex=0.3)
          scatter.smooth(x=X, y=Z, log="y", family="gaussian", ylab="Max. Coverage", xlab="Island width (bp)", col=rgb(0,0,1, alpha=0.4), pch=19, cex=0.3)
          title(main="MAX. COVERAGE AS A FUNCTION OF ISLAND WIDTH (bp)", sub=bamf)
          dev.off()
          return(cov.islands)
     }, param=param)
     names(cov.islands.CseqData) <- bam.files
     return(cov.islands.CseqData)
     gc()
     ###
}
# Example: cov.islands.CseqData <- fun.ReadBam.QC_Plots(bam.files=c("JHUSB01002_002_K4_SS02.sorted.bam", "JHUSB01002_003_K27_SS03.sorted.bam", "JHUSB01002_004_2o_SS04.sorted.bam", "JHUSB01002_005_H2AZ_SS05.sorted.bam", "JHUSB01002_001_INP_SS01.sorted.bam"), bam.files_path="/amber2/scratch/baylin/Hari/JHUSB01002/JHUSB01002_000_analysis/sorted_bam", chr.to.analyze=NULL, plot.dir = "/amber2/scratch/baylin/Hari/JHUSB01002/JHUSB01002_000_analysis/BED_Files/new/Plots")


##################################################################
# Create RData object containing above functions
##################################################################
#save.image("/Users/harieaswaran/Work Related Files/Cluster scripts/ChIP-SeqAnalysisScripts/R_Scripts/LibraryOfFunctions_RData/ChIP-SeqLibraryOfFunctions.RData")

heatmap.3 <- function(x,
                      Rowv = TRUE, Colv = if (symm) "Rowv" else TRUE,
                      distfun = dist,
                      hclustfun = hclust,
                      dendrogram = c("both","row", "column", "none"),
                      symm = FALSE,
                      scale = c("none","row", "column"),
                      na.rm = TRUE,
                      revC = identical(Colv,"Rowv"),
                      add.expr,
                      breaks,
                      symbreaks = max(x < 0, na.rm = TRUE) || scale != "none",
                      col = "heat.colors",
                      colsep,
                      rowsep,
                      sepcolor = "white",
                      sepwidth = c(0.05, 0.05),
                      cellnote,
                      notecex = 1,
                      notecol = "cyan",
                      na.color = par("bg"),
                      trace = c("none", "column","row", "both"),
                      tracecol = "cyan",
                      hline = median(breaks),
                      vline = median(breaks),
                      linecol = tracecol,
                      margins = c(5,5),
                      ColSideColors,
                      RowSideColors,
                      side.height.fraction=0.3,
                      cexRow = 0.2 + 1/log10(nr),
                      cexCol = 0.2 + 1/log10(nc),
                      labRow = NULL,
                      labCol = NULL,
                      key = TRUE,
                      keysize = 1.5,
                      density.info = c("none", "histogram", "density"),
                      denscol = tracecol,
                      symkey = max(x < 0, na.rm = TRUE) || symbreaks,
                      densadj = 0.25,
                      main = NULL,
                      xlab = NULL,
                      ylab = NULL,
                      lmat = NULL,
                      lhei = NULL,
                      lwid = NULL,
                      NumColSideColors = 1,
                      NumRowSideColors = 1,
                      KeyValueName="Value",...){
     
     invalid <- function (x) {
          if (missing(x) || is.null(x) || length(x) == 0)
               return(TRUE)
          if (is.list(x))
               return(all(sapply(x, invalid)))
          else if (is.vector(x))
               return(all(is.na(x)))
          else return(FALSE)
     }
     
     x <- as.matrix(x)
     scale01 <- function(x, low = min(x), high = max(x)) {
          x <- (x - low)/(high - low)
          x
     }
     retval <- list()
     scale <- if (symm && missing(scale))
          "none"
     else match.arg(scale)
     dendrogram <- match.arg(dendrogram)
     trace <- match.arg(trace)
     density.info <- match.arg(density.info)
     if (length(col) == 1 && is.character(col))
          col <- get(col, mode = "function")
     if (!missing(breaks) && (scale != "none"))
          warning("Using scale=\"row\" or scale=\"column\" when breaks are",
                  "specified can produce unpredictable results.", "Please consider using only one or the other.")
     if (is.null(Rowv) || is.na(Rowv))
          Rowv <- FALSE
     if (is.null(Colv) || is.na(Colv))
          Colv <- FALSE
     else if (Colv == "Rowv" && !isTRUE(Rowv))
          Colv <- FALSE
     if (length(di <- dim(x)) != 2 || !is.numeric(x))
          stop("`x' must be a numeric matrix")
     nr <- di[1]
     nc <- di[2]
     if (nr <= 1 || nc <= 1)
          stop("`x' must have at least 2 rows and 2 columns")
     if (!is.numeric(margins) || length(margins) != 2)
          stop("`margins' must be a numeric vector of length 2")
     if (missing(cellnote))
          cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
     if (!inherits(Rowv, "dendrogram")) {
          if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in%
                                                            c("both", "row"))) {
               if (is.logical(Colv) && (Colv))
                    dendrogram <- "column"
               else dedrogram <- "none"
               warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
                       dendrogram, "'. Omitting row dendogram.")
          }
     }
     if (!inherits(Colv, "dendrogram")) {
          if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in%
                                                            c("both", "column"))) {
               if (is.logical(Rowv) && (Rowv))
                    dendrogram <- "row"
               else dendrogram <- "none"
               warning("Discrepancy: Colv is FALSE, while dendrogram is `",
                       dendrogram, "'. Omitting column dendogram.")
          }
     }
     if (inherits(Rowv, "dendrogram")) {
          ddr <- Rowv
          rowInd <- order.dendrogram(ddr)
     }
     else if (is.integer(Rowv)) {
          hcr <- hclustfun(distfun(x))
          ddr <- as.dendrogram(hcr)
          ddr <- reorder(ddr, Rowv)
          rowInd <- order.dendrogram(ddr)
          if (nr != length(rowInd))
               stop("row dendrogram ordering gave index of wrong length")
     }
     else if (isTRUE(Rowv)) {
          Rowv <- rowMeans(x, na.rm = na.rm)
          hcr <- hclustfun(distfun(x))
          ddr <- as.dendrogram(hcr)
          ddr <- reorder(ddr, Rowv)
          rowInd <- order.dendrogram(ddr)
          if (nr != length(rowInd))
               stop("row dendrogram ordering gave index of wrong length")
     }
     else {
          rowInd <- nr:1
     }
     if (inherits(Colv, "dendrogram")) {
          ddc <- Colv
          colInd <- order.dendrogram(ddc)
     }
     else if (identical(Colv, "Rowv")) {
          if (nr != nc)
               stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
          if (exists("ddr")) {
               ddc <- ddr
               colInd <- order.dendrogram(ddc)
          }
          else colInd <- rowInd
     }
     else if (is.integer(Colv)) {
          hcc <- hclustfun(distfun(if (symm)
               x
               else t(x)))
          ddc <- as.dendrogram(hcc)
          ddc <- reorder(ddc, Colv)
          colInd <- order.dendrogram(ddc)
          if (nc != length(colInd))
               stop("column dendrogram ordering gave index of wrong length")
     }
     else if (isTRUE(Colv)) {
          Colv <- colMeans(x, na.rm = na.rm)
          hcc <- hclustfun(distfun(if (symm)
               x
               else t(x)))
          ddc <- as.dendrogram(hcc)
          ddc <- reorder(ddc, Colv)
          colInd <- order.dendrogram(ddc)
          if (nc != length(colInd))
               stop("column dendrogram ordering gave index of wrong length")
     }
     else {
          colInd <- 1:nc
     }
     retval$rowInd <- rowInd
     retval$colInd <- colInd
     retval$call <- match.call()
     x <- x[rowInd, colInd]
     x.unscaled <- x
     cellnote <- cellnote[rowInd, colInd]
     if (is.null(labRow))
          labRow <- if (is.null(rownames(x)))
               (1:nr)[rowInd]
     else rownames(x)
     else labRow <- labRow[rowInd]
     if (is.null(labCol))
          labCol <- if (is.null(colnames(x)))
               (1:nc)[colInd]
     else colnames(x)
     else labCol <- labCol[colInd]
     if (scale == "row") {
          retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
          x <- sweep(x, 1, rm)
          retval$rowSDs <- sx <- apply(x, 1, sd, na.rm = na.rm)
          x <- sweep(x, 1, sx, "/")
     }
     else if (scale == "column") {
          retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
          x <- sweep(x, 2, rm)
          retval$colSDs <- sx <- apply(x, 2, sd, na.rm = na.rm)
          x <- sweep(x, 2, sx, "/")
     }
     if (missing(breaks) || is.null(breaks) || length(breaks) < 1) {
          if (missing(col) || is.function(col))
               breaks <- 16
          else breaks <- length(col) + 1
     }
     if (length(breaks) == 1) {
          if (!symbreaks)
               breaks <- seq(min(x, na.rm = na.rm), max(x, na.rm = na.rm),
                             length = breaks)
          else {
               extreme <- max(abs(x), na.rm = TRUE)
               breaks <- seq(-extreme, extreme, length = breaks)
          }
     }
     nbr <- length(breaks)
     ncol <- length(breaks) - 1
     if (class(col) == "function")
          col <- col(ncol)
     min.breaks <- min(breaks)
     max.breaks <- max(breaks)
     x[x < min.breaks] <- min.breaks
     x[x > max.breaks] <- max.breaks
     if (missing(lhei) || is.null(lhei))
          lhei <- c(keysize, 4)
     if (missing(lwid) || is.null(lwid))
          lwid <- c(keysize, 4)
     if (missing(lmat) || is.null(lmat)) {
          lmat <- rbind(4:3, 2:1)
          
          if (!missing(ColSideColors)) {
               #if (!is.matrix(ColSideColors))
               #stop("'ColSideColors' must be a matrix")
               if (!is.character(ColSideColors) || nrow(ColSideColors) != nc)
                    stop("'ColSideColors' must be a matrix of nrow(x) rows")
               lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 1)
               #lhei <- c(lhei[1], 0.2, lhei[2])
               lhei=c(lhei[1], side.height.fraction*NumColSideColors, lhei[2])
          }
          
          if (!missing(RowSideColors)) {
               #if (!is.matrix(RowSideColors))
               #stop("'RowSideColors' must be a matrix")
               if (!is.character(RowSideColors) || ncol(RowSideColors) != nr)
                    stop("'RowSideColors' must be a matrix of ncol(x) columns")
               lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1), lmat[,2] + 1)
               #lwid <- c(lwid[1], 0.2, lwid[2])
               lwid <- c(lwid[1], side.height.fraction*NumRowSideColors, lwid[2])
          }
          lmat[is.na(lmat)] <- 0
     }
     
     if (length(lhei) != nrow(lmat))
          stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
     if (length(lwid) != ncol(lmat))
          stop("lwid must have length = ncol(lmat) =", ncol(lmat))
     op <- par(no.readonly = TRUE)
     on.exit(par(op))
     
     layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
     
     if (!missing(RowSideColors)) {
          if (!is.matrix(RowSideColors)){
               par(mar = c(margins[1], 0, 0, 0.5))
               image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
          } else {
               par(mar = c(margins[1], 0, 0, 0.5))
               rsc = t(RowSideColors[,rowInd, drop=F])
               rsc.colors = matrix()
               rsc.names = names(table(rsc))
               rsc.i = 1
               for (rsc.name in rsc.names) {
                    rsc.colors[rsc.i] = rsc.name
                    rsc[rsc == rsc.name] = rsc.i
                    rsc.i = rsc.i + 1
               }
               rsc = matrix(as.numeric(rsc), nrow = dim(rsc)[1])
               image(t(rsc), col = as.vector(rsc.colors), axes = FALSE)
               if (length(rownames(RowSideColors)) > 0) {
                    #        axis(1, 0:(dim(rsc)[2] - 1)/(dim(rsc)[2] - 1), rownames(RowSideColors), las = 2, tick = FALSE)
                    axis(1, 0:(dim(rsc)[2] - 1)/max(1,(dim(rsc)[2] - 1)), rownames(RowSideColors), las = 2, tick = FALSE)
               }
          }
     }
     
     if (!missing(ColSideColors)) {
          
          if (!is.matrix(ColSideColors)){
               par(mar = c(0.5, 0, 0, margins[2]))
               image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
          } else {
               par(mar = c(0.5, 0, 0, margins[2]))
               csc = ColSideColors[colInd, , drop=F]
               csc.colors = matrix()
               csc.names = names(table(csc))
               csc.i = 1
               for (csc.name in csc.names) {
                    csc.colors[csc.i] = csc.name
                    csc[csc == csc.name] = csc.i
                    csc.i = csc.i + 1
               }
               csc = matrix(as.numeric(csc), nrow = dim(csc)[1])
               image(csc, col = as.vector(csc.colors), axes = FALSE)
               if (length(colnames(ColSideColors)) > 0) {
                    axis(2, 0:(dim(csc)[2] - 1)/max(1,(dim(csc)[2] - 1)), colnames(ColSideColors), las = 2, tick = FALSE)
               }
          }
     }
     
     par(mar = c(margins[1], 0, 0, margins[2]))
     x <- t(x)
     cellnote <- t(cellnote)
     if (revC) {
          iy <- nr:1
          if (exists("ddr"))
               ddr <- rev(ddr)
          x <- x[, iy]
          cellnote <- cellnote[, iy]
     }
     else iy <- 1:nr
     image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, breaks = breaks, ...)
     retval$carpet <- x
     if (exists("ddr"))
          retval$rowDendrogram <- ddr
     if (exists("ddc"))
          retval$colDendrogram <- ddc
     retval$breaks <- breaks
     retval$col <- col
     if (!invalid(na.color) & any(is.na(x))) { # load library(gplots)
          mmat <- ifelse(is.na(x), 1, NA)
          image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "",
                col = na.color, add = TRUE)
     }
     axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0,
          cex.axis = cexCol)
     if (!is.null(xlab))
          mtext(xlab, side = 1, line = margins[1] - 1.25)
     axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0,
          cex.axis = cexRow)
     if (!is.null(ylab))
          mtext(ylab, side = 4, line = margins[2] - 1.25)
     if (!missing(add.expr))
          eval(substitute(add.expr))
     if (!missing(colsep))
          for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0, length(csep)), xright = csep + 0.5 + sepwidth[1], ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
     if (!missing(rowsep))
          for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, col = sepcolor, border = sepcolor)
     min.scale <- min(breaks)
     max.scale <- max(breaks)
     x.scaled <- scale01(t(x), min.scale, max.scale)
     if (trace %in% c("both", "column")) {
          retval$vline <- vline
          vline.vals <- scale01(vline, min.scale, max.scale)
          for (i in colInd) {
               if (!is.null(vline)) {
                    abline(v = i - 0.5 + vline.vals, col = linecol,
                           lty = 2)
               }
               xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
               xv <- c(xv[1], xv)
               yv <- 1:length(xv) - 0.5
               lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
          }
     }
     if (trace %in% c("both", "row")) {
          retval$hline <- hline
          hline.vals <- scale01(hline, min.scale, max.scale)
          for (i in rowInd) {
               if (!is.null(hline)) {
                    abline(h = i + hline, col = linecol, lty = 2)
               }
               yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
               yv <- rev(c(yv[1], yv))
               xv <- length(yv):1 - 0.5
               lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
          }
     }
     if (!missing(cellnote))
          text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote),
               col = notecol, cex = notecex)
     par(mar = c(margins[1], 0, 0, 0))
     if (dendrogram %in% c("both", "row")) {
          plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
     }
     else plot.new()
     par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
     if (dendrogram %in% c("both", "column")) {
          plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
     }
     else plot.new()
     if (!is.null(main))
          title(main, cex.main = 1.5 * op[["cex.main"]])
     if (key) {
          #    par(mar = c(5, 4, 2, 1), cex = 0.75)
          par(mar = c(5, 4, 2, 1), cex = 1)
          tmpbreaks <- breaks
          if (symkey) {
               max.raw <- max(abs(c(x, breaks)), na.rm = TRUE)
               min.raw <- -max.raw
               tmpbreaks[1] <- -max(abs(x), na.rm = TRUE)
               tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm = TRUE)
          }
          else {
               #      min.raw <- min(x, na.rm = TRUE)
               #      max.raw <- max(x, na.rm = TRUE)
               if(length(breaks)>10){ # if custom breaks, then set to minimum to minimum of breaks, etc
                    min.raw <- min(breaks, na.rm = TRUE)
                    max.raw <- max(breaks, na.rm=TRUE)
                    print("min/max.raw1 hit")
               }else{ # else, set min, max to min/max of values
                    min.raw <- min(x, na.rm = TRUE)
                    max.raw <- max(x, na.rm = TRUE)
                    print("min/max.raw2 hit")
               }
               
          }
          
          z <- seq(min.raw, max.raw, length = length(col))
          image(z = matrix(z, ncol = 1), col = col, breaks = tmpbreaks,
                xaxt = "n", yaxt = "n")
          par(usr = c(0, 1, 0, 1))
          lv <- pretty(breaks)
          xv <- scale01(as.numeric(lv), min.raw, max.raw)
          axis(1, at = xv, labels = lv)
          if (scale == "row")
               mtext(side = 1, "Row Z-Score", line = 2)
          else if (scale == "column")
               mtext(side = 1, "Column Z-Score", line = 2)
          else mtext(side = 1, KeyValueName, line = 2)
          if (density.info == "density") {
               dens <- density(x, adjust = densadj, na.rm = TRUE)
               omit <- dens$x < min(breaks) | dens$x > max(breaks)
               dens$x <- dens$x[-omit]
               dens$y <- dens$y[-omit]
               dens$x <- scale01(dens$x, min.raw, max.raw)
               lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol,
                     lwd = 1)
               axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
               title("Color Key\nand Density Plot")
               #       par(cex = 0.5)
               par(cex = 1)
               mtext(side = 2, "Density", line = 2)
          }
          else if (density.info == "histogram") {
               h <- hist(x, plot = FALSE, breaks = breaks)
               hx <- scale01(breaks, min.raw, max.raw)
               hy <- c(h$counts, h$counts[length(h$counts)])
               lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s",
                     col = denscol)
               axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
               title("Color Key\nand Histogram")
               par(cex = 0.5)
               mtext(side = 2, "Count", line = 2)
          }
          else title("Color Key")
     }
     else plot.new()
     retval$colorTable <- data.frame(low = retval$breaks[-length(retval$breaks)],
                                     high = retval$breaks[-1], color = retval$col)
     invisible(retval)
}

generate_combined_avg_plots <- function(sample.name, name.of.genelist, bin, outputDirectory){
     directory = outputDirectory
     directory
     files <- list.files(directory, pattern=".*R1.*.txt")
     files <- files[c(1,3,2,4,6,5)]
     files
     # read in average values into list
     datalist <- list()
     d<- data.matrix(read.delim(paste0(directory,files[1]), header=TRUE))
     class(d)
     for(i in 1:length(files)){
          tmp <- data.matrix(read.table(paste0(directory,files[i]), header=TRUE))
          datalist[[i]]  <- tmp
          colnames(datalist[[i]]) <- gsub("_R1.*", "", colnames(datalist[[i]]))
     }
     class(datalist[[5]])
     datalist
     int <- as.numeric(unlist(read.delim(paste0(dir, "Michelle/Required_Files/Annotations/chipseq.heatmap/int.txt"), sep="\t", header=FALSE)))
     length(int)
     int

     # Average plot
     filename=paste0(directory, "combined_avg_plots_", name.of.genelist, "_", sample.name, ".jpeg")
     filename
     jpeg(filename, height=600, width=900, quality=100)
     tmpmat <- rbind(c(0,1,1,1,2))
     layout(tmpmat, widths=c(.1,.1))
     #PLOT1
     par(mar=c(5,7,5,0))
     plot(0~0, type="n", xlim=range(int), ylim=c(0,range(datalist, na.rm=T)[2]), xlab="Dist. from TSS", ylab="Seq Counts", cex.lab=1.8, cex.axis=1.8, cex.main=2, main=paste0("Combined Average Plot for ", sample.name, " timepoints (", name.of.genelist, ")"))
     plot.col <- vector()
     plot.lty <- vector()
     sample.name <- vector()
     custom.color <- c("gray0", "gray45", "gray80", "goldenrod4", "darkorange", "goldenrod2")
     for(u in 1:length(datalist)){
          sample.name <- c(sample.name, colnames(datalist[[u]]))
          # Start loop average plot
          plot.col[u] <- u
          if(!any(datalist[[u]] %in% NA)){
#                lines(int, datalist[[u]], lty=1, lwd=4, col=u)
               lines(int, datalist[[u]], lty=1, lwd=4, col=custom.color[u])
               plot.lty[u] <- 1
#           } else{lines(int, rep(0, num.bins), lty=3, lwd=4, col=u)
          } else{lines(int, rep(0, num.bins), lty=3, lwd=4, col=custom.color[u])
                 text(x=int[1], y=0, labels=c("coverage data for atleast two genes not present"), adj=0, cex=0.8)
                 plot.lty[u] <- 3}
     } # end plotting average plot
     #PLOT2 (dummy plot)
     barplot(t(c(1,2)), main="", col=NA, border="NA", axes=FALSE, names.arg=rep('',2), xpd=TRUE)
#      legend("left", fill=c(1:length(datalist)), legend=sample.name, border="white", box.col="white", cex=2)
     legend("left", fill=custom.color, legend=sample.name, border="white", box.col="white", cex=2)
     dev.off()
}
