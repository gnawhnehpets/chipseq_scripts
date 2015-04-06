#!/bin/bash

##########
# Purpose:
##########
# Normalizing BED

# take files from: /home/steve/.gvfs/onc-analysis\$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/michelle/BED_files/Coverage_TSS_0bp_bin/unnormalizedBED_0bp_bin/)
# These contain all reads for each sample --< each line = 1 read
# wc -l  C10DDNMT1_R1.bam.sorted.bam_FILTERED.bed   #the bed file converted from bam. Each line is a read. so the number of total lines will be the total number of raw reads. Call this totalRawReads
# 10^8 / totalRawReads   #figure out the scaling factor. Call this scalingFactor
# Create new bed file with normalized values (from the bed file that is summarized into coverage bins)
# awk '{ printf ("%s\t%s\t%s\t%s\n", $1, $2, $3, $4* scalingFactor);}' C10DDNMT1_R1.bam.sorted.bam_FILTERED.bed_Coverage_200bp.bed  > C10DDNMT1_R1.bam.sorted.bam_FILTERED.bed_Coverage_200bp_normalized.bed



#FILES=*bed
FILES="../../Coverage_TSS_0bp_bin/unnormalizedBED_0bp_bin/*bed"
echo ${FILES}

for f in ${FILES}
do
    echo "FILE: $f"
    TOTALRAWREADS=$(wc -l < $f | grep -o "[0-9]\+")
    echo "total_raw_reads: ${TOTALRAWREADS}"
    SCALINGFACTOR=`awk 'BEGIN{printf("%0.4f", 100000000 / '$TOTALRAWREADS')}'`
    echo "scaling_factor:  ${SCALINGFACTOR}"
    echo "INPUT: ../unnormalizedBED_10bp_bin/"${f/*bin\//}"_Coverage_10bp.bed"
    echo "OUTPUT: "${f/*bin\//}"_Coverage_10bp.bed_normalized_bash.BED"
    awk '{print($1"\t"$2"\t"$3"\t"int($4 * "${SCALINGFACTOR}"))}' ../unnormalizedBED_10bp_bin/"${f/*bin\//}"_Coverage_TSS_10bp.bed > ${f/*bin\//}"_10bp_normalized_bash.BED" 
done