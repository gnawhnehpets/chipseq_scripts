#!/bin/bash
# Normalizing BED

#FILES=(/home/steve/.gvfs/onc-analysis\$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/michelle/BED_files/*)
#for f in "${FILES[@]}"; do echo "$f"; done 

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
#    awk '{print($1"\t"$2"\t"$3"\t"int($4 * "${SCALINGFACTOR}"))}' $f > $f"_normalized.bed"
     awk '{print($1"\t"$2"\t"$3"\t"int($4 * "${SCALINGFACTOR}"))}' ../unnormalizedBED_10bp_bin/"${f/*bin\//}"_Coverage_TSS_10bp.bed > ${f/*bin\//}"_10bp_normalized_bash.BED" 
    
done