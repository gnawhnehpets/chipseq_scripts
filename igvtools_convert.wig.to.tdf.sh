#!/bin/bash;

# usage: bash convert.wig.to.tdf.sh
# check default shell (Ubuntu points to dash)
# > readlink -f $(which sh)
# run it with ./convert.wig.to.tdf.sh or bash convert.wig.to.tdf.sh
# running with sh ***.sh will not work because hashbang line will be
# ignored and script will be interpreted by dash (note: in ubuntu)
for input in ./*wig
do
    echo "Input: $input"
    newextension=tdf
    output="${input/wig/$newextension}"
    echo "Output: $output"
     /home/steve/Desktop/IGVTools/igvtools toTDF -z 5 $input "./iteration3/"$output hg19
    
#    /home/steve/Desktop/IGVTools/igvtools toTDF -z 5 C10DDNMT1_R1.bam.sorted.bam_FILTERED.bed_normalized.bed_output.wig C10DDNMT1_R1.bam.sorted.bam_FILTERED.bed_normalized.bed_output.tdf hg19
#    /home/steve/Desktop/IGVTools/igvtools toTDF -z 5 $file C10DDNMT1_R1.bam.sorted.bam_FILTERED.bed_normalized.bed_output.tdf hg19
    
#WORKS
#    echo ${input/wig/tdf}
#    echo "New filename: " ${input/.wig/.tdf}
done

#FILES="/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/Michelle/BED_files/Coverage_TSS_200bp_bin/bedgraph_to_wig_200bp_bin/*"
#for f in $files
#do 
#    echo "File: $f"
#done