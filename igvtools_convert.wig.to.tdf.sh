#!/bin/bash;
#######
# usage
#######
# usage: bash convert.wig.to.tdf.sh
# check default shell (Ubuntu points to dash)
# > readlink -f $(which sh)
# run it with ./convert.wig.to.tdf.sh or bash convert.wig.to.tdf.sh
# running with sh ***.sh will not work because hashbang line will be
# ignored and script will be interpreted by dash (note: in ubuntu)

#######################
# igvtools toTDF usage:
#######################
# igvtools toTDF [options] [inputFile] [outputFile] [genome]

##########
# options:
##########
# -z num = specifies the maximum level to precompute. The default value is 7 and is sufficient for most files. To reduce file size at the expense of IGV performance, this value can be reduced
# -f list = a common delimited list specifying window functions to use when reducing the dat to precomputed tiels. Possible values are min, max, and mean. By default only the mean is calculated.
# -p file = specifies a 'bed' file to be used to map probe identifiers to locations. This option useful when preprocessing .gct files. The bed file should contain 4 columns:
#            chr  start  stop   name
#           where the name is the probe name in the .gct file

##########
# example:
##########
# /home/steve/Desktop/IGVTools/igvtools toTDF -z 5 C10DDNMT1_R1.bam.sorted.bam_FILTERED.bed_normalized.bed_output.wig C10DDNMT1_R1.bam.sorted.bam_FILTERED.bed_normalized.bed_output.tdf hg19

for input in ./*wig
do
    echo "Input: $input"
    newextension=tdf
    output="${input/wig/$newextension}"
    echo "Output: $output"
     /home/steve/Desktop/IGVTools/igvtools toTDF -z 5 $input "./iteration3/"$output hg19

done
