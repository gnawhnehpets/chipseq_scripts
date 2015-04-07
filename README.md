\#scripts used for chipseq analysis

prerequisites - BED files normalized using scaling factor: (100000000/total read count)*read count per bin

\######################

\# Normalization method

\######################

wc -l original.bed.files # equal to # of total reads

TOTALRAWREADS=$(wc -l original.bed.file)
SCALINGFACTOR=100000000/${TOTALRAWREADS}

Multiply each binned values by the SCALINGFACTOR

See normalizeBED3.sh script for automation
