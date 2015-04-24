#!/bin/bash

##########
# Purpose:
##########

# After creating ChIP-seq plots, this script will consolidate the files by copying all files per sample type into a separate folder; copies are made so the originals are left alone.

########
# usage:
########

# Changes should not be necessary, but make sure the folder name correspond to the sample type

#FILES=(/home/steve/.gvfs/onc-analysis\$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/michelle/BED_files/*)
#for f in "${FILES[@]}"; do echo "$f"; done 

#FILES=*bed
#FILES=$(grep -rI "TEXTSEARCH") #find all text files with 'textsearch' in file

mkdir -p "DNMT1"
mkdir -p "EZH2"
mkdir -p "H3"
mkdir -p "H3K4"
mkdir -p "H3K27"
mkdir -p "Input"

FILES=$(find . -name "*.txt")
echo ${FILES}"\n"

find . -name "*DNMT1*raw*.jpeg" -exec cp {} ./DNMT1/ \;
find . -name "*DNMT1*.txt" -exec cp {} ./DNMT1/ \;

find . -name "*EZH2*raw*.jpeg" -exec cp {} ./EZH2/ \;
find . -name "*EZH2*.txt" -exec cp {} ./EZH2/ \;

find . -name "*H3_*raw*.jpeg" -exec cp {} ./H3/ \;
find . -name "*H3_*.txt" -exec cp {} ./H3/ \;

find . -name "*H3K4*raw*.jpeg" -exec cp {} ./H3K4/ \;
find . -name "*H3K4*.txt" -exec cp {} ./H3K4/ \;

find . -name "*H3K27*raw*.jpeg" -exec cp {} ./H3K27/ \;
find . -name "*H3K27*.txt" -exec cp {} ./H3K27/ \;

find . -name "*Input*raw*.jpeg" -exec cp {} ./Input/ \;
find . -name "*Input*.txt" -exec cp {} ./Input/ \;