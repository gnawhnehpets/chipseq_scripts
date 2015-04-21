#!/usr/bin/bash

#######
#Usage:
#######
# place script into directory containing chipseq_script_commandprompt.R output and run
# bash consolidate_and_create.combined.avg.plots.sh

bash "/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/Michelle/Rscripts/consolidate_files.sh"

Rscript "/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/Michelle/Rscripts/combine_average_plots_commandprompt.R" "${PWD##*/}" 200

#######
#Usage:
#######

# bash consolidate_and_create.combined.avg.plots.sh genelist.name

#bash "/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/Michelle/Rscripts/consolidate_files.sh"

#Rscript "/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/Michelle/Rscripts/combine_average_plots_commandprompt.R" $1 200