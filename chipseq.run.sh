#!/usr/bin/bash

#######
#Usage:
#######
# place script into output_directory for genelist
# bash chipseq.run.sh

LOG="./log.txt"
rm -f $LOG

printf "##############################%s\n"
printf "# chipseq_script_commandprompt%s\n"
printf "##############################%s\n"

Rscript "/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/Michelle/Rscripts/chipseq_script_commandprompt.R" "${PWD##*/}" normalized 200 TRUE > $LOG

printf "######################%s\n"
printf "# consolidate_files.sh%s\n"
printf "######################%s\n"

bash "/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/Michelle/Rscripts/consolidate_files.sh" >> $LOG

printf "###################%s\n"
printf "# combine_avg_plots%s\n"
printf "###################%s\n"

Rscript "/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/Michelle/Rscripts/combine_average_plots_commandprompt.R" "${PWD##*/}" 200 >> $LOG

#######
#Usage:
#######

# bash consolidate_and_create.combined.avg.plots.sh genelist.name

#bash "/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/Michelle/Rscripts/consolidate_files.sh"

#Rscript "/home/steve/.gvfs/onc-analysis$ on onc-cbio2.win.ad.jhu.edu/users/shwang26/Michelle/Rscripts/combine_average_plots_commandprompt.R" $1 200