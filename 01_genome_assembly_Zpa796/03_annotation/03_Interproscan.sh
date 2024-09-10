#!/bin/bash
 
#SBATCH --job-name=interpro #Give your job a name.
#SBATCH --nodes=1#Only increase for openmpi jobs.
#SBATCH --ntasks=1 #e.g if set to 2, could run two softwares in the script at the same time.
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4 #Multithreading.
#SBATCH --time=96:00:00 #Time for a job to run given as hh:mm:ss.
#SBATCH --mem=32G #Total Memory per node to use for the job
#sbatch --gridOptions="--time=24:00:00 --partition norm,b1"
#SBATCH --error=job.%J.err #Std Error write standard error to this file
#SBATCH --output=job.%J.out #Std Out write standard output to this file
#SBATCH --mail-type=FAIL #Notify user by email when certain event types occur (BEGIN, END, FAIL, REQUEUE)
#SBATCH --mail-user=rojas@evolbio.mpg.de #Email for notifications from previous line
#SBATCH --partition=global #Request a specific partition for the resource allocation.


################################
# Danilo Pereira - Kiel - 2023 #
################################
#
#
#
####################################################################################
# This script will annotate functions to the various protein sequences             #
#                                                                                  #
# Files needed:                                                                    #
#             * Target species protein fasta                                       #
#                                                                                  #
# Tools manuals: https://interproscan-docs.readthedocs.io/en/latest/HowToRun.html  #
#                                                                                  #
####################################################################################
#

##............................................
##      PROJECT ORGANIZATION
## Define your variables
##............................................

PROJECT=Zpa796.bd2
WORKDIR=/home/rojas/Zpasserinii/1_PacBio_assembly/6_braker_RNA/BRAKER.bbduk2/interproscan
PROTEINS=/home/rojas/Zpasserinii/1_PacBio_assembly/6_braker_RNA/BRAKER.bbduk2/braker/augustus.hints.aa

######################################################################################################################
#                                                    Start                                                           #
######################################################################################################################

################
# interproscan #
################
cd ${WORKDIR}

## copy data
cp ${PROTEINS} .
mv *.aa ${PROJECT}.fasta

sed -i 's/*//g' ${PROJECT}.fasta

# If protein sequences end in "*", it should removed
#sed -i 's/*//g' /home/pereira/2020_POSTDOC_MAIN/$Species_is/3_analysis/11_interproscan/Asp_NRRL3357_protein.simpleHeader.fasta

## run
module load java/x64/11u1q
module load python/3.7.1
echo -e '#!/bin/bash' > interpro.run.sh
echo "/data/biosoftware/InterProScan/interproscan-5.48-83.0/interproscan.sh -appl SignalP_EUK-4.1 -appl TMHMM-2.0c -dp --goterms --pathways -iprlookup -f TSV, GFF3 -i ${PROJECT}.fasta" >> interpro.run.sh

# wallace
sbatch --job-name=${PROJECT} --nodes=1 --ntasks=1 --ntasks-per-node=1 --cpus-per-task=1 --time=48:00:00 --mem=64G --error=job.%J.err --output=job.%J.out --mail-type=ALL --mail-user=rojas@evolbio.mpg.de --partition=standard < interpro.run.sh

######################################################################################################################
#                                                      END                                                           #
######################################################################################################################

