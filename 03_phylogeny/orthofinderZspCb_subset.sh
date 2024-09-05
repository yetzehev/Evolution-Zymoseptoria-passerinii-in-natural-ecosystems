#!/bin/bash
 
#SBATCH --job-name=orthof #Give your job a name.
#SBATCH --nodes=1#Only increase for openmpi jobs.
#SBATCH --ntasks=1 #e.g if set to 2, could run two softwares in the script at the same time.
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4 #Multithreading.
#SBATCH --time=96:00:00 #Time for a job to run given as hh:mm:ss.
#SBATCH --mem=48G #Total Memory per node to use for the job
#sbatch --gridOptions="--time=24:00:00 --partition norm,b1"
#SBATCH --error=job.%J.err #Std Error write standard error to this file
#SBATCH --output=job.%J.out #Std Out write standard output to this file
#SBATCH --mail-type=FAIL #Notify user by email when certain event types occur (BEGIN, END, FAIL, REQUEUE)
#SBATCH --mail-user=rojas@evolbio.mpg.de #Email for notifications from previous line
#SBATCH --partition=global #Request a specific partition for the resource allocation.


####################################################################
#
#	
#  ORTHOFINDER SUBSET
# 
#
####################################################################


##............................................
##	PROJECT ORGANIZATION
## Define your variables
##............................................
## Variables are defined in capital letters
PROJECT=Zsp_subset
WORKDIR=/home/rojas/septoriasis_model/02_comp_genomics/03_phylogeny/05_orthologous
INPUTS=/home/rojas/septoriasis_model/02_comp_genomics/03_phylogeny/03_annotationDB/proteomes_phylogeny
SUBSET=${WORKDIR}/meta/subset_2.phylogeny_IDs.txt


## Note activate conda environment before running the script
#source /data/modules/python/python-anaconda3/etc/profile.d/conda.sh
#conda activate orthofinder-2.5.5
#orthofinder

## Set environment
cd ${WORKDIR}

mkdir ${PROJECT}
PROTEOME=${WORKDIR}/${PROJECT}

##.........................................................
## Make a soft link for the proteome inferres with augustus
##.........................................................
cd ${PROTEOME}

while read ISOLATE;do
## C. beticola

cp ${INPUTS}/${ISOLATE}.aa.fasta .
##.........................................................
## Add species name to he header
##.........................................................

sed 's%^>\(.*\)%>'${ISOLATE}'|\1%I' ${ISOLATE}.aa.fasta > ${ISOLATE}.augustus.rehead.fasta
rm ${ISOLATE}.aa.fasta

done<${SUBSET}
 
 
 
 
cd ${WORKDIR}
orthofinder -f ${PROTEOME}
