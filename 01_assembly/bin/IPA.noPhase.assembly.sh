#!/bin/bash
 
#SBATCH --job-name=IPAnP #Give your job a name.
#SBATCH --nodes=1#Only increase for openmpi jobs.
#SBATCH --ntasks=1 #e.g if set to 2, could run two softwares in the script at the same time.
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8 #Multithreading.
#SBATCH --time=96:00:00 #Time for a job to run given as hh:mm:ss.
#SBATCH --mem=40G #Total Memory per node to use for the job
#SBATCH --error=job.%J.err #Std Error write standard error to this file
#SBATCH --output=job.%J.out #Std Out write standard output to this file
#SBATCH --mail-type=FAIL #Notify user by email when certain event types occur (BEGIN, END, FAIL, REQUEUE)
#SBATCH --mail-user=rojas@evolbio.mpg.de #Email for notifications from previous line
#SBATCH --partition=standard #Request a specific partition for the resource allocation.

## Activate conda environment
#source /data/modules/python/python-anaconda3/etc/profile.d/conda.sh
#conda activate ipa


#Project organization
project="Zpa796_assembly"
home="/home/rojas/Zpasserinii/1_PacBio_assembly"
HiFi="/home/rojas/Zpasserinii/data/PacBio/HiFi/m64078e_221226_215108.hifi_reads.fastq.gz"
bin="$home/IPA"

ipa local --no-phase --nthreads 8 --njobs 1 --run-dir ./noPhase -i $HiFi --resume
