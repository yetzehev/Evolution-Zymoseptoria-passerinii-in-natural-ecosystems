#!/bin/bash
 
#SBATCH --job-name=busco #Give your job a name.
#SBATCH --nodes=1#Only increase for openmpi jobs.
#SBATCH --ntasks=1 #e.g if set to 2, could run two softwares in the script at the same time.
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=6 #Multithreading.
#SBATCH --time=96:00:00 #Time for a job to run given as hh:mm:ss.
#SBATCH --mem=64G #Total Memory per node to use for the job
#sbatch --gridOptions="--time=24:00:00 --partition norm,b1"
#SBATCH --error=job.%J.err #Std Error write standard error to this file
#SBATCH --output=job.%J.out #Std Out write standard output to this file
#SBATCH --mail-type=FAIL #Notify user by email when certain event types occur (BEGIN, END, FAIL, REQUEUE)
#SBATCH --mail-user=rojas@evolbio.mpg.de #Email for notifications from previous line
#SBATCH --partition=standard #global #Request a specific partition for the resource allocation.


####################################################################
#
#       
# 	Run BUSCO for Zpa796 IPA non-phased
# 
#
####################################################################


##............................................
##      PROJECT ORGANIZATION
## Define your variables
##............................................
## Variables are defined in capital letters

WORKDIR=/home/rojas/Zpasserinii/1_PacBio_assembly
BIN=/home/rojas/Zpasserinii/1_PacBio_assembly/1_assembly/bin
OUT=/home/rojas/Zpasserinii/1_PacBio_assembly/3_BUSCO_Z796
GENOME=/home/rojas/Zpasserinii/1_PacBio_assembly/1_assemblies_edited/IPAnp/assembly/Z796IPnp.fasta



##............................................
##       
## 	SOFTWARE
##............................................
## activate environment
# source /data/modules/python/python-anaconda3/etc/profile.d/conda.sh

cd ${WORKDIR}
# conda activate busco_local
#busco -i  /home/rojas/Zpasserinii/1_PacBio_assembly/1_assemblies_edited/IPAnp/assembly/Z796IPnp.fna -m genome -c 6 -o busco_run6 --offline -f -l /home/rojas/databases/busco/lineages/fungi_odb10 --out_path $OUT 
busco -i  /home/rojas/Zpasserinii/1_PacBio_assembly/1_assemblies_edited/IPAnp/assembly/Z796IPnp.fna -m genome -c 6 -o busco_asc --offline -f -l /home/rojas/databases/busco/lineages/ascomycota_odb10 --out_path $OUT
