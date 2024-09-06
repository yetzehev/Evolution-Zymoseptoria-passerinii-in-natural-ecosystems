#!/bin/bash
 
#SBATCH --job-name=editAssm #Give your job a name.
#SBATCH --nodes=1#Only increase for openmpi jobs.
#SBATCH --ntasks=1 #e.g if set to 2, could run two softwares in the script at the same time.
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4 #Multithreading.
#SBATCH --time=24:00:00 #Time for a job to run given as hh:mm:ss.
#SBATCH --mem=12G #Total Memory per node to use for the job
#sbatch --gridOptions="--time=24:00:00 --partition norm,b1"
#SBATCH --error=job.%J.err #Std Error write standard error to this file
#SBATCH --output=job.%J.out #Std Out write standard output to this file
#SBATCH --mail-type=FAIL #Notify user by email when certain event types occur (BEGIN, END, FAIL, REQUEUE)
#SBATCH --mail-user=rojas@evolbio.mpg.de #Email for notifications from previous line
#SBATCH --partition=fast #Request a specific partition for the resource allocation.


####################################################################
#
#       PREPROCESING - BRAKER (DB + RNA)
# 
# 
#
####################################################################


#............................................
#PROJECT ORGANIZATION
#Define your variables
#............................................
WORKDIR=/home/rojas/Zpasserinii/1_PacBio_assembly/6_braker_RNA/trimming
READS=/home/rojas/Zpasserinii/data/SSLB_RNA
PROJECT=ANOT.TRIM


#............................................
#SET YOUR ENVIRONMENT
# THIS SCRIPT INCLUDES MODIFICATIONS TO THE bbdul trimming. I added quality trimming
#............................................
cd $WORKDIR

#SOFTWARE
module load java/x64/8u121
export BBMAPDIR=/data/biosoftware/bbmap/bbmap


#Make a samples list
ls $READS |  grep fastq.gz | egrep 'R1|R2'|sed s/_R1_001.fastq.gz// |sed s/_R2_001.fastq.gz// | sort -r |uniq > $WORKDIR/meta/samplelist_$PROJECT.txt


#Set array values
TAG=($(cat $WORKDIR/meta/samplelist_$PROJECT.txt))
NJOB=($(seq 100 102)) ## It assigns a three digits index number per genotype.

for i in "${!NJOB[@]}"; do
# For trimming I'll use bbmap/bbduk
cd $WORKDIR
cd bin2

echo -e '#!/bin/bash' > ${NJOB[$i]}.trimming.sh
echo "$BBMAPDIR/bbduk.sh in1=$READS/${TAG[$i]}_R1_001.fastq.gz in2=$READS/${TAG[$i]}_R2_001.fastq.gz out1=$WORKDIR/out2/${TAG[$i]}_bbduk_1.fastq.gz out2=$WORKDIR/out2/${TAG[$i]}_bbduk_2.fastq.gz ref=$BBMAPDIR/resources/adapters.fa trimpolya=10 ktrim=r k=23 mink=11 hdist=1 qtrim=rl trimq=20 maq=20 tpe tbo &> $WORKDIR/out/${TAG[$i]}_${NJOB[$i]}.bbduk.log" >> ${NJOB[$i]}.trimming.sh

done


ls *trimming.sh > list.jobs.txt


#Create submitting jobs
# Create the submission script
echo -e '#!/bin/bash' > array_batch.sh
echo -e 'module load java/x64/8u121' >> array_batch.sh
echo -e '# get sample' >> array_batch.sh
echo -e 'printf -v SLURM_ARRAY_TASK_ID_three "%03d" $SLURM_ARRAY_TASK_ID' >> array_batch.sh # $SLURM_ARRAY_TASK_ID is a variable that takes a number among the array interval the user specify. But I need to make it 4-digits to use grep and avoid multiple matches, so this step.
echo -e 'sample_ID=$(cat list.jobs.txt | grep "$SLURM_ARRAY_TASK_ID.trimming.sh")' >> array_batch.sh
echo -e '# run the command' >> array_batch.sh
echo -e 'bash $sample_ID' >> array_batch.sh


sbatch --job-name=bbduk --array=100-102 --nodes=1 --ntasks=1 --ntasks-per-node=1 --cpus-per-task=1 --time=96:00:00 --mem=24G --error=job_%A_%a.err --output=job_name_%A_%a.out --mail-type=FAIL --mail-user=rojas@evolbio.mpg.de --partition=global array_batch.sh
