#!/bin/bash
 
#SBATCH --job-name=spades #Give your job a name.
#SBATCH --nodes=1#Only increase for openmpi jobs.
#SBATCH --ntasks=1 #e.g if set to 2, could run two softwares in the script at the same time.
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4 #Multithreading.
#SBATCH --time=24:00:00 #Time for a job to run given as hh:mm:ss.
#SBATCH --mem=32G #Total Memory per node to use for the job
#sbatch --gridOptions="--time=24:00:00 --partition norm,b1"
#SBATCH --error=job.%J.err #Std Error write standard error to this file
#SBATCH --output=job.%J.out #Std Out write standard output to this file
#SBATCH --mail-type=FAIL #Notify user by email when certain event types occur (BEGIN, END, FAIL, REQUEUE)
#SBATCH --mail-user=rojas@evolbio.mpg.de #Email for notifications from previous LINE
#SBATCH --partition=global #Request a specific partition for the resource allocation.


####################################################################
#
#	PREPROCESING - SPADES PIPELINE
# 
# 
#
####################################################################


#............................................
#PROJECT ORGANIZATION
#Define your variables
#............................................
WORKDIR=/home/rojas/septoriasis_model/02_comp_genomics/01_assemblies
BIN=/home/rojas/septoriasis_model/02_comp_genomics/01_assemblies/bin
ILLUMINAREADS=/home/rojas/septoriasis_model/00_SNVcalling/trimmed #Illumina reads directory
PROJECT="SSLB"
n=2 #cores

## Create working directories
cd $WORKDIR
mkdir rawfastq
mkdir quast
mkdir meta
mkdir out


## Create a soft link for the fastq files
cd rawfastq
ln -s $ILLUMINAREADS/*.fastq.gz .

#Make a list for the samples
cd $WORKDIR
ls $ILLUMINAREADS |  grep fastq | grep -v Zpa | grep -v fastqc | egrep 'R1|R2'|sed s/.bbduk_R1.fastq.gz// |sed s/.bbduk_R2.fastq.gz// | sort -r |uniq > ./meta/${PROJECT}_genomes_Zspp.txt

############################################################
##	Create a script for calling 
##	aligning and sorting each sample
############################################################
##Set array values
cd ${WORKDIR}
NAME=($(cat ./meta/${PROJECT}_genomes_Zspp.txt))
NJOB=($(seq 401 468)) ## As many paired samples: 112

## Open loop for creating individual scripts
cd ${BIN}
for i in "${!NJOB[@]}"; do


## Prepare your environment
## start assembly

echo -e '#!/bin/bash' >  ${NJOB[$i]}.spades.sh

echo "cd ${WORKDIR}" >>  ${NJOB[$i]}.spades.sh
echo "module load python/3.7.1" >>  ${NJOB[$i]}.spades.sh
echo "mkdir ${WORKDIR}/out/${NAME[$i]}" >>  ${NJOB[$i]}.spades.sh

echo "/data/biosoftware/spades/spades/spades.py -t 4 --careful -k 21,33,55,67,99,127 -o ${WORKDIR}/out/${NAME[$i]} --pe1-1 ${WORKDIR}/rawfastq/${NAME[$i]}.bbduk_R1.fastq.gz --pe1-2 ${WORKDIR}/rawfastq/${NAME[$i]}.bbduk_R2.fastq.gz" >>  ${NJOB[$i]}.spades.sh

echo "/data/biosoftware/quast/quast-4.6.3/quast.py ${WORKDIR}/out/${NAME[$i]}/scaffolds.fasta -t $n -o ${WORKDIR}/quast/${NAME[$i]}" >>  ${NJOB[$i]}.spades.sh

done

##############################################
## Submit an array
##############################################
#Make a list of jobs
ls *.spades.sh > list.jobs.txt

#Create submitting jobs
# Create the submission script
echo -e '#!/bin/bash' > array_batch.sh
echo -e '# get sample' >> array_batch.sh
echo -e 'printf -v SLURM_ARRAY_TASK_ID_three "%03d" $SLURM_ARRAY_TASK_ID' >> array_batch.sh # $SLURM_ARRAY_TASK_ID is a variable that takes a number among the array interval the user specify. But I need to make it 3-digits to use grep and avoid multiple matches, so this step.
echo -e 'sample_ID=$(cat list.jobs.txt | grep "$SLURM_ARRAY_TASK_ID.spades.sh")' >> array_batch.sh    #CHANGE NAME OF THE JOB
echo -e '# run the command' >> array_batch.sh
echo -e 'bash $sample_ID' >> array_batch.sh


#Submit job
sbatch --job-name=spades --array=401-468 --nodes=1 --ntasks=1 --ntasks-per-node=1 --cpus-per-task=1 --time=96:00:00 --mem=36G --error=job_call_%A_%a.err --output=job_call_%A_%a.out --mail-type=FAIL --mail-user=rojas@evolbio.mpg.de --partition=global array_batch.sh

