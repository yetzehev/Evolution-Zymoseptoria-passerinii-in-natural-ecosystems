#!/bin/bash
 
#SBATCH --job-name=combGVCF #Give your job a name. 
#SBATCH --nodes=1#Only increase for openmpi jobs.
#SBATCH --ntasks=1 #e.g if set to 2, could run two softwares in the script at the same time.
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4 #Multithreading.
#SBATCH --time=96:00:00 #Time for a job to run given as hh:mm:ss.
#SBATCH --mem=64G #Total Memory per node to use for the job
#SBATCH --error=job.%J.err #Std Error write standard error to this file
#SBATCH --output=job.%J.out #Std Out write standard output to this file
#SBATCH --mail-type=FAIL #Notify user by email when certain event types occur (BEGIN, END, FAIL, REQUEUE)
#SBATCH --mail-user=rojas@evolbio.mpg.de #Email for notifications from previous line
#SBATCH --partition=standard #Request a specific partition for the resource allocation.


####################################################################
#
#      Combine  GVCF
# 
####################################################################
##............................................
##      PROJECT ORGANIZATION
## Define your variables and substitute values 
##............................................

## Variables are defined in capital letters
PROJECT=SSLB_model
WORKDIR=/home/rojas/septoriasis_model/00_SNVcalling
META=${WORKDIR}/meta
BIN=${WORKDIR}/bin # This directory contains this script

## INPUTS
GVCF=${WORKDIR}/gvcf ## input GVCF files
REF=/home/rojas/ref_genomes/Zpa796/IPAnp/assembly/Z796IPnp.fasta ## reference genome

##OUTPUTS
mkdir ${WORKDIR}/out
OUT=${WORKDIR}/out

## Make a list with all the samples
ls ${GVCF} | grep "vcf" | grep -v ".idx" > ${META}/GVCFs.samplelist_${PROJECT}.txt
tmp=""

while read NAME 
do
        tmp="$tmp --variant ${GVCF}/${NAME}"
done < ${META}/GVCFs.samplelist_${PROJECT}.txt

#For GATK 4
#CombineGVCFs is meant to be used for merging of GVCFs that will eventually be input into GenotypeGVCFs.
gatk --java-options "-Xmx64g" CombineGVCFs \
   -R ${REF} \
   ${tmp} \
   -O ${OUT}/${PROJECT}_v4.2.5.comGty.vcf.gz \
   --read-filter MappingQualityReadFilter \
   --read-filter OverclippedReadFilter \
   
#GenotypeGVCF: perform joint genotyping on a single input

gatk --java-options "-Xmx64g" GenotypeGVCFs \
 -R ${REF} \
 -V $OUT/${PROJECT}_v4.2.5.comGty.vcf.gz  \
 -ploidy 1 \
 -O $OUT/${PROJECT}_v4.2.5.vcf.gz
 
 rm $OUT/${PROJECT}_v4.2.5.comGty.vcf.gz 
