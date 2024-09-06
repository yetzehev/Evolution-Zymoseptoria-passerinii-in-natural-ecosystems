#!/bin/bash

#SBATCH --job-name=mask
#SBATCH --ntasks=16
#SBATCH --nodes=1
#SBATCH --time=336:00:00
#SBATCH --mem=64G
#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out
#SBATCH --partition=global
#SBATCH --mail-type=FAIL #Notify user by email when certain event types occur (BEGIN, END, FAIL, REQUEUE)
#SBATCH --mail-user=rojas@evolbio.mpg.de #Email for notifications from previous line
#SBATCH --partition=global #Request a specific partition for the resource allocation.




#.................................................
#
#	Organize project directory
#
#
#.................................................
WorkDir="/home/rojas/Zpasserinii/1_PacBio_assembly/5_Repet/REPET3.IPnp" #Path to the directory that contain the pipeline 
Project="Z796IPnp" #The ID of your project
ncontigs="31" # Number of contigs for your assembly, you can obtain the number of contigs with the comands provided below


GENOMESFOLDER=/home/rojas/Zpasserinii/1_PacBio_assembly/5_Repet/REPET3.IPnp/assembly
TE_PREDICTION=/home/rojas/Zpasserinii/1_PacBio_assembly/5_Repet/REPET3.IPnp

cd $WorkDir
bedtools maskfasta -soft -fi $GENOMESFOLDER/$Project.fna -bed $TE_PREDICTION/$Project/FLFannot/${Project}_refTEs.gff -fo $GENOMESFOLDER/$Project.softmask.fasta
 
bedtools maskfasta -fi $GENOMESFOLDER/$Project.fna -bed $TE_PREDICTION/$Project/FLFannot/${Project}_refTEs.gff -fo $GENOMESFOLDER/$Project.hardmask.fasta
cd $WorkDir
