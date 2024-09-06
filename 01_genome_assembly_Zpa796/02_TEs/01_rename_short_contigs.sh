#!/bin/bash
 
#SBATCH --job-name=editAssm #Give your job a name.
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


####################################################################
#
#	PREPROCESING - REPET3 PIPELINE
# 
# 
#
####################################################################


#............................................
#PROJECT ORGANIZATION
#Define your variables
#............................................
WorkDir="/home/rojas/Zpasserinii/1_PacBio_assembly/5_Repet/REPET3.IPnp" #Path to the directory that contain the pipeline 
assemblyPath=" /home/rojas/Zpasserinii/1_PacBio_assembly/1_assembly/IPA/noPhase/assembly-results/final.p_ctg.fasta" # Path to the assembly. Input should be in fasta format
Project="Z796IPnp" #The ID of your project
ncontigs="31" # Number of contigs for your assembly, you can obtain the number of contigs with the comands provided below

#
#Note
#To obtain the number of contigs, you can use the following command
#grep '>' $grep '>' $assemblyPath | wc -l


#............................................
#Module 1.  Activate conda environment
#............................................

#Use seqkit from anaconda
#Activate anacon out of the script using the command 
 eval "$(/data/modules/python/python-miniconda3-qiime/bin/conda shell.bash hook)"

#Create an enviroment for the repet detection
#  $ conda create -n repet
#  $ conda install -n repet -c bioconda seqkit


# Activate environment

  $ conda activate repet


# To deactivate environment 
  $ conda deactivate




#............................................
#Module 2.  Rename contigs of my assembly
#............................................


#Rename assembly
cd $WorkDir
mkdir assembly
cd assembly
cp $assemblyPath .
mv *.fasta $Project.fna 

# Simplify scaffods header as suggested in BRAKER2 gitpagae
grep '>' $Project.fna > scaffolds.name.$Project.txt
sed -i.bak 's/>//' scaffolds.name.$Project.txt
cat scaffolds.name.$Project.txt | wc -l # gives the number of scafolds

#Change scaffolds name
for num in $( seq 1 $ncontigs )
do
Newheader=''$Project'_'$num''
echo $Newheader >> scaffolds.numbers.$Project.txt
done

# merge names and new header.
paste scaffolds.name.$Project.txt scaffolds.numbers.$Project.txt > scaffolds.name.numbers.$Project.txt

# replace with sed in the loop
while read line
do
linearray=( $line )
OldName=${linearray[0]}
NewName=${linearray[1]}
sed -i 's|'$OldName'|'$NewName'|g' $Project.fna 
done < scaffolds.name.numbers.$Project.txt 

echo "Renaming contigs finish"

rm scaffolds.name.$Project.txt 
rm scaffolds.numbers.$Project.txt
rm scaffolds.name.$Project.txt.bak

#............................................
#Module 2.  Short fasta file: each nucleic line has only 60 bps  (or less).
#............................................


seqkit seq $Project.fna  -w 60  > tmp.fna
mv tmp.fna $Project.fna

