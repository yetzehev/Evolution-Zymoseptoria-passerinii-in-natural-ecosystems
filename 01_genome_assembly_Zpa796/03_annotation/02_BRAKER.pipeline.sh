#!/bin/bash
 
#SBATCH --job-name=BRAKERbt3 #Give your job a name.
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
#	Braker pipeline Danilo Pereira
# RUN BRAKER2 - with flags from Stauber et al  
# The libraries and dependencies are installed on Wallace at  19.04.23
#
####################################################################
# braker2 version 2.1.6
# pipeline used BRAKER with proteins of any evolutionary distance
# what files from the output to use? https://github.com/Gaius-Augustus/BRAKER/issues/194

# SETUP
# OBS: I'm using an installation made by Kristian on wallace.
# Follow instructions in /data/biosoftware/braker2/READMEv2.1.6.txt 

# if white spaces on protein fasta, a warning is issued, but it can continue
#sed 's/ .*//' Sn15.2021_protein.fa > Sn15.2021.nowhite_protein.fa


## I M P O R T A N T
#1) Every time that you run BRAKER, the project name should be different or the files at the path should be deleted for reusing a project name
#ERROR MESSAGE
#data/biosoftware/braker2/deps/Augustus/config/species/PROJECT_NAME already exists. Choose another species name, delete this directory or use the existing species with the option --useexisting. Be aware that existing parameters will then be overwritten during training.
#

# 2) If Gene Mark fails!
#
# GeneMark requires a valid hidden key file in your home directory (~/.gm_key). The file expires after 200 days. Please check whether you have a valid key file before reporting an issue about this. 
# The key can be downloaded here:
#http://exon.gatech.edu/genemark/license_download.cgi

#............................................
#PROJECT ORGANIZATION
#Define your variables
#............................................
WorkDir=/home/rojas/Zpasserinii/1_PacBio_assembly/6_braker_RNA/BRAKER.bbduk2 #Path to the directory that contain the pipeline
RNAseq_bam=/home/rojas/Zpasserinii/1_PacBio_assembly/6_braker_RNA/trimming/bbduk.test2/bam/Zpa796.allReads.bbduk.hisat2.bam #BamFile with reads aligned against reference genome to annotate
Project="ZPA796.bt3" #The ID of your project, bbdk for reads trimmed with bbduk default
assemblyPath=/home/rojas/Zpasserinii/1_PacBio_assembly/5_Repet/REPET3.IPnp/assembly/Z796IPnp.softmask.fasta # Path to the assembly. Input should be in fasta format
targetSp=/groups/envgenom/Alice_Feurtey/Scripts/Assemblies_and_annotations/2_Annotations/Orthofinder/Zpa63_EVM.all.shorter_names.prot.fasta #Path to the fasta file that contains the protein sequence for your target species
orthodb="https://v100.orthodb.org/download/odb10_fungi_fasta.tar.gz" #Enter the URL address for orthodb
ncontigs="31" # Number of contigs for your assembly, you can obtain the number of contigs with the comands provided below
#
#
## Note
## To obtain the number of contigs, you can use the following command
## grep '>' $grep '>' $assemblyPath | wc -l
##
##...........................................
## Module 1.  Rename contigs of my assembly
##............................................

## Rename assembly
cd $WorkDir
mkdir assembly
cd assembly
cp $assemblyPath .
mv *.fasta $Project.fasta 


## Remove white spaces at contigs nahELLO me
sed -e 's/ /_/g' $Project.fasta  > tmp.fasta
mv tmp.fasta $Project.fasta

## Simplify scaffods header as suggested in BRAKER2 gitpagae
grep '>' $Project.fasta > scaffolds.name.$Project.txt
sed -i.bak 's/>//' scaffolds.name.$Project.txt
cat scaffolds.name.$Project.txt | wc -l # gives the number of scafolds

## Change scaffolds name
for num in $( seq 1 $ncontigs )
do
Newheader='contig_'$num''
echo $Newheader >> scaffolds.numbers.$Project.txt
done

## merge names and new header.
paste scaffolds.name.$Project.txt scaffolds.numbers.$Project.txt > scaffolds.name.numbers.$Project.txt

## replace with sed in the loop
while read line
do
linearray=( $line )
OldName=${linearray[0]}
NewName=${linearray[1]}
sed -i 's|'$OldName'|'$NewName'|g' $Project.fasta 
done < scaffolds.name.numbers.$Project.txt 

echo "Renaming contigs finish"

rm scaffolds.name.$Project.txt 
rm scaffolds.numbers.$Project.txt
rm scaffolds.name.$Project.txt.bak

cd $WorkDir

##................................................
## Module 2.  Download and create your database
##................................................
mkdir targetSp.proteins
cp $targetSp ./targetSp.proteins
mkdir ./db; cd  ./db

## get Fungi protein database from orthodb. Needs to be done once. Then, for each outgroup, add the proteins from target into the "clean" database.
wget --no-check-certificate $orthodb
tar -xf odb10_fungi_fasta.tar.gz

## After download, follow how to proceed: https://github.com/gatech-genemark/ProtHint#protein-database-preparation
## Copy the protein fasta file from the target fungi, and cat it together with the concatanated fasta from previous step. Basically, add the target species proteins into the database.
cat ./fungi/Rawdata/*.fs $WorkDir/targetSp.proteins/*.fasta  > orthoDb.$Project.fasta

rm odb10_fungi_fasta.tar.gz
cd $WorkDir


##................................................
## Module 3.  RUN BRAKER
##................................................


## exports and module for wallace
## cp /data/biosoftware/braker2/deps/gm_key_64 $HOME/.gm_key
module load perl/5.26.1
source /data/modules/python/python-anaconda3/etc/profile.d/conda.sh
conda activate
module load java/x64/8u121
export PATH=/data/biosoftware/braker2/BRAKER-2.1.6/scripts:$PATH

## create script for wallace using BRAKER-2.1.4
n=8
AUGUSTUS_BASE=/data/biosoftware/braker2/deps/Augustus/

export AUGUSTUS_CONFIG_PATH=$AUGUSTUS_BASE/config/
export AUGUSTUS_SCRIPTS_PATH=$AUGUSTUS_BASE/scripts
export AUGUSTUS_BIN_PATH=$AUGUSTUS_BASE/bin
export GENEMARK_PATH=/data/biosoftware/braker2/deps/gmes_linux_64/
export BAMTOOLS_PATH=/data/biosoftware/braker2/deps/bamtools/bin/usr/local/bin/
export DIAMOND_PATH=/data/biosoftware/braker2/deps/diamond/
export BLAST_PATH=/data/biosoftware/braker2/deps/ncbi-blast-2.11.0+/bin/
export PROTHINT_PATH=/data/biosoftware/braker2/deps/ProtHint/bin/
export SAMTOOLS_PATH=/data/biosoftware/braker2/deps/samtools/
export CDBTOOLS_PATH=/data/biosoftware/braker2/deps/cdbfasta/


##................................................
## Module 3.1 OPTIONAL!!! If you have RNA data, 
## you can add the info for improve annotation of the genome
##................................................
## The RNA reads should be trimmend, filter by quality and  aligned to the  $REF_genome
## from Module 3
## IMPORTANT!
##
##
cd ${WorkDir}
cd assembly
bam2hints --in=${RNAseq_bam} --out=./${Project}_hints.gff
sed 's/Z796IPnp/contig/g' ${Project}_hints.gff > tmp  #Change contigs header
mv tmp  ${Project}_hints.gff

## The contigs are renamed during BREAKER pipeline 
## Check the contigs of the file ${Project}_hints.gff, if these do not match the pattern
## contig_XX, rename them
## sed 's/OLDContig_ID/contig/g' *_hints.gff > tmp > $Project_hints.gff

cd $WorkDir


## References
REF_genome="${WorkDir}/assembly/${Project}.fasta"
BAM_hints="${WorkDir}/assembly/${Project}_hints.gff"
PROTEIN_database="${WorkDir}/db/orthoDb.$Project.fasta"

echo -e '#!/bin/bash' > sbatch_BRAKER2.sh
echo "braker.pl --fungus --cores $n --etpmode \
--species=$Project \
--genome=$REF_genome \
--hints=$BAM_hints \
--prot_seq=$PROTEIN_database \
--UTR=off --gff3 --alternatives-from-evidence=false" >> sbatch_BRAKER2.sh

sbatch --job-name=BRAKERbt2 --nodes=1 --ntasks=1 --ntasks-per-node=1 --cpus-per-task=8 --time=72:00:00 --mem=64G --error=job.%J.err --output=job.%J.out --mail-type=ALL --mail-user=rojas@evolbio.mpg.de --partition=highmem < sbatch_BRAKER2.sh
